function compute_train_set(FILES,varargin)
%Computes threshold for ISI and waveform correlations (JSD and
%normalized peak xcorr)
%
%


nparams=length(varargin);
bootstrap=0; % get additional estimates through bootstrapping trials
isibins=linspace(0,.05,200);

if mod(nparams,2)>0
	error('ephysPipeline:argChk','Parameters must be specified as parameter/value pairs!');
end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'labels'
			labels=varargin{i+1};
		case 'bootstrap'
			bootstrap=varargin{i+1};
	end
end

FILES.match.name={};
FILES.match.idx=[];
FILES.nomatch.name={};
FILES.match.clustidx=[];
FILES.nomatch.clustidx=[];

% collect two primary file lists, allow the user to multi-select at the end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DATA COLLECTION

counter=1;
if nargin<1 | isempty(FILES)

	% get the match directory

	flag=0;

	disp('Select root directory for MATCH examples');

	while ~flag

		[filelist clusterids]=recursive_file_find(pwd,'sua_channels','Choose directories FOR ONE CELL for MATCH');

		for j=1:length(filelist)
			FILES.match.name{end+1}=filelist{j};
			FILES.match.clustidx(end+1)=clusterids(j);
			FILES.match.idx(end+1)=counter;
		end

		response=[];

		while isempty(response)
			response=input('(D)one choosing root folders or (c)ontinue?  ','s');

			switch lower(response(1))


				case 'd'
					flag=1;
				case 'c'
					
				otherwise
					response=[];

			end
		end

		counter=counter+1;

	end

	% get the non-match directory

	selectfig=figure('Units','Pixels','Position',[0 0 1e3 600]);
	test=uicontrol(selectfig,'style','listbox','max',10,'min',0,...
		'Position',[20 20 900 500],'string',FILES.match.name);
	exitbutton=uicontrol(selectfig,'style','pushbutton',...
		'Position',[400 530 250 50],'string','DONE',...
		'callback','uiresume(gcbf)');
	disp('Click done after highlighting all MATCH files');
	uiwait(selectfig);

	selected_files=get(test,'value');
	close(selectfig);

	FILES.match.name=FILES.match.name(selected_files);
	FILES.match.idx=FILES.match.idx(selected_files);
	FILES.match.clustidx=FILES.match.clustidx(selected_files);

	flag=0;

	disp('Select root directory for NON-MATCH examples');

	while ~flag

		[filelist clusterids]=recursive_file_find(pwd,'sua_channels','Choose root directory for NON-MATCH');

		for j=1:length(filelist)
			FILES.nomatch.name{end+1}=filelist{j};
			FILES.nomatch.clustidx(end+1)=clusterids(j);
		end

		response=[];

		while isempty(response)
			response=input('(D)one choosing root folders or (c)ontinue?  ','s');

			switch lower(response(1))


				case 'd'
					flag=1;
				case 'c'
					
				otherwise
					response=[];

			end
		end

		counter=counter+1;

	end

	selectfig=figure('Units','Pixels','Position',[0 0 1e3 600]);
	test=uicontrol(selectfig,'style','listbox','max',10,'min',0,...
		'Position',[20 20 900 500],'string',FILES.nomatch.name);
	exitbutton=uicontrol(selectfig,'style','pushbutton',...
		'Position',[400 530 250 50],'string','DONE',...
		'callback','uiresume(gcbf)');
	disp('Click done after highlighting all NON-MATCH files');
	uiwait(selectfig);

	selected_files=get(test,'value');
	close(selectfig);

	FILES.nomatch.name=FILES.nomatch.name(selected_files);
	FILES.nomatch.clustidx=FILES.nomatch.clustidx(selected_files);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SCORE COMPUTATION

% assumes data is stored in two main directories, match and nomatch
% scores from match  are boostrapped to form the match distribution, 
% these are then compared with no match

% hopefully a region is occupied by the matches such that we can choose
% cutoff for JSD and xcorr score could use MAD rule-of-thumb

% now grab the relevant waveforms and ISIs

ncells=unique(FILES.match.idx);

% how many cells do we have for the matching condition

match_wavescore=[];
match_isiscore=[];

nomatch_wavescore=[];
nomatch_isiscore=[];

train_matrix=[];
class=[];


for i=1:length(ncells)

	% take the data from the first file, compare to all others,
	% add to the MATCH distribution

	idxs=find(FILES.match.idx==ncells(i));
	currclusters=FILES.match.clustidx(idxs);

	load(FILES.match.name{idxs(1)},'clust_spike_vec','clusterwindows','IFR');

	templatewave=mean(clusterwindows{currclusters(1)},2);
	templatewave=zscore(templatewave);
	templateifr=mean(IFR{1}{currclusters(1)});
	templateifr=zscore(templateifr);
	templateisi=get_isi(clust_spike_vec{1}{currclusters(1)});

	for j=2:length(idxs)

		% definitely need to bootstrap here, too few datapoints...

		[wavescore isiscore ifrscore]=get_scores(templatewave,templateisi,templateifr,...
			FILES.match.name{idxs(j)},currclusters(j),isibins,bootstrap);
		
		train_matrix=[train_matrix;[wavescore isiscore ifrscore]];

		class=[class;ones(size(wavescore))];

	end

	% include the match list, non-matching cells

	idxs=find(FILES.match.idx~=ncells(i));
	currclusters=FILES.match.clustidx(idxs);

	for j=1:length(idxs)

		[wavescore isiscore ifrscore]=get_scores(templatewave,templateisi,templateifr,...
			FILES.match.name{idxs(j)},currclusters(j),isibins,0);


		train_matrix=[train_matrix;[wavescore isiscore ifrscore]];
		class=[class;0];


	end


	for j=1:length(FILES.nomatch.name)

		[wavescore isiscore ifrscore]=get_scores(templatewave,templateisi,templateifr,...
			FILES.nomatch.name{(j)},FILES.nomatch.clustidx(j),isibins,0);

		train_matrix=[train_matrix;[wavescore isiscore ifrscore]];
		class=[class;0];


	end

end

% compute the nomatch scores, can now include match data (just don't compare
% within cells)

save('training_data.mat','train_matrix','class','FILES');

end

function [LIST CLUSTERID]=recursive_file_find(ROOTDIR,FILTER,DIALOG)

rootfolder=uigetdir(ROOTDIR,DIALOG);

subfolders_pre=genpath(rootfolder);
subfolders=regexp(subfolders_pre,pathsep,'split');

LIST={};
CLUSTERID=[];

for i=1:length(subfolders)

	filelisting_pre=dir(subfolders{i});
	filelisting={filelisting_pre(:).name};

	% read files in order, look for first file with cluster in the filename

	for j=1:length(filelisting)
		
		% also get the cluster number!

		if findstr(filelisting{j},'cluster') 
			[junk1 junk2 junk3 match]=regexp(filelisting{j},'cluster_([0-9]+)','match');
			CLUSTERID(end+1)=str2num(filelisting{j}(match{1}(1):match{1}(2)));
			break;
		end

	end


	for j=1:length(filelisting)
		
		if findstr(filelisting{j},FILTER)
			LIST{end+1}=fullfile(subfolders{i},filelisting{j});
			break;
		end

	end
end


end

function [WAVESCORE ISISCORE IFRSCORE]=...
		get_scores(TEMPLATEWAVE,TEMPLATEISI,TEMPLATEIFR,FILE,CLUSTER,ISIBINS,BOOTSTRAP)


load(FILE,'clust_spike_vec','clusterwindows','IFR');

if BOOTSTRAP>0

	ntrials=size(IFR{1}{CLUSTER},1);
	trialpool=1:ntrials;

	nspikes=size(clusterwindows{CLUSTER},2);

	parfor i=1:BOOTSTRAP

		trials=trialpool(randsample(ntrials,ntrials,1));

		tmpwave=mean(clusterwindows{CLUSTER}(:,trials),2);
		tmpwave=zscore(tmpwave);

		tmpifr=mean(IFR{1}{CLUSTER}(trials,:));
		tmpifr=zscore(tmpifr);

		[r,lags]=xcov(tmpwave,TEMPLATEWAVE,'coeff');

		WAVESCORE(i,1)=max(r);

		[r,lags]=xcov(tmpifr,TEMPLATEIFR,'none');
		r=r./(norm(TEMPLATEIFR)*norm(tmpifr));

		IFRSCORE(i,1)=max(r);

		tmpisi=get_isi(clust_spike_vec{1}{CLUSTER}(trials));	

		ISISCORE(i,1)=sqrt(kld(tmpisi,TEMPLATEISI,'jsd',1,'bins',ISIBINS));

	end

else

	scorewave=mean(clusterwindows{CLUSTER},2);
	scorewave=zscore(scorewave);
	scoreisi=get_isi(clust_spike_vec{1}{CLUSTER});
	scoreifr=mean(IFR{1}{CLUSTER});
	scoreifr=zscore(scoreifr);

	[r,lags]=xcov(scorewave,TEMPLATEWAVE,'coeff');

	WAVESCORE=max(r);

	[r,lags]=xcov(scoreifr,TEMPLATEIFR,'none');

	r=r./(norm(TEMPLATEIFR)*norm(scoreifr));
	IFRSCORE=max(r);

	ISISCORE=sqrt(kld(scoreisi,TEMPLATEISI,'jsd',1,'bins',ISIBINS));
end

end
