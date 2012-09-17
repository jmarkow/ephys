function compute_test_set(FILES,varargin)
%Computes threshold for ISI and waveform correlations (JSD and
%normalized peak xcorr)
%
%


nparams=length(varargin);
spike_fs=100e3;
isibins=linspace(0,.05,200);

if mod(nparams,2)>0
	error('ephysPipeline:argChk','Parameters must be specified as parameter/value pairs!');
end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'spike_fs'
			tail=varargin{i+1};
		case 'labels'
			labels=varargin{i+1};
	end
end

FILES.name={};
FILES.clustidx=[];
FILES.idx=[];
FILES.date=[];

% collect two primary file lists, allow the user to multi-select at the end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DATA COLLECTION

cellidx=0;
if nargin<1 | isempty(FILES)

	% get the match directory

	flag=0;

	while ~flag

		[filelist clusterids containers]=recursive_file_find(pwd,'sua_channels','Choose directories with cells to MATCH');

		% in this case we need to strip the last directory from the container

		for j=1:length(containers)
			tokens=regexp(containers{j},filesep,'split');

			newpath='';
			for k=1:length(tokens)-2
				newpath=[newpath tokens{k} filesep];
			end

			filedate{j}=datenum(regexp(tokens{end-1},'^\d+-\d+-\d+','match'),'yyyy-mm-dd');
			containers{j}=newpath;

		end

		% unique root directories
		% also collect the dates!

		uniq_dirs=unique(containers);

		for j=1:length(filelist)
			FILES.name{end+1}=filelist{j};
			FILES.clustidx(end+1)=clusterids(j);
			FILES.idx(end+1)=find(strcmp(containers{j},uniq_dirs))+cellidx;
			FILES.date(end+1)=filedate{j};
		end

		cellidx=cellidx+length(uniq_dirs);
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

	end

	selectfig=figure('Units','Pixels','Position',[0 0 1e3 600]);
	test=uicontrol(selectfig,'style','listbox','max',10,'min',0,...
		'Position',[20 20 900 500],'string',FILES.name);
	exitbutton=uicontrol(selectfig,'style','pushbutton',...
		'Position',[400 530 250 50],'string','DONE',...
		'callback','uiresume(gcbf)');
	disp('Click done after highlighting all MATCH files');
	uiwait(selectfig);

	selected_files=get(test,'value');
	close(selectfig);


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SCORE COMPUTATION

% assumes data is stored in two main directories, match and nomatch
% scores from match  are boostrapped to form the match distribution, 
% these are then compared with no match

% hopefully a region is occupied by the matches such that we can choose
% cutoff for JSD and xcorr score could use MAD rule-of-thumb

% now grab the relevant waveforms and ISIs

ncells=unique(FILES.idx);

% how many cells do we have for the matching condition

match_wavescore=[];
match_isiscore=[];

nomatch_wavescore=[];
nomatch_isiscore=[];

test_matrix=[];
class=[];

ATTRIBUTES.daysince=[];
ATTRIBUTES.date=[];
ATTRIBUTES.filename={};
ATTRIBUTES.cellidx=[];

for i=1:length(ncells)

	% take the data from the first file, compare to all others,
	% add to the MATCH distribution

	idxs=find(FILES.idx==ncells(i));
	currclusters=FILES.clustidx(idxs);

	load(FILES.name{idxs(1)},'clust_spike_vec','clusterwindows','IFR');

	templatewave=mean(clusterwindows{currclusters(1)},2);
	templatewave=zscore(templatewave);
	templateifr=mean(IFR{1}{currclusters(1)});
	templateifr=zscore(templateifr);
	templateisi=get_isi(clust_spike_vec{1}{currclusters(1)});

	for j=2:length(idxs)

		[wavescore isiscore ifrscore]=get_scores(templatewave,templateisi,templateifr,...
			FILES.name{idxs(j)},currclusters(j),isibins);
		
		test_matrix=[test_matrix;[wavescore isiscore ifrscore]];

		ATTRIBUTES.daysince(end+1)=daysdif(FILES.date(idxs(1)),FILES.date(idxs(j)));
		ATTRIBUTES.date(end+1)=FILES.date(idxs(j));
		ATTRIBUTES.filename{end+1}=FILES.name{idxs(j)};
		ATTRIBUTES.cellidx(end+1)=FILES.idx(idxs(j));

	end

end

% compute the nomatch scores, can now include match data (just don't compare
% within cells)

save('testing_data.mat','test_matrix','FILES','ATTRIBUTES');

end

function [LIST CLUSTERID CONTAINER]=recursive_file_find(ROOTDIR,FILTER,DIALOG)

rootfolder=uigetdir(ROOTDIR,DIALOG);

subfolders_pre=genpath(rootfolder);
subfolders=regexp(subfolders_pre,pathsep,'split');

LIST={};
CLUSTERID=[];
CONTAINER={};

for i=1:length(subfolders)

	filelisting_pre=dir(subfolders{i});
	filelisting={filelisting_pre(:).name};

	% read files in order, look for first file with cluster in the filename

	for j=1:length(filelisting)
		
		% also get the cluster number!

		if findstr(filelisting{j},'cluster') 
			
			[junk1 junk2 junk3 match]=regexp(filelisting{j},'cluster_([0-9]+)','match');
			CLUSTERID(end+1)=str2num(filelisting{j}(match{1}(1):match{1}(2)));
			CONTAINER{end+1}=subfolders{i};

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
		get_scores(TEMPLATEWAVE,TEMPLATEISI,TEMPLATEIFR,FILE,CLUSTER,ISIBINS)


load(FILE,'clust_spike_vec','clusterwindows','IFR');

scorewave=mean(clusterwindows{CLUSTER},2);
scorewave=zscore(scorewave);
scoreisi=get_isi(clust_spike_vec{1}{CLUSTER});
scoreifr=mean(IFR{1}{CLUSTER});
scoreifr=zscore(scoreifr);

[r,lags]=xcov(scorewave,TEMPLATEWAVE,'coeff');

WAVESCORE=max(r);

[r,lags]=xcov(TEMPLATEIFR,scoreifr,'none');

r=r./(norm(TEMPLATEIFR)*norm(scoreifr));
IFRSCORE=max(r);

ISISCORE=sqrt(kld(scoreisi,TEMPLATEISI,'jsd',1,'bins',ISIBINS));

end

