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
FILES.date={};

% collect two primary file lists, allow the user to multi-select at the end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DATA COLLECTION

cellidx=0;
if nargin<1 | isempty(FILES)

	% get the match directory

	flag=0;

	while ~flag

		[filelist]=recursive_file_find(pwd,'sua_channels','Choose directories with cells to MATCH');

		% in this case we need to strip the last directory from the container

		for j=1:length(filelist)

			% strip the last two directories to get the "container"

			[path,name,ext]=fileparts(filelist{j});

			tokens=regexp(path,filesep,'split');

			newpath='';

			for k=1:length(tokens)-2
				newpath=[ newpath tokens{k} filesep ];
			end

			containers{j}=newpath;

			fid=fopen(fullfile(path,'..','cellinfo.txt'),'r');
			
			readdata=textscan(fid,'%s%[^\n]','commentstyle','#',...
				'delimiter','\t','MultipleDelimsAsOne',1);	
	
			fclose(fid);

			filedate{j}=readdata{2}{find(strcmpi(readdata{1},'date:'))};
			clusterids(j)=str2num(readdata{2}{find(strcmpi(readdata{1},'cluster:'))});

		end

		% unique root directories

		uniq_dirs=unique(containers);

		for j=1:length(filelist)

			FILES.name{end+1}=filelist{j};
			FILES.clustidx(end+1)=clusterids(j);
			FILES.idx(end+1)=find(strcmp(containers{j},uniq_dirs))+cellidx;
			FILES.date{end+1}=filedate{j};
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

	FILES.name=FILES.name(selected_files);
	FILES.clustidx=FILES.clustidx(selected_files);
	FILES.idx=FILES.idx(selected_files);
	FILES.date=FILES.date(selected_files);

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
ATTRIBUTES.date={};
ATTRIBUTES.filename={};
ATTRIBUTES.cellidx=[];
ATTRIBUTES.ifr={};

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

		load(FILES.name{idxs(j)},'clust_spike_vec','clusterwindows','IFR');

		scorewave=mean(clusterwindows{currclusters(j)},2);
		scorewave=zscore(scorewave);
		scoreisi=get_isi(clust_spike_vec{1}{currclusters(j)});
		scoreifr=mean(IFR{1}{currclusters(j)});
		scoreifr=zscore(scoreifr);

		[wavescore isiscore ifrscore]=get_scores(templatewave,templateisi,templateifr,...
			scorewave,scoreisi,scoreifr,isibins);

		test_matrix=[test_matrix;[wavescore isiscore ifrscore]];

		ATTRIBUTES.daysince(end+1)=daysdif(FILES.date{idxs(1)},FILES.date{idxs(j)});
		ATTRIBUTES.date{end+1}=FILES.date{idxs(j)};
		ATTRIBUTES.filename{end+1}=FILES.name{idxs(j)};
		ATTRIBUTES.cellidx(end+1)=FILES.idx(idxs(j));
		ATTRIBUTES.ifr{end+1}=downsample(smooth(scoreifr,.01*25e3),10); % probably want to smooth the estimate, trying raw first...
	end

end

% compute the nomatch scores, can now include match data (just don't compare
% within cells)

save('testing_data.mat','test_matrix','FILES','ATTRIBUTES');

end

function [LIST]=recursive_file_find(ROOTDIR,FILTER,DIALOG)

rootfolder=uigetdir(ROOTDIR,DIALOG);

subfolders_pre=genpath(rootfolder);
subfolders=regexp(subfolders_pre,pathsep,'split');

LIST={};
CLUSTERID=[];
CONTAINER={};

for i=1:length(subfolders)

	filelisting_pre=dir(subfolders{i});
	filelisting={filelisting_pre(:).name};

	for j=1:length(filelisting)

		if findstr(filelisting{j},FILTER)
			LIST{end+1}=fullfile(subfolders{i},filelisting{j});
			break;
		end

	end
end


end

function [WAVESCORE ISISCORE IFRSCORE]=...
	get_scores(TEMPLATEWAVE,TEMPLATEISI,TEMPLATEIFR,...
	SCOREWAVE,SCOREISI,SCOREIFR,ISIBINS)


[r,lags]=xcov(SCOREWAVE,TEMPLATEWAVE,'coeff');

WAVESCORE=max(r);

[r,lags]=xcov(TEMPLATEIFR,SCOREIFR,'none');

r=r./(norm(TEMPLATEIFR)*norm(SCOREIFR));
IFRSCORE=max(r);

ISISCORE=sqrt(kld(SCOREISI,TEMPLATEISI,'jsd',1,'bins',ISIBINS));

end

