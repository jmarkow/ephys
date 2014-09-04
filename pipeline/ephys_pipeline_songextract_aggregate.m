function ephys_songextract_aggregate(FILEDIR,SAVEDIR,CONFIG,SLEEPSTATUS)
%ephys_songextract_aggregate.m is part of the ephys pipeline
%and acts as a wrapper script for computing fields, multi-unit
%and (soon) single-unit rasters
%
%	ephys_songextract_aggregate(FILEDIR,SAVEDIR)
%
%	FILEDIR
%	directory with .mat files to process
%
%	SAVEDIR
%	base directory to store results
%

% do we want to run ephys scripts after aggregating?

% CONFIG FILE LIVES HERE!

if nargin<4
	SLEEPSTATUS=0;
end

if isdeployed
	SLEEPSTATUS=str2double(SLEEPSTATUS);
end

if nargin<3
	error('Need config file to continue!');
end

parameters=ephys_pipeline_readconfig(CONFIG);

% read config file one line at a time

listing=dir(fullfile(FILEDIR,'*.mat'));

% where do we want to store the output?

if length(listing)<1
	disp(['No files found in ' FILEDIR ' bailing...' ]);
	return;
end

% pre-allocate the matrix, store as a single

is_legacy=check_legacy(fullfile(FILEDIR,listing(1).name));

if is_legacy
	load(fullfile(FILEDIR,listing(1).name),'ephys_data','channels','fs');
	
	ephys.labels=channels;
	ephys.ports='';
	ephys.data=ephys_data;
	ephys.fs=fs;

	audio.fs=ephys.fs;

	clearvars ephys_data channels fs;
else
	load(fullfile(FILEDIR,listing(1).name),'ephys','audio','playback');
end

EPHYS.fs=ephys.fs;
AUDIO.fs=audio.fs;

EPHYS.t=ephys.t;
AUDIO.t=audio.t;

if exist('playback','var') & ~isempty(playback.data)	
	PLAYBACK.t=playback.t;
	PLAYBACK.fs=playback.fs;
end

% number of samples and channels

[samples,channels]=size(ephys.data);
ntrials=length(listing);

disp(['Aggregating trials in: ' FILEDIR]); 
disp(['Number of trials: ' num2str(ntrials)]);

if SLEEPSTATUS

	disp(['Identified sleep data...']);
	parameters.trial_win=parameters.trial_win_sleep;
end

dircount=1;

% extract in a sliding window, generate histograms and trigger mua, sua and lfp processing

savefun=@(filename,agg_ephys,agg_audio,agg_ttl,agg_playback,agg_file_datenum) ...
	save(filename,'agg_ephys','agg_audio','agg_ttl','agg_playback','agg_file_datenum','-v7.3');

for i=0:parameters.trial_win:ntrials
	
	if ~exist(fullfile(FILEDIR,[SAVEDIR '_' num2str(dircount)]),'dir')
		mkdir(fullfile(FILEDIR,[SAVEDIR '_' num2str(dircount)]));
	else
		if exist(fullfile(FILEDIR,[SAVEDIR '_' num2str(dircount)],'.done_extraction'),'file')
			disp(['Directory ' num2str(dircount) ' finished, continuing']);
			dircount=dircount+1;
			continue;
		end
	end

	currtrials=i+1:i+parameters.trial_win;
	currtrials(currtrials>ntrials)=[];

	all_labels=[];
	all_ports='';

	disp(['Identifying channels...']);

	for j=1:length(currtrials)

		if is_legacy
			load(fullfile(FILEDIR,listing(currtrials(j)).name),'channels');
			ephys.labels=channels;
			ephys.ports=repmat('A',[1 length(channels)]);
		else
			load(fullfile(FILEDIR,listing(currtrials(j)).name),'ephys');
		end

		% check for any inconsistency in channel labels

		for k=1:length(ephys.labels)

			label_chk=ephys.labels(k)==all_labels;
			port_chk=ephys.ports(k)==all_ports;

			% loop and if any channels are not included in the channel_label vector, include!

			if ~any(label_chk&port_chk)	
				all_labels=[all_labels ephys.labels(k)];
				all_ports=[all_ports ephys.ports(k)];
			end
		end
	end

	%[all_labels=sort(all_labels);

	disp(['Found channels:  ']);

	for j=1:length(all_labels)
		fprintf(1,'%i%s ',all_labels(j),all_ports(j));
	end

	fprintf(1,'\n');

	EPHYS.data=zeros(samples,length(currtrials),length(all_labels),'single');
	AUDIO.data=zeros(samples,length(currtrials));
	PLAYBACK.data=zeros(samples,length(currtrials));
	FILE_DATENUM=zeros(1,length(currtrials));
	TTL.data=zeros(samples,length(currtrials));

	EPHYS.labels=all_labels;
	EPHYS.ports=all_ports;

	for j=1:length(currtrials)

		disp(['Processing trial ' num2str(currtrials(j))]);

		if is_legacy
			load(fullfile(FILEDIR,listing(currtrials(j)).name),'ephys_data','mic_data','channels','fs','start_datenum','ttl_data');
			
			ephys.data=ephys_data;
			ephys.fs=fs;
			ephys.labels=channels;

			audio.data=mic_data;
			audio.fs=fs;

			ttl.data=ttl_data;
			ttl.fs=fs;

			file_datenum=start_datenum;

			clearvars ephys_data mic_data channels fs start_datenum ttl_data;

		else
			load(fullfile(FILEDIR,listing(currtrials(j)).name),'ephys','audio','file_datenum','ttl','playback');
		end

		AUDIO.data(:,j)=audio.data;
		
		if ~exist('ttl','var') | ~isfield(ttl,'data') | isempty(ttl.data) 
			ttl.data=ones(size(audio.data)).*NaN;
		end

		if ~exist('playback','var') |  ~isfield(playback,'data') | isempty(playback.data) 
			playback.data=ones(size(audio.data)).*NaN;
		end

		TTL.data(:,j)=ttl.data;
		PLAYBACK.data(:,j)=playback.data;

		for k=1:length(ephys.labels)

			ch_idx=find(ephys.labels(k)==all_labels);

			if isempty(ch_idx)
				continue;
			end

			EPHYS.data(:,j,ch_idx)=single(ephys.data(:,k));

		end

		if exist('file_datenum','var')
			FILE_DATENUM(j)=file_datenum;
		else

			% attempt to find the datenum

			tokens=regexpi(listing(currtrials(j)).name,parameters.delimiter,'split');
			FILE_DATENUM(j)=datenum([tokens{end-3} tokens{end-2}],'yymmddHHMMSS');

		end

	end

	% clear any unused data

	if sum(all(isnan(TTL.data)))==length(currtrials)
		TTL.data=[];
	end

	if sum(all(isnan(PLAYBACK.data)))==length(currtrials)
		PLAYBACK.data=[];
	end

	savefun(fullfile(FILEDIR,[SAVEDIR '_' num2str(dircount)],'aggregated_data.mat'),EPHYS,AUDIO,TTL,PLAYBACK,FILE_DATENUM);

	if SLEEPSTATUS

		% write out a file to indicate that this is sleep data

		fid=fopen(fullfile(FILEDIR,[SAVEDIR '_' num2str(dircount)],'SLEEP_DATA'),'w');
		fclose(fid);

	else
		histogram=ephys_visual_histogram(AUDIO.data,'fs',AUDIO.fs);
		save(fullfile(FILEDIR,[SAVEDIR '_' num2str(dircount)],'histogram.mat'),'histogram','-v7.3');
	end

	% touch signals to multi-unit, single-unit, and fields processing

	if ntrials>parameters.trial_min

		% multi-unit signal

		fid=fopen(fullfile(FILEDIR,[SAVEDIR '_' num2str(dircount)],'.mua_signal'),'w');
		fclose(fid);

		% signal-unit signal

		fid=fopen(fullfile(FILEDIR,[SAVEDIR '_' num2str(dircount)],'.sua_signal'),'w');
		fclose(fid);

		% lfp signal

		fid=fopen(fullfile(FILEDIR,[SAVEDIR '_' num2str(dircount)],'.lfp_signal'),'w');
		fclose(fid);

	end

	if i+parameters.trial_win<=ntrials

		% don't process again if we have a full window

		fid=fopen(fullfile(FILEDIR,[SAVEDIR '_' num2str(dircount)],'.done_extraction'),'w');
		fclose(fid);

	end

	dircount=dircount+1;

end

clearvars EPHYS AUDIO;

% if there aren't too many directories, go ahead and concatenate the data

if (dircount-1)<=parameters.cat_limit

	disp('Concatenating directory data...');

	ephys_cat_collect_data(FILEDIR);
	
	load(fullfile(FILEDIR,'ephys_cat','cat_data.mat'),'cat_audio','cat_ephys','cat_file_datenum');

	disp('Computing firing rate across all trials on all channels...');

	ephys_cat_ratetracker(cat_ephys,cat_audio,cat_file_datenum,'savedir',fullfile(FILEDIR,'ephys_cat'),...
		'channels',cat_ephys.labels);

	clearvars cat_audio cat_ephys cat_file_datenum;

	load(fullfile(FILEDIR,'ephys_cat','cat_fr.mat'),'fr');

	disp('Generating plots...');

	ephys_visual_cat_mutrack(fr,'savedir',fullfile(FILEDIR,'ephys_cat'));

end
