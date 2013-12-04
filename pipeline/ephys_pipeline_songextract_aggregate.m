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

load(fullfile(FILEDIR,listing(1).name),'ephys_data','channels','fs');

% number of samples and channels

[samples,channels]=size(ephys_data);
ntrials=length(listing);

disp(['Number of trials ' num2str(ntrials)]);

if SLEEPSTATUS
	parameters.trial_win=parameters.trial_win_sleep;
end

dircount=1;

% extract in a sliding window, generate histograms and trigger mua, sua and lfp processing

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

	channel_labels=[];

	for j=1:length(currtrials)
		load(fullfile(FILEDIR,listing(currtrials(j)).name),'channels');

		% check for any inconsistency in channel labels

		for k=1:length(channels)

			% loop and if any channels are not included in the channel_label vector, include!

			if ~any(channels(k)==channel_labels)
				channel_labels=[channel_labels channels(k)];
			end
		end
	end

	channel_labels=sort(channel_labels);
	disp(['Found channels ' num2str(channel_labels)]);

	EPHYS_DATA=zeros(samples,length(currtrials),length(channel_labels),'single');
	MIC_DATA=zeros(samples,length(currtrials));
	START_DATENUM=zeros(1,length(currtrials));
	TTL_DATA=zeros(samples,length(currtrials));

	for j=1:length(currtrials)

		load(fullfile(FILEDIR,listing(currtrials(j)).name),'ephys_data','mic_data','channels','fs','start_datenum','ttl_data');
		MIC_DATA(:,j)=mic_data;
		
		if ~exist('ttl_data','var') | isempty(ttl_data)
			ttl_data=zeros(size(mic_data));
		end

		TTL_DATA(:,j)=ttl_data;

		for k=1:length(channels)

			ch_idx=find(channels(k)==channel_labels);

			if isempty(ch_idx)
				continue;
			end

			EPHYS_DATA(:,j,ch_idx)=single(ephys_data(:,k));

		end

		if exist('start_datenum','var')
			START_DATENUM(j)=start_datenum;
		else

			% attempt to find the datenum

			tokens=regexpi(listing(currtrials(j)).name,parameters.delimiter,'split');

			if length(tokens)<6
				continue;
			end

			% fourth is date

			if strcmpi(tokens{5}(1:3),'ttl')
				idx=6;
			else
				idx=5;
			end

			START_DATENUM(j)=datenum([tokens{idx} tokens{idx+1}],'yymmddHHMMSS');

		end

	end

	CHANNELS=channel_labels;
	save(fullfile(FILEDIR,[SAVEDIR '_' num2str(dircount)],'aggregated_data.mat'),...
		'EPHYS_DATA','fs','MIC_DATA','CHANNELS','START_DATENUM','TTL_DATA','-v7.3');

	if SLEEPSTATUS

		% write out a file to indicate that this is sleep data

		fid=fopen(fullfile(FILEDIR,[SAVEDIR '_' num2str(dircount)],'SLEEP_DATA'),'w');
		fclose(fid);
	else
		histogram=ephys_visual_histogram(MIC_DATA,'savedir',fullfile(FILEDIR,[SAVEDIR '_' num2str(dircount)]));
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
