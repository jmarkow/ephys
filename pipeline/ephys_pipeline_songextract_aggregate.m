function ephys_songextract_aggregate(FILEDIR,SAVEDIR,CONFIG)
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

if ~exist(fullfile(FILEDIR,SAVEDIR),'dir')
	mkdir(fullfile(FILEDIR,SAVEDIR));
end

% pre-allocate the matrix, store as a single

load(fullfile(FILEDIR,listing(1).name),'ephys_data','channels','fs');

% number of samples and channels

%channel_labels=channels;
[samples,channels]=size(ephys_data);
ntrials=length(listing);

disp(['Number of trials ' num2str(ntrials)]);

% assume the first trial has a consistent number of channels


% matrix is samples,trials,channels

channel_labels=[];

% check for all possible channels across the whole day, a matrix will be filled with zeros if the channel
% gets knocked out somehow...

for i=1:ntrials
	load(fullfile(FILEDIR,listing(i).name),'channels');

	% check for any inconsistency in channel labels

	for j=1:length(channels)

		% loop and if any channels are not included in the channel_label vector, include!

		if ~any(channels(j)==channel_labels)
			channel_labels=[channel_labels channels(j)];
		end
	end
end

channel_labels=sort(channel_labels);
disp(['Found channels ' num2str(channel_labels)]);

% set the hard max to something like 600, we can always turn off later to extract all data
% > 600 seems to lean pretty hard on the RAM (>3-4 GB)

if ntrials>parameters.trial_hard_max
	disp(['N Trials exceeds hardmax, truncating to ' num2str(parameters.trial_hard_max)]);
	ntrials=parameters.trial_hard_max;
end

EPHYS_DATA=zeros(samples,ntrials,length(channel_labels),'single');
MIC_DATA=zeros(samples,ntrials);

% we may need a hard max in case the number of trials exceeds ~600 trials (starts to eat up an inordinate amount of memory)
% we can always override this if need be

% also causes the pipeline to spend an inordinate amount of time with file i/o

for i=1:ntrials
	
	load(fullfile(FILEDIR,listing(i).name),'ephys_data','mic_data','channels','fs');
	MIC_DATA(:,i)=mic_data;

	for j=1:length(channels)
		
		ch_idx=find(channels(j)==channel_labels);

		if isempty(ch_idx)
			continue;
		end

		EPHYS_DATA(:,i,ch_idx)=single(ephys_data(:,j));
	
	end
end

CHANNELS=channel_labels;
save(fullfile(FILEDIR,SAVEDIR,'aggregated_data.mat'),'EPHYS_DATA','MIC_DATA','CHANNELS');

% aggregator also computes histogram for now...

if ntrials>parameters.trial_max
	disp(['Exceeded trial max, truncating to ' num2str(parameters.trial_max)]);
	MIC_DATA=MIC_DATA(:,1:parameters.trial_max);
end

histogram=ephys_visual_histogram(MIC_DATA,'savedir',fullfile(FILEDIR,SAVEDIR));

% touch signals to multi-unit, single-unit, and fields processing

% add a script to detect if spikes exist with SNR>1.1

if ntrials>parameters.trial_min

	% multi-unit signal

	fid=fopen(fullfile(FILEDIR,SAVEDIR,'.mua_signal'),'w');
	fclose(fid);

	% signal-unit signal

	fid=fopen(fullfile(FILEDIR,SAVEDIR,'.sua_signal'),'w');
	fclose(fid);

	% lfp signal

	fid=fopen(fullfile(FILEDIR,SAVEDIR,'.lfp_signal'),'w');
	fclose(fid);

end
