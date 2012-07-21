function ephys_lfp_standalone(PROCDIR,CONFIG)
%
%
% runs the single-unit portion of the pipeline standalone
%
%

parameters=ephys_pipeline_readconfig(CONFIG);

fprintf('Parameters\n\n%-10s%-10s%-10s\n\n','Sigma_t','Align','FS');
fprintf('%-10d%-10s%-10d\n\n\n',parameters.sigma_t,parameters.spike_align,parameters.fs);

% load in agg data and histogram

disp('Computing single unit SNR and rasters');
disp(['Loading data from ' PROCDIR]);

try
	load(fullfile(PROCDIR,'aggregated_data.mat'),'EPHYS_DATA','CHANNELS');
	load(fullfile(PROCDIR,'histogram.mat'),'HISTOGRAM');
catch
	try
		disp('Pausing for 60 seconds and will retry');
		pause(60);
		load(fullfile(PROCDIR,'aggregated_data.mat'),'EPHYS_DATA','CHANNELS');
		load(fullfile(PROCDIR,'histogram.mat'),'HISTOGRAM');
	catch

		disp('Could not properly load files, bailing!');
		return;
	end
end

[samples,ntrials,nchannels]=size(EPHYS_DATA);

if ntrials>parameters.trial_max
	disp(['Exceeded trial max, truncating to ' num2str(parameters.trial_max)]);
	EPHYS_DATA=EPHYS_DATA(:,1:parameters.trial_max,:);
end

[samples,ntrials,nchannels]=size(EPHYS_DATA);

delete(fullfile(PROCDIR,'snr_channel_*'));
disp('Computing SNR on all channels');
SNR=ephys_pipeline_candidate_su(EPHYS_DATA,CHANNELS,'savedir',PROCDIR);

candidate_channels=CHANNELS(SNR>=parameters.snr_cutoff); % channels with average SNR over cutoff for all trials
%	get processed

for i=1:length(candidate_channels)
	ephys_visual_sua(EPHYS_DATA,HISTOGRAM,CHANNELS,...
		'channels',candidate_channels(i),'sort',1,'auto_clust',1,...
		'savedir',PROCDIR,'align',parameters.spike_align,'sigma_t',parameters.sigma_t,'sr',parameters.fs);
end


