function ephys_lfp_standalone(PROCDIR,CONFIG)
%
%
% runs the single-unit portion of the pipeline standalone
%
%

parameters=ephys_pipeline_readconfig(CONFIG);

fprintf('Parameters\n\n%-10s%-10s%-10s%-10s\n\n','Sigma_t','Align','FS','Comp.');
fprintf('%-10d%-10s%-10d%-10s\n\n\n',parameters.sigma_t,parameters.spike_align,parameters.spike_fs,parameters.mode_selection);

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
SNR=ephys_pipeline_candidate_su(EPHYS_DATA,CHANNELS,'savedir',PROCDIR,'snr_threshold',parameters.snr_cutoff);

candidate_channels=CHANNELS(SNR>=parameters.snr_cutoff); % channels with average SNR over cutoff for all trials
%	get processed

for i=1:length(candidate_channels)
	ephys_visual_sua(EPHYS_DATA,HISTOGRAM,CHANNELS,...
		'channels',candidate_channels(i),'sort',1,'auto_clust',1,...
		'savedir',PROCDIR,'align',parameters.spike_align,'sigma_t',parameters.sigma_t,'sr',parameters.fs,...
		'jitter',parameters.jitter,'wavelet_method',parameters.wavelet_sort,'wavelet_mpca',parameters.wavelet_mpca,...
		'wavelet_coeffs',parameters.wavelet_coeffs,'clust_choice',parameters.mode_selection,...
		'interpolate_fs',parameters.spike_fs,'filt_type',parameters.spike_filt_type,'freq_range',parameters.spike_freq_range);
end


