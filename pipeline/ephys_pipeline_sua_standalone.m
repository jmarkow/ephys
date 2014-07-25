function ephys_lfp_standalone(PROCDIR,CONFIG)
%
%
% runs the single-unit portion of the pipeline standalone
%
%

parameters=ephys_pipeline_readconfig(CONFIG);

fprintf('Parameters\n\n%-10s%-10s%-10s%-10s\n\n','Sigma_t','Align','FS');
fprintf('%-10d%-10s%-10d%-10s\n\n\n',parameters.sigma_t,parameters.spike_align,parameters.spike_fs);

% load in agg data and histogram

disp('Computing single unit SNR and rasters');
disp(['Loading data from ' PROCDIR]);

START_DATENUM=[];

[agg_ephys,histogram,agg_file_datenum]=ephys_pipeline_dataload(PROCDIR,parameters.agg_filename);
[samples,ntrials,nchannels]=size(agg_ephys.data);

if ntrials>parameters.trial_max
	disp(['Exceeded trial max, truncating to ' num2str(parameters.trial_max)]);
	ephys.data=ephys.data(:,1:parameters.trial_max,:);
end

[samples,ntrials,nchannels]=size(agg_ephys.data);
delete(fullfile(PROCDIR,'snr_channel_*'));

disp('Computing SNR on all channels');
SNR=ephys_pipeline_candidate_su(agg_ephys,'savedir',PROCDIR,'snr_threshold',parameters.snr_cutoff);

% check for contiguous windows with SNR>SNR_min

candidate_channels=[];
trials=[];

for i=1:nchannels
	
	snrvec=SNR(:,i);

	strvec=[0;snrvec>parameters.snr_cutoff;0]';

	startidx=strfind(strvec,[0 1]);
	stopidx=strfind(strvec,[1 0])-1;
	
	[maxregion loc]=max(stopidx-startidx);

	if length(maxregion>1)
		maxregion=maxregion(1);
		loc=loc(1);
	end

	if maxregion>parameters.snr_trials
		
		candidate_channels=[candidate_channels;agg_ephys.labels(i)];
		trials=[trials;[startidx(loc) stopidx(loc)]];
	
	end

end

% uncommnent next line to sort everything ...
% candidate_channels=setdiff(CHANNELS,parameters.ref_channel);

for i=1:length(candidate_channels)
	try 
		ephys_visual_sua(agg_ephys,histogram,...
			'channels',candidate_channels(i),'sort',1,'auto_clust',1,...
			'savedir',PROCDIR,'align',parameters.spike_align,'sigma_t',parameters.sigma_t,...
			'jitter',parameters.jitter,'interpolate_f',parameters.spike_interpolate_f,'filt_type',parameters.spike_filt_type,...
			'freq_range',parameters.spike_freq_range,'snr_cutoff',parameters.unit_snr_cutoff,...
			'lratio_cutoff',parameters.unit_lratio_cutoff,'isod_cutoff',parameters.unit_isod_cutoff,...
			'isi_cutoff',parameters.unit_isi_cutoff,'trial_timestamps',agg_file_datenum,...
			'sort_f',parameters.spike_sort_f,'maxnoisetraces',parameters.spike_exact_maxnoisetraces,...
			'filt_order',parameters.spike_filt_order,...
			'pcs',parameters.spike_pcs,'cluststart',parameters.spike_cluststart,...
			'car_exclude',parameters.spike_car_exclude,'noise',parameters.spike_noise_method,'car_trim',parameters.spike_car_trim,...
			'spike_window',parameters.spike_window,'spikelimit',parameters.spike_spikelimit,'modelselection',...
			parameters.spike_modelselection,'smem',parameters.spike_smem,'garbage',parameters.spike_garbage,...
			'spikeworkers',parameters.spike_workers);
	catch
		disp(['Could not process channel ' num2str(candidate_channels(i)) ', continuing...']);
	end

end

