function ephys_lfp_standalone(PROCDIR,CONFIG)
%
%
% runs the multi-unit portion of the pipeline standalone
%
%

parameters=ephys_pipeline_readconfig(CONFIG);

fprintf('Parameters\n\n%-10s%-10s%-10s%-10s%-10s%-10s%-10s%-10s\n\n','NFFT','N','Overlap','Method','W','FS','Scale','Proc FS');
fprintf('%-10d%-10d%-10d%-10s%-10d%-10d%-10s%-10d\n\n\n',...
	parameters.lfp_nfft,parameters.lfp_n,parameters.lfp_overlap,parameters.lfp_method,...
	parameters.lfp_w,parameters.fs,parameters.lfp_scale,parameters.proc_fs);

% load in agg data and histogram

disp('Computing LFP spectrograms and rasters...');
disp(['Loading data from ' PROCDIR]);

if ~exist(fullfile(PROCDIR,'histogram.mat'),'file')
	histogram=[];
end

[agg_ephys,histogram,agg_file_datenum]=ephys_pipeline_dataload(PROCDIR,parameters.agg_filename);
[samples,ntrials,nchannels]=size(agg_ephys.data);

if ntrials>parameters.trial_max
	disp(['Exceeded trial max, truncating to ' num2str(parameters.trial_max)]);
	agg_ephys.data=agg_ephys.data(:,1:parameters.trial_max,:);
end

[samples,ntrials,nchannels]=size(agg_ephys.data);

% check for peak in the mua that exceeds 4*std seems to be a good rule of thumb

disp('Computing spectrograms...');

% only compute spectrograms if the extraction is sufficiently long

% hops start at n/2 and increment n-overlap

samplemin=parameters.lfp_n/2+(parameters.lfp_minhops*(parameters.lfp_n-parameters.lfp_overlap));


disp('Computing LFP amplitudes...');

for i=1:length(parameters.freq_range)
	ephys_visual_lfp_amp(agg_ephys,histogram,...
		'freq_range',parameters.freq_range{i},'savedir',PROCDIR,'fs',parameters.fs,'proc_fs',parameters.proc_fs,...
		'notch',parameters.lfp_notch,'notch_bw',parameters.lfp_notch_bw,'ripple',parameters.lfp_ripple,'attenuation',...
		parameters.lfp_attenuation);
end

if exist(fullfile(PROCDIR,'SLEEP_DATA'),'file')
	disp('Sleep data, skipping spectrograms...');
	return;
end

if samples>samplemin & ~exist(fullfile(PROCDIR,'SLEEP_DATA'),'file')

	if lower(parameters.lfp_method(1))=='c'

		disp('Computing LFP contours...');
		ephys_visual_lfp_tf_contour(agg_ephys,histogram,...
			'savedir',PROCDIR,'lfp_min_f',...
			parameters.lfp_min_f,'lfp_max_f',parameters.lfp_max_f,'lfp_n',...
			parameters.lfp_n,'lfp_nfft',parameters.lfp_nfft,'lfp_overlap',...
			parameters.lfp_overlap,'scale',parameters.lfp_scale,...
			'fs',parameters.fs,'proc_fs',parameters.proc_fs,'padding',parameters.lfp_padding,...
			'notch',parameters.lfp_notch,'notch_bw',parameters.lfp_notch_bw,'ripple',...
			parameters.lfp_ripple,'attenuation',parameters.lfp_attenuation);
	
	else

		disp('Computing LFP spectrograms...');
		ephys_visual_lfp_tf(agg_ephys,histogram,...
			'method',parameters.lfp_method,'savedir',PROCDIR,'lfp_min_f',...
			parameters.lfp_min_f,'lfp_max_f',parameters.lfp_max_f,'lfp_n',...
			parameters.lfp_n,'lfp_nfft',parameters.lfp_nfft,'lfp_overlap',...
			parameters.lfp_overlap,'scale',parameters.lfp_scale,...
			'lfp_ntapers',parameters.lfp_ntapers,'lfp_w',parameters.lfp_w,...
			'fs',parameters.fs,'proc_fs',parameters.proc_fs,'padding',parameters.lfp_padding,...
			'clipping',parameters.lfp_clipping,'notch',parameters.lfp_notch,'notch_bw',parameters.lfp_notch_bw,...
			'ripple',parameters.lfp_ripple,'attenuation',parameters.lfp_attenuation);

	end

else
	disp([' Number of samples ' num2str(samples) ' less than samplemin ' num2str(samplemin)]);
end

