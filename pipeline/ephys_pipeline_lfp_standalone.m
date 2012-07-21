function ephys_lfp_standalone(PROCDIR,CONFIG)
%
%
% runs the multi-unit portion of the pipeline standalone
%
%

parameters=ephys_pipeline_readconfig(CONFIG);

fprintf('Parameters\n\n%-10s%-10s%-10s%-10s%-10s%-10s%-10s%-10s\n\n','NFFT','N','Overlap','Method','W','FS','Scale','Downsampling');
fprintf('%-10d%-10d%-10d%-10s%-10d%-10d%-10s%-10d\n\n\n',...
	parameters.lfp_nfft,parameters.lfp_n,parameters.lfp_overlap,parameters.lfp_method,...
	parameters.lfp_w,parameters.fs,parameters.lfp_scale,parameters.lfp_downsampling);

% load in agg data and histogram

disp('Computing LFP spectrograms and rasters...');
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

% check for peak in the mua that exceeds 4*std seems to be a good rule of thumb

disp('Computing spectrograms...');

% only compute spectrograms if the extraction is sufficiently long

% hops start at n/2 and increment n-overlap

samplemin=parameters.lfp_n/2+(parameters.lfp_minhops*(parameters.lfp_n-parameters.lfp_overlap));

if samples>samplemin
	ephys_visual_lfp_tf(EPHYS_DATA,HISTOGRAM,CHANNELS,'method',parameters.lfp_method,'savedir',PROCDIR,'lfp_min_f',parameters.lfp_min_f,'lfp_max_f',...
		parameters.lfp_max_f,'lfp_n',parameters.lfp_n,'lfp_nfft',parameters.lfp_nfft,'lfp_overlap',parameters.lfp_overlap,'scale',parameters.lfp_scale,...
		'lfp_ntapers',parameters.lfp_ntapers,'lfp_w',parameters.lfp_w,'SR',parameters.fs);
else
	disp([' Number of samples ' num2str(samples) ' less than samplemin ' num2str(samplemin)]);
end

disp('Computing LFP amplitudes...');

for i=1:length(parameters.freq_range)
	ephys_visual_lfp_amp(EPHYS_DATA,HISTOGRAM,CHANNELS,'freq_range',parameters.freq_range{i},'savedir',PROCDIR,'SR',parameters.fs);	
end
