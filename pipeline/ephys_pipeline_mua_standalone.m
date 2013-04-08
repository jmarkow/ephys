function ephys_mua_standalone(PROCDIR,CONFIG)
%
%
% runs the multi-unit portion of the pipeline standalone
%
%

parameters=ephys_pipeline_readconfig(CONFIG);

fprintf('Parameters\n\n%-10s%-10s\n\n','FS','Downsampling');
fprintf('%-10d%-10d\n\n\n',parameters.fs,parameters.mua_downsampling);

% load in agg data and histogram

disp('Computing multi-unit rasters...');
disp(['Loading data from ' PROCDIR]);

if ~exist(fullfile(PROCDIR,'histogram.mat'),'file')
	HISTOGRAM=[];
end

try
	load(fullfile(PROCDIR,'aggregated_data.mat'),'EPHYS_DATA','CHANNELS');
	if exist(fullfile(PROCDIR,'histogram.mat'),'file')
		load(fullfile(PROCDIR,'histogram.mat'),'HISTOGRAM');
	end
catch
	try
		disp('Pausing for 60 seconds and will retry');
		pause(60);
		load(fullfile(PROCDIR,'aggregated_data.mat'),'EPHYS_DATA','CHANNELS');
		if exist(fullfile(PROCDIR,'histogram.mat'),'file')
			load(fullfile(PROCDIR,'histogram.mat'),'HISTOGRAM');
		end
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

mua=ephys_visual_mua(EPHYS_DATA,HISTOGRAM,CHANNELS,'savedir',PROCDIR,'fs',parameters.fs);

% check for peak in the mua that exceeds 4*std seems to be a good rule of thumb

if exist(fullfile(PROCDIR,'SLEEP_DATA'),'file')
	disp('Sleep directory, skipping projection neuron checks');
	return;
end

delete(fullfile(PROCDIR,'proj_channel_*'));
trial_edges=[1:parameters.trial_window:ntrials ntrials];

% slide a window window_trials long and detect peaks

for i=1:length(CHANNELS)

	% sum across mult-unit, check for 4.5-5*std+mean peaks, if they exist include in the candidate
	% list, also mark as potential projection neurons

	% change to include sliding window across n trials

	for j=1:length(trial_edges)-1

		sumvec=sum(mua.image(trial_edges(j):trial_edges(j+1),:,i));
		hits=find(sumvec>parameters.proj_cutoff*std(sumvec)+mean(sumvec));

		if ~isempty(hits)
			disp(['Possible proj. cell detected on ' num2str(CHANNELS(i)) ... 
				' trials ' num2str([trial_edges(j:j+1)]) ]);
			fid=fopen(fullfile(PROCDIR,['proj_channel_' num2str(CHANNELS(i)) ...
				'_trials_' num2str(trial_edges(j)) '_' num2str(trial_edges(j+1))	]),'w');
			fclose(fid);
		end
	end

end

