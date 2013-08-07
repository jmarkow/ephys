function ephys_visual_sua(EPHYS_DATA,HISTOGRAM,CHANNELS,varargin)
%generates song-aligned single-unit rasters
%
%	ephys_visual_sua(EPHYS_DATA,HISTOGRAM,CHANNELS,varargin)
%
%	EPHYS_DATA
%	sound-aligned voltage traces from extracted_data.mat (should be the variable ephys_data)
%
%	HISTOGRAM
%	contour histogram returned by ephys_visual_histogram.m (or loaded from histogram.mat)
%
%	CHANNELS
%	channel labels (i.e. the channel that corresponds to a given element in the cell array
%	ephys_data) from extracted_data.mat
%
%	the following may be specified as parameter/value pairs:
%
%		tetrode_channels
%		channels may act as a putative n-trode (can specify an arbitrary number of channels),
%		defaults to manual cluster cutting
%
%		car_exclude
%		carelectrodes to exclude from noise estimate
%
%		fs
%		data sampling rate (default: 25e3)
%		
%		noise
%		noise rejection method ('car' for common average 'nn' for nearest neighbor, or 'none',
%		default: 'none')
%
%		freq_range
%		vector with two elements to specify the frequency range (one element specifies low pass, default: [500 8e3])
%
%		filt_type
%		filter type (default: 'bandpass', options 'low','high','bandpass')
%
%		filt_order
%		filter order (Butterworth, default: 6)
%
%		savedir
%		directory to store results (default: pwd)
%
%		min_f
%		lowermost frequency to display for contour histogram (default: 1e3)
%
%		max_f
%		uppermost frequency to display for contour histogram (default: 10e3)
%
%		auto_clust
%		perform automated cluster cutting via fitting a GMM through EM (default: 1)
%
%		sigma_t
%		multiple of variance estimate for automatic threshold setting (uses the Quiroga formulation, default: 4)
%
%		singletrials
%		number of singletrials to plot
%
%		savedir
%		directory to store results (default: pwd)
%
%		subtrials
%		vector of trials to include in the analysis (default: all trials)
%
%		modelselection
%		method of choosing the number of components for GMM clustering (default: icl, options, 
%		'BIC' for Bayes Information, 'mml' for minimum message length, 'icl' for BIC-ICL)
%
%		isi_cutoff
%		cutoff in ISI violation probability for candidate clusters (p<5 milliseconds, default: .05)
%
%		lratio_cutoff
%		cutoff in l-ratio for candidate clusters (default: .2)
%
%		isod_cutoff
%		cutoff in isolation distance for candidate clusters (default: 20)
%
%		snr_cutoff
%		cutoff in SNR (p2p of mean waveform/(6*noise_estimate)) for candidate clusters (default: 1.1)
%
%		spike_window
%		seconds before and after threshold crossing to store (in seconds, default: [.0005 .0005])
%
%		clust_start
%		vector of number of components to use in GMM clustering (default: [1:10])
%
%		pcs
%		number of principal components to use in GMM clustering (default: 2)
%
%		garbage
%		include a garbage-collecting uniform density in GMM clustering (0 or 1, default: 1)
%
%		smem
%		use split-and-merge algorithm for GMM clustering (0 or 1, default: 1)
%	
%		maxnoisetraces
%		maximum number of noise traces to use in spike whitening (default: 1e6)
%
% see also ephys_visual_mua.m,ephys_visual_lfp_amp.m,ephys_visual_lfp_tf.m,ephys_spike_cluster_auto.m,ephys_spike_clustergui_tetrode.m,ephys_spike_detect.m


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETER COLLECTION %%%%%%%%%%%%%%%%%

% TODO:  remove multi-channel support

if nargin<3
	error('ephysPipeline:suavis:notenoughparams','Need 3 arguments to continue, see documentation');
end

nparams=length(varargin);

if mod(nparams,2)>0
	error('ephysPipeline:argChk','Parameters must be specified as parameter/value pairs!');
end

channelboundary=[];
fs=25e3;
noise='none'; % none, nn for nearest neighbor, or car for common average
car_exclude=4;
savedir=pwd;
min_f=1; % min frequency to show for song histogram
max_f=10e3; % max frequency
hist_colors='jet'; % colormap for histogram
figtitle='';
freq_range=[800 11e3]; % bandpassing <10e3 distorted results, reasoning that >800 Hz is fine for spikes < 1ms long
filt_type='bandpass'; % high,low or bandpass
filt_order=6;
filt_name='b';
spikesort=1; % do we want to sort?
auto_clust=1; % 0 for manual cluster cutting (GUI), 1 for automated clustering
tetrode_channels=[];
sigma_t=4; % multiple of noise estimate for spike threshold (generally 3-4, using Quiroga's method)
jitter=10; % max jitter in samples for spike re-alignment (4 is reasonable
singletrials=5; % number of random single trials to plot per cluster
subtrials=[];
align='min'; % how to align spike waveforms can be min, max or com for center-of-mass
interpolate_fs=200e3; % 200 has worked best here
channels=CHANNELS;
smooth_rate=1e3;
sigma=.0025;
car_trim=40;
decomp_level=7;

savename=''; % add if doing multiple manual sorts, will append a name to the filename
isi_cutoff=.01; % percentage of ISI values <.001
lratio_cutoff=.2; % l-ratio cutoff from Reddish
isod_cutoff=20; % isolation distance defined by Harris et al.
snr_cutoff=8; % SNR definition from Ludwig et al. 2009 (J. Neurophys)
spike_window=[.0005 .0005];
trial_timestamps=[];
sort_fs=25e3;
cluststart=1:8;
pcs=2;
spikelimit=[];
garbage=1;
smem=1;
spikeworkers=1;
spikecut=1;
modelselection='icl';
maxnoisetraces=1e6;
wavelet_denoise=0;
noisewhiten=1;

% remove eps generation, too slow here...

colors={'b','r','g','c','m','y','k','r','g','b'};

for i=1:2:nparams
	switch lower(varargin{i})
		case 'fs'
			fs=varargin{i+1};
		case 'noise'
			noise=varargin{i+1};
		case 'sigma_t'
			sigma_t=varargin{i+1};
		case 'filtering'
			filtering=varargin{i+1};
		case 'savedir'
			savedir=varargin{i+1};
		case 'min_f'
			min_f=varargin{i+1};
		case 'max_f'
			max_f=varargin{i+1};
		case 'hist_colors'
			hist_colors=varargin{i+1};
		case 'car_exclude'
			car_exclude=varargin{i+1};
		case 'figtitle'
			figtitle=varargin{i+1};
		case 'filt_type'
			filt_type=varargin{i+1};
		case 'filt_name'
			filt_name=varargin{i+1};
		case 'decomp_level'
			decomp_level=varargin{i+1};
		case 'freq_range'
			freq_range=varargin{i+1};
		case 'spikesort'
			spikesort=varargin{i+1};
		case 'channels'
			channels=varargin{i+1};
		case 'tetrode_channels'
			tetrode_channels=varargin{i+1};
		case 'auto_clust'
			auto_clust=varargin{i+1};
		case 'sigma_t'
			sigma_t=varargin{i+1};
		case 'subtrials'
			subtrials=varargin{i+1};
		case 'align'
			align=varargin{i+1};
		case 'jitter'
			jitter=varargin{i+1};
		case 'interpolate_fs'
			interpolate_fs=varargin{i+1};
		case 'isi_cutoff'
			isi_cutoff=varargin{i+1};
		case 'isod_cutoff'
			isod_cutoff=varargin{i+1};
		case 'lratio_cutoff'
			lratio_cutoff=varargin{i+1};
		case 'snr_cutoff'
			snr_cutoff=varargin{i+1};
		case 'savename'
			savename=varargin{i+1};
		case 'spike_window'
			spike_window=varargin{i+1};
		case 'trial_timestamps'
			trial_timestamps=varargin{i+1};
		case 'sort_fs'
			sort_fs=varargin{i+1};
		case 'method'
			method=varargin{i+1};
		case 'maxnoisetraces'
			maxnoisetaces=varargin{i+1};
		case 'filt_order'
			filt_order=varargin{i+1};
		case 'cluststart'
			cluststart=varargin{i+1};
		case 'car_trim'
			car_trim=varargin{i+1};
		case 'pcs'
			pcs=varargin{i+1};
		case 'spikelimit'
			spikelimit=varargin{i+1};
		case 'spikecut'
			spikecut=varargin{i+1};
		case 'pcs'
			pcs=varargin{i+1};
		case 'garbage'
			garbage=varargin{i+1};
		case 'smem'
			smem=varargin{i+1};
		case 'spikeworkers'
			spikeworkers=varargin{i+1};
		case 'modelselection'
			modelselection=varargin{i+1};
		case 'wavelet_denoise'
			wavelet_denoise=varargin{i+1};
		case 'noisewhiten'
			noisewhiten=varargin{i+1};
		case 'decomp_level'
			decomp_level=varargin{i+1};
	end
end


[samples,ntrials,ncarelectrodes]=size(EPHYS_DATA);
proc_data=zeros(samples,ntrials,length(channels));

if isempty(subtrials)
	subtrials=1:ntrials;
end

if isempty(figtitle)
	figtitle=['Trials ' num2str(subtrials(1)) ' to ' num2str(subtrials(end))];
end


TIME=[1:samples]./fs; % time vector for plotting

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SIGNAL CONDITIONING %%%%%%%%%%%%%%%%

proc_data=ephys_denoise_signal(EPHYS_DATA,CHANNELS,channels,'method',noise,'car_exclude',car_exclude,'car_trim',car_trim);
proc_data=ephys_condition_signal(proc_data,'s','freq_range',...
	freq_range,'filt_type',filt_type,'filt_order',filt_order,'filt_name',filt_name,...
	'wavelet_denoise',wavelet_denoise,'decomp_level',decomp_level);

if ~isempty(tetrode_channels)
	tetrode_data=ephys_denoise_signal(EPHYS_DATA,CHANNELS,tetrode_channels,'method',noise,'car_exclude',car_exclude,'car_trim',car_trim);
	tetrode_data=ephys_condition_signal(tetrode_data,'s','freq_range',freq_range,'filt_type',filt_type,'filt_order',...
		filt_order,'filt_name',filt_name,'wavelet_denoise',wavelet_denoise,'decomp_level',decomp_level);
else
	tetrode_data=[];
end

clear EPHYS_DATA;

proc_data=proc_data(:,subtrials,:);

if ~isempty(tetrode_channels)
	tetrode_data=tetrode_data(:,subtrials,:);
end
[samples,ntrials,newchannels]=size(proc_data);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SPIKE DETECTION %%%%%%%%%%%%%%%%%%%%

disp('Entering spike detection...');
disp(['Alignment method:  ' align]);
disp(['Interpolate FS:  ' num2str(interpolate_fs)]);
% need a noise cutoff...let's start with 3*std or Quiroga's measure

disp(['Electrode ' num2str(channels)])

if ~isempty(tetrode_channels)
	disp(['Will use tetrode channels ' num2str(tetrode_channels) ' for sorting']);
end

% collect spikes

sort_data=cat(3,proc_data(:,:,1),tetrode_data);
totalspikes=0;
nchannels=size(sort_data,3);

for j=1:nchannels
	spikeless{j}=[];
end

for j=1:ntrials
	
	%spikes=[];
	%spikeless_data=[];

	spikethreshold=sigma_t*median(abs(sort_data(:,j,1))/.6745);
	%spikethreshold=10;
	% get the threshold crossings (based on first channel)

	spikes(j)=ephys_spike_detect(squeeze(sort_data(:,j,:)),spikethreshold,'fs',fs,'visualize','n','align',align,...
		'jitter',jitter,'window',spike_window);

	% get the spikeless data

	tmp=ephys_spike_removespikes(sort_data(:,j,1),spikes(j));
	spikeless{1}=[spikeless{1};tmp];

	% remove spikes from any other channels

	for k=2:nchannels
		tmp_thresh=sigma_t*median(abs(sort_data(:,j,k))/.6745);
		tmp_spikes=ephys_spike_detect(squeeze(sort_data(:,j,k)),tmp_thresh,'fs',fs,'visualize','n','align',align,...
			'window',spike_window);
		tmp=ephys_spike_removespikes(sort_data(:,j,k),tmp_spikes);
		spikeless{k}=[spikeless{k};tmp];
	end

	threshold(j)=spikethreshold;
	totalspikes=totalspikes+length(spikes(j).times);

end

disp([ num2str(totalspikes) ' total spikes']);

if totalspikes<spikelimit
	warning('ephyspipeline:sua:spikelimitnotexceeded','Not enough spikes to continue processing...');
	return;
end

clear sort_data;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% simplest way to draw raster...plot([sample;sample]
% vector of times [120 130 200;120 130 200], then next vector specifies bottom and top of tickmark
% e.g. [1.5 1.5 1.5;.5 .5 .5] for trial 1, etc.

% collapse spike windows into a single matrix, get cluster ids back and then 
% change color according to ID

disp('Generating figures...');
disp(['Will save to directory:  ' savedir]);

% scale pixels by time
% now we have clusterid and the trial for all spikes, we need to color them appropriately

[path,name,ext]=fileparts(fullfile(savedir));
savefilename=[ name '_sua_freqrange_' num2str(freq_range) '_electrode_' ];

if ~auto_clust & ~isempty(savename)
	savefilename=[ savefilename savename '_'];
end

savedir=fullfile(savedir,'sua');

if isdeployed
	savedir=fullfile(savedir,'pipeline');
elseif auto_clust
	savedir=fullfile(savedir,'auto');
else
	savedir=fullfile(savedir,'manual');
end

if ~exist(savedir,'dir')
	mkdir(fullfile(savedir));
end

disp(['Channel ' num2str(channels)]);

if spikesort
	if auto_clust
		[cluster.windows cluster.times cluster.trials cluster.isi cluster.stats... 
			cluster.outliers cluster.spikedata cluster.model]=...
			ephys_spike_cluster_auto(spikes,spikeless,...
			'fs',fs,'interpolate_fs',interpolate_fs,'proc_fs',sort_fs,...
			'maxnoisetraces',maxnoisetraces,'cluststart',cluststart,'pcs',pcs,...
			'workers',spikeworkers,'garbage',garbage,'smem',smem,'modelselection',...
			modelselection,'align',align,'noisewhiten',noisewhiten);
	else
		error('GUI sorting non-functional ATM, stay tuned...');
		[cluster.windows cluster.times cluster.trials cluster.isi cluster.stats]=...
			ephys_spike_clustergui_tetrode(spikewindows,spiketimes,spikeless,...
			'fs',fs,'wavelet_method',wavelet_method,'wavelet_mpca',wavelet_mpca,'interpolate_fs',interpolate_fs,...
			'template_cutoff',outlier_cutoff,'ranklimit',ranklimit,'garbage',garbage,'smem',smem);
		cluster.spikedata=cluster.windows;
	end

	cluster.parameters.fs=fs;
	cluster.parameters.interpolate_fs=interpolate_fs;
	cluster.parameters.sort_fs=sort_fs;
	cluster.parameters.threshold=threshold;
	cluster.parameters.tetrode_channels=tetrode_channels;
	cluster.parameters.spike_window=spike_window;

else
	clusterid=[];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% IFR, SMOOTH RATE %%%%%%%%%%%%%%%%%%%

uniq_clusters=1:length(cluster.windows);
nclust=length(uniq_clusters);

if ~isempty(cluster.windows)

	% cycle through each cluster id

	for j=1:nclust

		IFR{j}=zeros(ntrials,samples);

		for k=1:ntrials

			clusterspikes=cluster.times{j}(cluster.trials{j}==k);

			% IFR will use the current sampling rate

			IFR{j}(k,:)=ephys_ifr(round(clusterspikes),samples,fs);

		end
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nclust>1
	nplots=nclust;
else
	nplots=1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% STATS PLOTTING %%%%%%%%%%%%%%%%%%%%%

% plot the spike stats

savefilename_stats=[ savefilename num2str(channels)...
	'_stats_cluster_'];

% delete all old files

delete(fullfile(savedir,[ '*_sua_freqrange*electrode_' num2str(channels) '*']));

% delete old candidate files

delete(fullfile(savedir,[ 'candidate_unit_ch' num2str(channels) '*']));
delete(fullfile(savedir,[ '.candidate_unit_ch' num2str(channels) '*']));
delete(fullfile(savedir,[ '.sua_channels ' num2str(channels) '*']));

if ~isempty(cluster.windows)

	% delete old stats figures if they exist

	delete(fullfile(savedir,[savefilename_stats '*.png']));
	delete(fullfile(savedir,[savefilename_stats '*.eps']));

	% compute stats and generate stats figures

	cluster=ephys_cluster_autostats(cluster,spikeless,...
		'savefilename_stats',savefilename_stats,'savedir',savedir,'ntrials',ntrials,...
		'channelboundary',channelboundary,'spikecut',spikecut,'snr_cutoff',snr_cutoff,...
		'isod_cutoff',isod_cutoff,'lratio_cutoff',lratio_cutoff,'isi_cutoff',isi_cutoff,...
		'channels',channels);
	
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RASTER PLOTTING %%%%%%%%%%%%%%%%%%%%%

savefilename_raster=[ savefilename num2str(channels)...
	'_raster_cluster_'];

ax=[];

% delete any old single trial plots

if exist(fullfile(savedir,'singletrials',['ch' num2str(channels) ]),'dir')
	rmdir(fullfile(savedir,'singletrials',['ch' num2str(channels) ]),'s');
end

singleunit_raster(TIME,HISTOGRAM,cluster,proc_data,'figtitle',figtitle,'savefilename_raster',...
	savefilename_raster,'channels',channels,'savedir',savedir,'ntrials',ntrials,'min_f',min_f,...
	'max_f',max_f,'hist_colors',hist_colors,'singletrials',singletrials);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if spikesort & auto_clust 
	save(fullfile(savedir,['sua_channels ' num2str(channels) '.mat']),...
		'cluster','CHANNELS','channels','TIME','proc_data','freq_range',...
		'spikeless','IFR','subtrials','trial_timestamps');
elseif spikesort
	save(fullfile(savedir,['sua_channels ' num2str(channels) '.mat']),...
		'cluster','threshold','CHANNELS','channels','TIME','proc_data','freq_range',...
		'spikeless','IFR','subtrials','trial_timestamps');
else
	save(fullfile(savedir,['sua_channels ' num2str(channels) '.mat']),...
		'threshold','CHANNELS','channels','TIME','proc_data','freq_range','fs','trial_timestamps');
end

