function [cluster spikeless]=ephys_spikesort(EPHYS_DATA,varargin)
%spike sorting without automated visualizations (as in ephys_visual_sua)
%
%	cluster=ephys_spikesort(EPHYS_DATA,varargin)
%
%	EPHYS_DATA
%	Data for spike sorting (samples x trials)
%
%	the following may be specified as parameter/value pairs:
%
%		interpolate_f
%		interpolation factor used for spike alignment (cubic splines, default: 8, i.e. upsample by factor of 8)
%
%		sort_f
%		downsample factor after alignment for clustering (default: interpolate_f, i.e. spikes upsampled then downsampled to original fs)
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
%		auto_clust
%		perform automated cluster cutting via fitting a GMM through EM (default: 1)
%
%		sigma_t
%		multiple of variance estimate for automatic threshold setting (uses the Quiroga formulation, default: 4)
%
%		modelselection
%		method of choosing the number of components for GMM clustering (default: icl, options, 
%		'BIC' for Bayes Information, 'mml' for minimum message length, 'icl' for BIC-ICL)
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
%		use split-and-merge algorithm for GMM clustering (0, 1, or 2 forr FREE smem, default: 1)
%	
%		maxnoisetraces
%		maximum number of noise traces to use in spike whitening (default: 1e6)
%
%		align_method
%		method for spike alignment ('min','max', or 'com', default: 'min');
%		
%		noise
%		denoising method ('car' for common average, 'none' for none, default: 'none');
%
%		jitter
%		maximum jitter for spike realignment in samples (default: 10)
%
%		noisewhiten
%		whiten noise using Cholesky decomposition (0 or 1, default: 1)
%
%		wavelet_denoise
%		denoise using wavelets (0 or 1, default: 0)
%
%		decomp_level
%		parameter for wavelet denoising (default: 7)
%
%
% see also ephys_visual_sua.m,ephys_spike_cluster_auto.m,ephys_spike_clustergui.m,ephys_spike_detect.m


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETER COLLECTION %%%%%%%%%%%%%%%%%

% TODO:  remove multi-channel support

if nargin<1
	error('ephysPipeline:suavis:notenoughparams','Need 1 arguments to continue, see documentation');
end

if isvector(EPHYS_DATA)
	EPHYS_DATA=EPHYS_DATA(:);
end

nparams=length(varargin);

if mod(nparams,2)>0
	error('ephysPipeline:argChk','Parameters must be specified as parameter/value pairs!');
end

fs=25e3;
noise='none'; % none, nn for nearest neighbor, or car for common average
car_exclude=[];

% 300 Hz E high-pass, see Quiroga et al. 2013

freq_range=[400]; % bandpassing <10e3 distorted results, reasoning that >800 Hz is fine for spikes < 1ms long
filt_type='high'; % high,low or bandpass
filt_order=3;
filt_name='e';
auto_clust=1; % 0 for manual cluster cutting (GUI), 1 for automated clustering
noise='none'; % none, nn for nearest neighbor, or car for common average

tetrode_channels=[];
sigma_t=4; % multiple of noise estimate for spike threshold (generally 3-4, using Quiroga's method)
jitter=10; % max jitter in samples for spike re-alignment (4 is reasonable
subtrials=[];
align_method='min'; % how to align spike waveforms can be min, max or com for center-of-mass
method='n';
interpolate_f=8; % interpolate factor
sort_f=[]; % if empty, downsamples back to original fs

smooth_rate=1e3;
car_trim=40;
decomp_level=7;

spike_window=[.0005 .0005];
cluststart=1:8;
pcs=2;
garbage=1;
smem=1;

spikeworkers=1;
modelselection='icl';
maxnoisetraces=1e6;
wavelet_denoise=0;
noisewhiten=1;

% remove eps generation, too slow here...

for i=1:2:nparams
	switch lower(varargin{i})
		case 'fs'
			fs=varargin{i+1};
		case 'noise'
			noise=varargin{i+1};
		case 'sigma_t'
			sigma_t=varargin{i+1};
		case 'car_exclude'
			car_exclude=varargin{i+1};
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
		case 'auto_clust'
			auto_clust=varargin{i+1};
		case 'align_method'
			align_method=varargin{i+1};
		case 'jitter'
			jitter=varargin{i+1};
		case 'interpolate_f'
			interpolate_f=varargin{i+1};
		case 'spike_window'
			spike_window=varargin{i+1};
		case 'sort_f'
			sort_f=varargin{i+1};
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

if isempty(sort_f)
	sort_f=interpolate_f;
	disp(['Setting sort downsample factor to spike upsample factor:  ' num2str(sort_f)]);
end

interpolate_fs=fs*interpolate_f;
sort_fs=interpolate_fs/sort_f;

[samples,ntrials,ncarelectrodes]=size(EPHYS_DATA);
if ncarelectrodes==1 & strcmp(noise,'car')
	disp('Turning off CAR, number of electrodes is 1');
	noise='none';
	car_exclude=[];
end

channels=1:ncarelectrodes;
proc_data=zeros(samples,ntrials,length(channels));

TIME=[1:samples]./fs; % time vector for plotting

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SIGNAL CONDITIONING %%%%%%%%%%%%%%%%

proc_data=ephys_denoise_signal(EPHYS_DATA,channels,channels,'method',noise,'car_exclude',car_exclude,'car_trim',car_trim);
proc_data=ephys_condition_signal(proc_data,'s','freq_range',...
	freq_range,'filt_type',filt_type,'filt_order',filt_order,'filt_name',filt_name,...
	'wavelet_denoise',wavelet_denoise,'decomp_level',decomp_level);


if ~isempty(tetrode_channels)
	tetrode_data=ephys_denoise_signal(EPHYS_DATA,channels,tetrode_channels,'method',noise,'car_exclude',car_exclude,'car_trim',car_trim);
	tetrode_data=ephys_condition_signal(tetrode_data,'s','freq_range',freq_range,'filt_type',filt_type,'filt_order',...
		filt_order,'filt_name',filt_name,'wavelet_denoise',wavelet_denoise,'decomp_level',decomp_level);
else
	tetrode_data=[];
end

clear EPHYS_DATA;

if ~isempty(tetrode_channels)
	tetrode_data=tetrode_data(:,subtrials,:);
end

[samples,ntrials,newchannels]=size(proc_data);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SPIKE DETECTION %%%%%%%%%%%%%%%%%%%%

disp('Entering spike detection...');
disp(['Alignment method:  ' align_method]);
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

	spikethreshold=sigma_t*median(abs(sort_data(:,j,1))/.6745);
	%spikethreshold=10;
	% get the threshold crossings (based on first channel)

	spikes(j)=ephys_spike_detect(squeeze(sort_data(:,j,:)),spikethreshold,'fs',fs,'visualize','n','align_method',align_method,...
		'jitter',jitter,'window',spike_window,'method',method);

	% get the spikeless data

	tmp=ephys_spike_removespikes(sort_data(:,j,1),spikes(j));
	spikeless{1}=[spikeless{1};tmp];

	% remove spikes from any other channels

	for k=2:nchannels
		tmp_thresh=sigma_t*median(abs(sort_data(:,j,k))/.6745);
		tmp_spikes=ephys_spike_detect(squeeze(sort_data(:,j,k)),tmp_thresh,'fs',fs,'visualize','n','align_method',align_method,...
			'window',spike_window);
		tmp=ephys_spike_removespikes(sort_data(:,j,k),tmp_spikes);
		spikeless{k}=[spikeless{k};tmp];
	end

	threshold(j)=spikethreshold;
	totalspikes=totalspikes+length(spikes(j).times);

end

disp([ num2str(totalspikes) ' total spikes']);
clear sort_data;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


disp(['Channel ' num2str(channels)]);

cluster.parameters.fs=fs;
cluster.parameters.interpolate_fs=interpolate_fs;
cluster.parameters.sort_fs=sort_fs;
cluster.parameters.threshold=threshold;
cluster.parameters.tetrode_channels=tetrode_channels;
cluster.parameters.spike_window=spike_window;
cluster.parameters.align_method=align_method;

if auto_clust
	[cluster.windows cluster.times cluster.trials cluster.isi cluster.stats... 
		cluster.outliers cluster.spikedata cluster.model]=...
		ephys_spike_cluster_auto(spikes,spikeless,...
			'fs',fs,'interpolate_fs',interpolate_fs,'proc_fs',sort_fs,...
			'maxnoisetraces',maxnoisetraces,'cluststart',cluststart,'pcs',pcs,...
			'workers',spikeworkers,'garbage',garbage,'smem',smem,'modelselection',...
			modelselection,'align_method',align_method,'noisewhiten',noisewhiten);
else
	[cluster.windows cluster.times cluster.trials cluster.isi cluster.stats...
		cluster.outliers cluster.spikedata cluster.model]=...
		ephys_spike_clustergui(spikes,spikeless,cluster.parameters,...
			'fs',fs,'interpolate_fs',interpolate_fs,'proc_fs',sort_fs,...
			'maxnoisetraces',maxnoisetraces,'cluststart',cluststart,'pcs',pcs,...
			'workers',spikeworkers,'garbage',garbage,'smem',smem,'modelselection',...
			modelselection,'align_method',align_method,'noisewhiten',noisewhiten);
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

cluster.IFR=IFR;
