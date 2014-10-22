function [fr]=ephys_cat_ratetracker(STORE_EPHYS,STORE_AUDIO,STORE_DATENUM,varargin)
%generates song-aligned single-unit rasters
%
%	ephys_visual_sua(EPHYS.data,HISTOGRAM,EPHYS.labels,varargin)
%
%	EPHYS
%	structure with the following fields
%
%	EPHYS.data
%	sound-aligned voltage traces from extracted_data.mat (should be the variable ephys_data)
%
%	EPHYS.labels
%	channel labels (i.e. the channel that corresponds to a given element in the cell array
%	ephys_data) from extracted_data.mat
%
%	HISTOGRAM
%	contour histogram returned by ephys_visual_histogram.m (or loaded from histogram.mat)
%
%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETER COLLECTION %%%%%%%%%%%%%%%%%

nparams=length(varargin);

if mod(nparams,2)>0
	error('ephysPipeline:argChk','Parameters must be specified as parameter/value pairs!');
end

noise='none'; % none, nn for nearest neighbor, or car for common average
car_exclude=[];

lfp_bands=[1 15;20 40;40 60;60 100;100 200;200 400;500 1e3];

fr_sigma_t=3; % multiple of noise estimate for spike threshold (generally 3-4, using Quiroga's method)
subtrials=[];
channels=[];
proc_fs=500;
fr_sigma=.005;
bound=.5;
audio_range=[1.7e3 4e3];
audio_thresh=.01;
audio_scale=.2;
ephys_range=[400 4.5e3];
savedir='';
savefile='cat_fr.mat';
rms=[];

% remove eps generation, too slow here...

for i=1:2:nparams
	switch lower(varargin{i})
		case 'noise'
			noise=varargin{i+1};
		case 'savedir'
			savedir=varargin{i+1};
		case 'car_exclude'
			car_exclude=varargin{i+1};
		case 'figtitle'
			figtitle=varargin{i+1};
		case 'filt_type'
			filt_type=varargin{i+1};
		case 'filt_name'
			filt_name=varargin{i+1};
		case 'freq_range'
			freq_range=varargin{i+1};
		case 'channels'
			channels=varargin{i+1};
		case 'fr_sigma_t'
			fr_sigma_t=varargin{i+1};
		case 'subtrials'
			subtrials=varargin{i+1};
		case 'trial_timestamps'
			trial_timestamps=varargin{i+1};
		case 'method'
			method=varargin{i+1};
		case 'filt_order'
			filt_order=varargin{i+1};
		case 'car_trim'
			car_trim=varargin{i+1};
		case 'lfp_bands'
			lfp_bands=varargin{i+1};
		case 'ephys_range'
			ephys_vange=varargin{i+1};
		case 'bound'
			bound=varargin{i+1};
		case 'savedir'
			savedir=varargin{i+1};
		case 'savefile'
			savefile=varargin{i+1};
		case 'rms'
			rms=varargin{i+1};
	end
end

[nsamples,ntrials,nchannels]=size(STORE_EPHYS.data);

if isempty(channels)
	channels=STORE_EPHYS.labels;
end

disp(['Will process channels: ' num2str(channels)]);

% decide which trials to include for threshold detection, check for silence after song

if ~isempty(STORE_AUDIO)
	
	[b,a]=ellip(5,.2,40,[audio_range]/(STORE_AUDIO.fs/2),'bandpass');
	STORE_AUDIO.data=filtfilt(b,a,double(STORE_AUDIO.data));

	audio_scale_smps=round(audio_scale*STORE_AUDIO.fs);
	coeffs=ones(1,audio_scale_smps)*1/audio_scale_smps;

	score=sqrt(filter(coeffs,1,STORE_AUDIO.data(end-bound*STORE_AUDIO.fs:end,:).^2));
	include=sum(score>audio_thresh)==0;

else	
	include=1:size(STORE_EPHYS.data,2);
end

if ~isempty(bound)
	disp(['Threshold bound:  ' num2str(bound)]);
	thresh_time=(nsamples-round(bound*STORE_EPHYS.fs)):nsamples;
else
	disp('No bound presented, using all trial data to compute threshold');
	thresh_time=1:nsamples;
end

for i=1:length(channels)

	disp(['Computing spike rate for channel ' num2str(channels(i)) '...']);

	if ~isempty(rms)
		tmprms=rms.silentpts.robust(find(rms.labels==channels(i)),:);
	else
		tmprms=[];
	end

	[ifr,trialrms,threshold,proc_data]=ephys_murate(STORE_EPHYS,'channels',channels(i),'sigma_t',...
		fr_sigma_t,'noise','car','bound',bound,'thresh_trials',find(include),'thresh_time',thresh_time,'rms',tmprms);

	kernedges=[-3*fr_sigma:1/STORE_EPHYS.fs:3*fr_sigma];
	kernel=normpdf(kernedges,0,fr_sigma);
	kernel=kernel./sum(kernel);

	downfact=STORE_EPHYS.fs/proc_fs;
	smooth_ifr=filter(kernel,1,ifr');

	fr(i).channel=channels(i);
	fr(i).time_series=single(downsample(smooth_ifr,downfact));	
	fr(i).threshold=single(threshold);
	fr(i).rms.edge=single(trialrms.edge);
	fr(i).rms.song=single(trialrms.song);
	fr(i).rms.time_series=single(downsample(trialrms.time_series,downfact));
	fr(i).rms.fs=proc_fs;
	fr(i).fs=proc_fs;
	fr(i).datenums=STORE_DATENUM;

end

if ~isempty(savedir)
	mkdir(savedir);
	save(fullfile(savedir,savefile),'fr','-v7.3');
end
