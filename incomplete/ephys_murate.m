function [IFR,RMS,THRESHOLD,PROC_DATA]=ephys_murate(EPHYS,varargin)
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

% TODO:  remove multi-channel support

if isvector(EPHYS.data)
	EPHYS.data=EPHYS.data(:);
end

nparams=length(varargin);

if mod(nparams,2)>0
	error('ephysPipeline:argChk','Parameters must be specified as parameter/value pairs!');
end

channelboundary=[];
noise='none'; % none, nn for nearest neighbor, or car for common average
car_exclude=[];
savedir=pwd;

min_f=1; % min frequency to show for song histogram
max_f=10e3; % max frequency
hist_colors='jet'; % colormap for histogram

figtitle='';
freq_range=[400 4e3]; % bandpassing <10e3 distorted results, reasoning that >800 Hz is fine for spikes < 1ms long
downfilt=5e3;

filt_type='bandpass'; % high,low or bandpass
filt_order=6;
filt_name='e';

sigma_t=3; % multiple of noise estimate for spike threshold (generally 3-4, using Quiroga's method)
subtrials=[];
channels=EPHYS.labels;

bound=.2;
car_trim=0;
decomp_level=7;

interpolate_f=8; % interpolate factor
sort_f=1; % if empty, downsamples back to original fs

savename=''; % add if doing multiple manual sorts, will append a name to the filename
trial_timestamps=[];
proc_fs=10e3;
downfilt=5e3;
thresh_trials=[];
thresh_time=[];
rms=[];
threshold_smoothing=35;
fr_smoothing=.015;

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
		case 'sigma_t'
			sigma_t=varargin{i+1};
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
		case 'bound'
			bound=varargin{i+1};
		case 'proc_fs'
			proc_fs=varargin{i+1};
		case 'thresh_trials'
			thresh_trials=varargin{i+1};
		case 'thresh_time'
			thresh_time=varargin{i+1};
		case 'rms'
			rms=varargin{i+1};
		case 'threshold_smoothing'
			threshold_smoothing=varargin{i+1};
		case 'fr_smoothing'
			fr_smoothing=varargin{i+1};
	end
end

[nsamples,ntrials,nchannels]=size(EPHYS.data);

if isempty(subtrials)
	subtrials=1:ntrials;
end

if isempty(thresh_trials)
	thresh_trials=1:ntrials;
end

if isempty(thresh_time)
	thresh_time=1:nsamples;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SIGNAL CONDITIONING %%%%%%%%%%%%%%%%

PROC_DATA=ephys_denoise_signal(EPHYS.data,EPHYS.labels,channels,'method',noise,'car_exclude',car_exclude,'car_trim',car_trim);
clear EPHYS.data;

PROC_DATA=ephys_condition_signal(PROC_DATA,'s','freq_range',...
	freq_range,'filt_type',filt_type,'filt_order',filt_order,'filt_name',filt_name,'fs',proc_fs);
PROC_DATA=PROC_DATA(:,subtrials,:);
[nsamples,ntrials,newchannels]=size(PROC_DATA);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SPIKE DETECTION %%%%%%%%%%%%%%%%%%%%

disp('Entering spike detection...');
disp(['Electrode ' num2str(channels)])

% collect spikes
totalspikes=0;

% median filter the spikethreshold

if isempty(rms)
	
	threshdata=PROC_DATA(thresh_time,thresh_trials,1);
	spikethreshold=sigma_t*median(abs(threshdata(:))/.6745);
	disp(['Spike threshold:  ' num2str(spikethreshold)]);
	spikethreshold=repmat(spikethreshold,[1 ntrials]);

else
	if ~isempty(threshold_smoothing)
		spikethreshold=medfilt1(sigma_t*rms,min(threshold_smoothing,ntrials));
	else
		spikethreshold=sigma_t*rms;
	end
end

RMS.edge=zeros(1,ntrials);
RMS.song=zeros(1,ntrials);
RMS.time_series=zeros(nsamples,ntrials);
THRESHOLD=zeros(1,ntrials);

[nblanks formatstring]=fb_progressbar(100);
fprintf(1,['Progress:  ' blanks(nblanks)]);

for i=1:ntrials
	
	fprintf(1,formatstring,round((i/ntrials)*100));	

	spikes(i)=ephys_spike_detect(squeeze(PROC_DATA(:,i,:)),spikethreshold(i),'fs',EPHYS.fs,'visualize','n');

	%tmp=ephys_spike_removespikes(PROC_DATA(:,i,1),spikes(i));

	THRESHOLD(i)=spikethreshold(i);

	if ~isempty(bound)
		RMS.edge(i)=sqrt(mean(PROC_DATA(thresh_time,i,:).^2));
		RMS.song(i)=sqrt(mean(PROC_DATA(bound*EPHYS.fs:end-bound*EPHYS.fs,i,:).^2));
	else
		RMS.edge(i)=sqrt(mean(PROC_DATA(:,i,:).^2));
		RMS.song(i)=RMS.edge(i);
	end

	RMS.time_series(:,i)=sqrt(smooth(PROC_DATA(:,i,:).^2,round(fr_smoothing*EPHYS.fs)));
	totalspikes=totalspikes+length(spikes(i).times);

end

fprintf(1,'\n');
disp([ num2str(totalspikes) ' total spikes']);

IFR=zeros(ntrials,nsamples);

for i=1:ntrials

	tmp=spikes(i).times;

	% IFR will use the current sampling rate

	IFR(i,:)=ephys_ifr(tmp,nsamples,EPHYS.fs);
end

