function IFR=ephys_murate(EPHYS,varargin)
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
freq_range=[400 6e3]; % bandpassing <10e3 distorted results, reasoning that >800 Hz is fine for spikes < 1ms long
filt_type='bandpass'; % high,low or bandpass
filt_order=6;
filt_name='e';

sigma_t=3; % multiple of noise estimate for spike threshold (generally 3-4, using Quiroga's method)
subtrials=[];
channels=EPHYS.labels;

car_trim=40;
decomp_level=7;

interpolate_f=8; % interpolate factor
sort_f=1; % if empty, downsamples back to original fs

savename=''; % add if doing multiple manual sorts, will append a name to the filename
trial_timestamps=[];

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

	end
end

[samples,ntrials,nchannels]=size(EPHYS.data);

if isempty(subtrials)
	subtrials=1:ntrials;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SIGNAL CONDITIONING %%%%%%%%%%%%%%%%

proc_data=ephys_denoise_signal(EPHYS.data,EPHYS.labels,channels,'method',noise,'car_exclude',car_exclude,'car_trim',car_trim);
proc_data=ephys_condition_signal(proc_data,'s','freq_range',...
	freq_range,'filt_type',filt_type,'filt_order',filt_order,'filt_name',filt_name,'fs',EPHYS.fs);
clear EPHYS.data;

size(proc_data)
proc_data=proc_data(:,subtrials,:);
[samples,ntrials,newchannels]=size(proc_data);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SPIKE DETECTION %%%%%%%%%%%%%%%%%%%%

disp('Entering spike detection...');
disp(['Electrode ' num2str(channels)])

% collect spikes
totalspikes=0;
threshdata=proc_data(end-.2*EPHYS.fs:end,:,1);

%spikethreshold=sigma_t*sqrt(mean(threshdata(:).^2))
spikethreshold=sigma_t*median(abs(threshdata(:))/.6745);
disp(['Spike threshold:  ' num2str(spikethreshold)]);

[nblanks formatstring]=fb_progressbar(100);
fprintf(1,['Progress:  ' blanks(nblanks)]);

for i=1:ntrials
	fprintf(1,formatstring,round((i/ntrials)*100));	

    %spikethreshold=50;
	%spikethreshold=sigma_t*median(abs(proc_data(end-.1*EPHYS.fs:end,:,i,1))/.6745)
	%spikethreshold=sigma_t*std(proc_data(:,i,1))
	%spikethreshold=sigma_t*sqrt(mean(proc_data(end-.1*EPHYS.fs:end,i,1).^2));
	%spikethreshold=sigma_t*sqrt(mean(proc_data(:,i,1).^2))
	% get the threshold crossings (based on first channel)

	spikes(i)=ephys_spike_detect(squeeze(proc_data(:,i,:)),spikethreshold,'fs',EPHYS.fs,'visualize','n');

	% get the spikeless data

	tmp=ephys_spike_removespikes(proc_data(:,i,1),spikes(i));

	threshold(i)=spikethreshold;
	totalspikes=totalspikes+length(spikes(i).times);

end

fprintf(1,'\n');
disp([ num2str(totalspikes) ' total spikes']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% simplest way to draw raster...plot([sample;sample]
% vector of times [120 130 200;120 130 200], then next vector specifies bottom and top of tickmark
% e.g. [1.5 1.5 1.5;.5 .5 .5] for trial 1, etc.

% collapse spike windows into a single matrix, get cluster ids back and then 
% change color according to ID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% IFR, SMOOTH RATE %%%%%%%%%%%%%%%%%%%

IFR=zeros(ntrials,samples);

for i=1:ntrials

	tmp=spikes(i).times;

	% IFR will use the current sampling rate

	IFR(i,:)=ephys_ifr(tmp,samples,EPHYS.fs);
end

