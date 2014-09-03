function [MEANVAL,T]=ephys_bandpower(EPHYS,varargin)
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
freq_range=[1 15;20 40;40 60;60 100;100 200;200 400;500 1e3];
subtrials=[];
channels=EPHYS.labels;

car_trim=40;
winsize=4096;
winoverlap=4096-244;
winfft=4096;
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
		case 'subtrials'
			subtrials=varargin{i+1};
		case 'method'
			method=varargin{i+1};
		case 'car_trim'
			car_trim=varargin{i+1};
		case 'winsize'
			winsize=varargin{i+1};
		case 'winoverlap'
			winoverlap=varargin{i+1};
		case 'winfft'
			winfft=varargin{i+1};
	end
end

[samples,ntrials,nchannels]=size(EPHYS.data);

if isempty(subtrials)
	subtrials=1:ntrials;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SIGNAL CONDITIONING %%%%%%%%%%%%%%%%

proc_data=ephys_denoise_signal(EPHYS.data,EPHYS.labels,channels,'method',noise,'car_exclude',car_exclude,'car_trim',car_trim);
clear EPHYS.data;

proc_data=proc_data(:,subtrials,:);
[samples,ntrials,newchannels]=size(proc_data);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SPIKE DETECTION %%%%%%%%%%%%%%%%%%%%

totalspikes=0;

[nblanks formatstring]=fb_progressbar(100);
fprintf(1,['Progress:  ' blanks(nblanks)]);

nbands=size(freq_range,1);

[t,f]=getspecgram_dim(samples,winsize,winoverlap,winfft,EPHYS.fs);
MEANVAL=zeros(length(t),ntrials,nbands);

for i=1:ntrials

	fprintf(1,formatstring,round((i/ntrials)*100));	
	[s,f,t]=spectrogram(proc_data(:,i,1),winsize,winoverlap,winfft,EPHYS.fs);

	% get frequency edges
	
	for j=1:nbands

		minf=max(find(f<=freq_range(j,1)));
		maxf=min(find(f>=freq_range(j,2)));

		if isempty(minf), minf=1; end
		if isempty(maxf), maxf=length(f); end

		MEANVAL(:,i,j)=mean(s(minf:maxf,:));

	end
end

T=t;
fprintf(1,'\n');
