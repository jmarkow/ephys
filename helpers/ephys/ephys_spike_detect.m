function [SPIKES]=ephys_spike_detect(DATA,THRESH,varargin)
%ephys_spike_detect.m performs spike detection on a vector with a pre-determined
%threshold
%
% [SPIKES SPIKES_PB]=spike_detect(DATA,fs,traces,THRESH)
%
%
% DATA
% sample x trace matrix of voltage recrordings (threshold crossings detected on the first column, the rest are slaved)
%
% fs
% sampling rate of the recording
%
% THRESH
% threshold for detecting spikes
%
% the following can be specified as a parameter/value pair
%
% censor
% length between spikes in S, i.e. the censor period (default .001)
%
% window
% two element vector specifying the distance before and after the spike to store
% (default [.001 .001])
%
% method
% string specifying whether to use (p)ositive threshold crossings, (n), or (b) 
%
% visualize
% generate a figure to visualize spike detection results
%
% realign
% after spike detection, realign ('y' or 'n', defeault: 'y') to peak
%
% align
% alignment method (only applicable if realign=='y') ('min' for absolute minimum,
% 'max' for absolute maximum and 'com' for center of mass about the minimum)
%
% OUTPUT
%
% SPIKES
% array of structures (number of traces)
%
% SPIKES.pos.times 
% pos going spike times (in samples) 
%
% SPIKES.pos.values
% pos going spike values
%
% SPIKES.pos.window
% pos going spike windows
% 
% SPIKES.neg, same as pos
%
% SPIKES, same as pos
%
% SPIKES_PB
% array of structures (number of traces)
%
% SPIKES_PB.pos.trains
% sample x trace matrix with binned spikes (each bin is 1 sample)
%

%DATA=DATA(:);
SPIKES=[];

if nargin<1
	error('Need the input data to continue!');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%

% input argument collection

censor=.75e-3; % minimum time between spikes, i.e. censor period
	       % per Hill, Mehta and Kleinfeld (2011), 750 microseconds
window=[.0004 .0004]; % how large of a window to grab, seconds before and after spike
method='b';
visualize='y';
fs=25e3;
realign='y';
jitter=4; % how much jitter do we allow before tossing out a spike (in samples of original fs)?
align='min';

%%%%%%%%%%%%%%%%%%%%%%%%%%

nparams=length(varargin);

if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs!');
end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'censor'
			censor=varargin{i+1};
		case 'window'
			window=varargin{i+1};
		case 'method'
			method=varargin{i+1};
		case 'visualize'
			visualize=varargin{i+1};
		case 'fs'
			fs=varargin{i+1};
		case 'jitter'
			jitter=varargin{i+1};	
		case 'align'
			align=varargin{i+1};
	end
end

[samples,traces,channels]=size(DATA);
if samples==1
	warning('Either data is in wrong format or only contains 1 sample!');
end

% specify the spike window in terms of samples

% the frame to grab (in the original sampling rate) around the threshold crossing
% collect jitter samples around the frame for upsample and realignment

SPIKES.frame=window;
SPIKES.jitter=jitter;

frame=round(window*fs);
frame=frame+jitter;
frame_length=length([-frame(1):frame(2)]);
timepoints=-frame(1):frame(2);

% collect the pos-going and neg-going spikes

if lower(method(1))=='b' || lower(method(1))=='a'	
	spike_times=find(abs(DATA(:,1))>THRESH);
elseif lower(method(1))=='p'
	spike_times=find(DATA(:,1)>THRESH);
else
	spike_times=find(DATA(:,1)<-THRESH);
end

nspikes=length(spike_times);

% censor period

counter=2;
while counter<=length(spike_times)
	dtime=spike_times(counter)-spike_times(counter-1);
	if dtime<censor*fs
		spike_times(counter)=[];
	else
		counter=counter+1;
	end
end

SPIKES.times=zeros(1,nspikes);
SPIKES.windows=zeros(frame_length,nspikes,traces);
SPIKES.storewindows=[];
SPIKES.storetimes=[];

counter=1;

for j=1:length(spike_times)

	if spike_times(j)-frame(1)>0 && spike_times(j)+frame(2)<length(DATA(:,1))

		% find the absolute minimum (or max) and use as the spike peak for alignment

		tmp_time=spike_times(j);
		tmp_window=DATA(tmp_time-frame(1):tmp_time+frame(2),:);

		% find the absolute min in the window

		switch lower(align)
			case 'max'
				[val loc]=max(tmp_window(:,1));
			otherwise
				[val loc]=min(tmp_window(:,1));
		end

		peak_time=tmp_time-frame(1)+(loc(1)-1);

		% need to grab new time based on absolute peak, grab extra samples for jitter

		if peak_time-frame(1)>0 && peak_time+frame(2)<length(DATA(:,1))
			tmp_window=DATA(peak_time-frame(1):peak_time+frame(2),:);
		else
			continue;
		end

		SPIKES.times(counter)=peak_time;
		SPIKES.windows(:,counter,:)=tmp_window;

		counter=counter+1;

	end
end

% how much of array is left unused...

SPIKES.times(counter:nspikes)=[];
SPIKES.windows(:,counter:nspikes,:)=[];
SPIKES.fs=fs;
SPIKES.censor=censor;
SPIKES.window_time=timepoints./fs;
% make sure we haven't made any alignments that violate the censor period

counter=2;
while counter<=length(SPIKES.times)
	dtime=SPIKES.times(counter)-SPIKES.times(counter-1);
	if dtime<censor*fs
		SPIKES.times(counter)=[];
		SPIKES.windows(:,counter,:)=[];
	else
		counter=counter+1;
	end
end


% visualize the voltage trace, threshold(s) and spikes

if lower(visualize(1))=='y'

	nsamples=length(DATA);
	figure();
	plot([1:nsamples]./fs,DATA(:,1),'b');hold on
	ylabel({'Voltage (in V)';['Threshold (in V):  ' num2str(THRESH)]},'FontSize',13,'FontName','Helvetica');
	xlabel('T (in s)','FontSize',13,'FontName','Helvetica');
	plot([1:nsamples]./fs,ones(nsamples,1).*THRESH,'r');
	plot([1:nsamples]./fs,ones(nsamples,1).*-THRESH,'r');

	plot(SPIKES.times/fs,DATA(SPIKES.times,1),'b*','markersize',10);
	set(gca,'FontSize',11,'FontName','Helvetica')
	box off
	axis tight;

end

