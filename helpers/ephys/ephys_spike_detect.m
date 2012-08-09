function [SPIKES_PP SPIKES_PB]=ephys_spike_detect(DATA,THRESH,varargin)
%ephys_spike_detect.m performs spike detection on a vector with a pre-determined
%threshold
%
% [SPIKES_PP SPIKES_PB]=spike_detect(DATA,fs,traces,THRESH)
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
% shadow
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
% SPIKES_PP
% array of structures (number of traces)
%
% SPIKES_PP.pos.times 
% pos going spike times (in samples) 
%
% SPIKES_PP.pos.values
% pos going spike values
%
% SPIKES_PP.pos.window
% pos going spike windows
% 
% SPIKES_PP.neg, same as pos
%
% SPIKES_PP.abs, same as pos
%
% SPIKES_PB
% array of structures (number of traces)
%
% SPIKES_PB.pos.trains
% sample x trace matrix with binned spikes (each bin is 1 sample)
%

%DATA=DATA(:);
SPIKES_PP=[];
SPIKES_PB=[];

if nargin<1
	error('Need the input data to continue!');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%

% input argument collection

shadow=.75e-3; % minimum time between spikes, i.e. censor period
	       % per Hill, Mehta and Kleinfeld (2011), 750 microseconds
window=[.0004 .0004]; % how large of a window to grab, seconds before and after spike
method='b';
visualize='y';
fs=25e3;
realign='y';
interpolate=1; % do we want to interpolate for realignment and subsequent sorting (default 'y');
interpolate_fs=50e3; % what fs should we intepolate to? (50e3 has worked in my hands, consider going higher for low SNR)
align='com'; % you'll want to use COM here, others seem a bit unreliable
jitter=4; % how much jitter do we allow before tossing out a spike (in samples of original fs)?
peak_frac=.5; % fraction of peak to use as cutoff for COM calculation (i.e. all samples below peak_frac*peak are included)
peak_width=5; % how many samples about the peak to include in COM (interpolated space, 5-8 is reasonable here for 50e3 fs)

%%%%%%%%%%%%%%%%%%%%%%%%%%

nparams=length(varargin);

if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs!');
end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'shadow'
			shadow=varargin{i+1};
		case 'window'
			window=varargin{i+1};
		case 'sigma'
			sigma=varargin{i+1};
		case 'method'
			method=varargin{i+1};
		case 'visualize'
			visualize=varargin{i+1};
		case 'fs'
			fs=varargin{i+1};
		case 'interpolate'
			interpolate=varargin{i+1};
		case 'interpolate_fs'
			interpolate_fs=varargin{i+1};
		case 'align'
			align=varargin{i+1};
		case 'jitter'
			jitter=varargin{i+1};
		case 'tetrode_data'
			tetrode_data=varargin{i+1};
	end
end

[samples,traces]=size(DATA);

%%%%%%%%%%%%%%%%%%%%%%%%%%%

% pre-allocate the binned output if the user wants it!

if nargout>1
	if lower(method(1))=='b' || lower(method(1))=='a'
		SPIKES_PB.abs.trains=zeros(size(DATA(:,1)),'int8');
	end
end

% specify the spike window in terms of samples

% the frame to grab (in the original sampling rate) around the threshold crossing

frame=round(window*fs);
frame=frame+jitter;
frame_length=length([-frame(1):frame(2)]);

if interpolate

	% interpolation factor

	expansion=interpolate_fs/fs;

	% window for interpolated spike

	spike_window=round(window*interpolate_fs);

	% the frame will be expanded in the interpolation, adjust the center accordingly

	frame_center=floor(median([1:frame_length*expansion]));
else
	expansion=1;
	spike_window=round(window*fs);
	frame_center=floor(median([1:frame_length]));
end

spike_window_length=length(-spike_window(1):spike_window(2));
spike_window_center=floor(median([1:spike_window_length]));

% collect the pos-going and neg-going spikes

[pos_times]=find(DATA>THRESH);
[neg_times]=find(DATA<-THRESH);

if lower(method(1))=='b' || lower(method(1))=='a'
	
	abs_times=[pos_times(:);neg_times(:)];
	abs_isneg=[zeros(size(pos_times));ones(size(neg_times))];
	[abs_times idx]=sort(abs_times);
	abs_isneg=abs_isneg(idx);

elseif lower(method(1))=='p'
	abs_times=pos_times;
	abs_isneg=zeros(size(pos_times));
else
	abs_times=neg_times;
	abs_isneg=ones(size(neg_times));
end

nspikes=length(abs_times);

% censor period

counter=2;
while counter<=length(abs_times)
	dtime=abs_times(counter)-abs_times(counter-1);
	if dtime<shadow*fs
		abs_times(counter)=[];
		abs_isneg(counter)=[];
	else
		counter=counter+1;
	end
end

SPIKES_PP.abs.times=zeros(1,nspikes);
SPIKES_PP.abs.values=zeros(1,nspikes);
SPIKES_PP.abs.windows=zeros(spike_window_length,nspikes,traces);

counter=1;
for j=1:length(abs_times)

	if abs_times(j)-frame(1)>0 && abs_times(j)+frame(2)<length(DATA(:,1))

		% find the absolute minimum (or max) and use as the spike peak for alignment

		tmp_time=abs_times(j);
		isneg=abs_isneg(j);
		tmp_window=DATA(tmp_time-frame(1):tmp_time+frame(2),:);

		% find the absolute min in the window

		[val loc]=min(tmp_window(:,1));
		peak_time=tmp_time-frame(1)+(loc(1)-1);

		% need to grab new time based on absolute peak, grab extra samples for jitter

		if peak_time-frame(1)>0 && peak_time+frame(2)<length(DATA(:,1))
			tmp_window=DATA(peak_time-frame(1):peak_time+frame(2),:);
		else
			continue;
		end
	
		% upsample the window with sinc interpolation, or spline

		[samples,channels]=size(tmp_window);

		timepoints=[1:samples];

		if interpolate	

			newtimepoints=linspace(1,samples,expansion*samples)';

			% sinc interpolation

			%interp_window=sinc(newtimepoints(:,ones(size(timepoints)))-...
			%	timepoints(:,ones(size(newtimepoints)))')*tmp_window(:);

			% spline interpolation

			for k=1:channels
				interp_window(:,k)=spline(timepoints,tmp_window(:,k),newtimepoints);
			end

		else

			newtimepoints=timepoints;
			
			for k=1:channels
				interp_window(:,k)=tmp_window(:,k);
			end

		end

		% align by com (center of mass), min or max
		% first realign to min or max, whichever is more reliable

		switch lower(align)

			case 'com'

				% the negative-going peak has been the most reliable, use it to compute COM
				% per sahani '99, try to capture the majority of the peak

				[val loc]=min(interp_window(:,1));
				compointsneg=find(interp_window(:,1)<=peak_frac*val);
				compoints=sort(compointsneg);

				% take only points about the peak

				%compoints=compoints(diff(compoints)==1);
				compoints((compoints>loc+(peak_width))|(compoints<(loc-peak_width)))=[];
			
				% toss out any points where alignment>jitter, assume to be outliers

				com=sum(compoints.*abs(interp_window(compoints,1)))/sum(abs(interp_window(compoints,1)));

				if isnan(com) || isempty(com)
					continue;
				end

				alignpoint=round(com);
				
				% if the alignpoint is too far from the center of the frame, dump the spike (likely outlier)

				if abs(alignpoint-frame_center)>=jitter*expansion
					continue;
				end

			
			case 'min'

				% just take the min

				[val loc]=min(interp_window(:,1));

				if abs(loc-frame_center)>=jitter*expansion
					continue;
				end

				alignpoint=loc;

			case 'max'

				% just take the max

				[val loc]=max(interp_window(:,1));

				if abs(loc-frame_center)>=jitter*expansion
					continue;
				end

				alignpoint=loc;

		end

		% new spike window

		new_spikewindow=interp_window(alignpoint-spike_window(1):alignpoint+spike_window(2),:);

		% get the spike time in the old sample space

		new_time=peak_time-frame(1)+(round(newtimepoints(alignpoint))-1);
		new_val=DATA(new_time,1);
		%new_time=peak_time-spike_window(1)+(alignpoint-1);

		SPIKES_PP.abs.times(counter)=new_time;
		SPIKES_PP.abs.values(counter)=new_val;
		SPIKES_PP.abs.windows(1:spike_window_length,counter,1:traces)=new_spikewindow;
		counter=counter+1;


	end
end

% how much of array is left unused...

SPIKES_PP.abs.times(counter:nspikes)=[];
SPIKES_PP.abs.values(counter:nspikes)=[];
SPIKES_PP.abs.windows(:,counter:nspikes,:)=[];
if interpolate
	SPIKES_PP.abs.fs=interpolate_fs;
else
	SPIKES_PP.abs.fs=fs;
end

% make sure we haven't made any alignments that violate the censor period

counter=2;
while counter<=length(SPIKES_PP.abs.times)
	dtime=SPIKES_PP.abs.times(counter)-SPIKES_PP.abs.times(counter-1);
	if dtime<shadow*fs
		SPIKES_PP.abs.times(counter)=[];
		SPIKES_PP.abs.values(counter)=[];
		SPIKES_PP.abs.windows(:,counter,:)=[];
	else
		counter=counter+1;
	end
end


% if the user wants binned point processes, give them to him/her!

if nargout>1
	% all of the spikes times are given as one, note that the bin size
	% here is ONE SAMPLE, hence we cannot have more than one spike in the bin

	if lower(method(1))=='b' || lower(method(1))=='a'
		SPIKES_PB.abs.trains(SPIKES_PP.abs.times)=1;
	end

end

% visualize the voltage trace, threshold(s) and spikes

if lower(visualize(1))=='y'

	nsamples=length(DATA);
	figure();
	plot([1:nsamples]./fs,DATA,'b');hold on
	ylabel({'Voltage (in V)';['Threshold (in V):  ' num2str(THRESH)]},'FontSize',13,'FontName','Helvetica');
	xlabel('T (in s)','FontSize',13,'FontName','Helvetica');
	plot([1:nsamples]./fs,ones(nsamples,1).*THRESH,'r');
	plot([1:nsamples]./fs,ones(nsamples,1).*-THRESH,'r');

	plot(SPIKES_PP.abs.times/fs,SPIKES_PP.abs.values,'b*','markersize',10);
	set(gca,'FontSize',11,'FontName','Helvetica')
	box off

end

