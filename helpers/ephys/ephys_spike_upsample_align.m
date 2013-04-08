function NEWSPIKES=ephys_spike_upsample_align(SPIKES,varargin)
%
%
%
%
%

visualize='y';
interpolate_fs=200e3; % what fs should we intepolate to? (50e3 has worked in my hands, consider going higher for low SNR)
align='min'; % you'll want to use COM here, others seem a bit unreliable
peak_frac=.6; % fraction of peak to use as cutoff for COM calculation (i.e. all samples below peak_frac*peak are included)
peak_width=4; % how many samples about the peak to include in COM (interpolated space, 5-8 is reasonable here for 50e3 fs)

%%%%%%%%%%%%%%%%%%%%%%%%%%

nparams=length(varargin);

if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs!');
end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'censor'
			censor=varargin{i+1};
		case 'method'
			method=varargin{i+1};
		case 'visualize'
			visualize=varargin{i+1};
		case 'interpolate_fs'
			interpolate_fs=varargin{i+1};
		case 'align'
			align=varargin{i+1};
		case 'tetrode_data'
			tetrode_data=varargin{i+1};
	end
end

% interpolation factor


% work out the interpolation window including jitter

% read out parameters from spikes

fs=SPIKES.fs;
window=SPIKES.frame;
jitter=SPIKES.jitter;
censor=SPIKES.censor;

expansion=interpolate_fs/fs;

frame=round(window*fs);
frame=frame+jitter;
frame_length=length([-frame(1):frame(2)]);
timepoints=-frame(1):frame(2);

interp_samples=frame_length*expansion;
newtimepoints=linspace(-frame(1),frame(2),...
	interp_samples); % points in time relative to center
newframepoints=linspace(1,frame_length,interp_samples); % points in time in samples (original window)

% window for interpolated spike

frame_center=max(find(newtimepoints<=0));
spike_window=round(window*interpolate_fs); % use the window without the pad for jitter
spike_window_length=length(-spike_window(1):spike_window(2));

% number of spikes

nspikes=length(SPIKES.times);
traces=size(SPIKES.windows,3);

NEWSPIKES.times=zeros(1,nspikes);
NEWSPIKES.windows=zeros(spike_window_length,nspikes,traces);

if isfield(SPIKES,'oldwindows')
	NEWSPIKES.oldwindows=zeros(size(NEWSPIKES.windows));
end

NEWSPIKES.window_time=[-spike_window(1):spike_window(2)]'./interpolate_fs;
NEWSPIKES.frame=SPIKES.frame;
NEWSPIKES.jitter=SPIKES.jitter;
NEWSPIKES.fs=interpolate_fs;
NEWSPIKES.original_fs=fs;
NEWSPIKES.censor=censor;
counter=1;

jitter=jitter*expansion;

for i=1:nspikes

	tmp_window=SPIKES.windows(:,i,:);
	peak_time=SPIKES.times(i);
	[samples,channels]=size(tmp_window);

	% spline interpolation

	for j=1:channels
		interp_window(:,j)=spline(timepoints,tmp_window(:,j),newtimepoints);	
	end

	if isfield(SPIKES,'oldwindows')
		for j=1:channels
			interp_window2(:,j)=spline(timepoints,SPIKES.oldwindows(:,i,j),newtimepoints);
		end
	else
		interp_window2=[];
	end

	% align by com (center of mass), min or max
	% first realign to min or max, whichever is more reliable

	switch lower(align)

		case 'com'

			% the negative-going peak has been the most reliable, use it to compute COM
			% per sahani '99, try to capture the majority of the peak

			% grab the peak after the threshold crossings	
			% first get contiguous region

			[val loc]=min(interp_window(:,1));

			% we can make these for loops go away, but honestly we're not taking a huge
			% performance hit ATM...

			% walk to and fro to the threshold, use this for the peak CoM measurement

			for j=loc:-1:1
				if interp_window(j,1)>-THRESH
					break;
				end
			end

			left_edge=k;

			for j=loc:length(interp_window)
				if interp_window(j,1)>-THRESH
					break;
				end
			end

			right_edge=k;

			peakind=left_edge:right_edge;

			masked_spike=THRESH-interp_window(:,1);

			idx=[1:length(masked_spike)]';

			mask=zeros(size(idx));
			mask(peakind)=1;

			masked_spike=masked_spike.*mask;

			% toss out any points where alignment>jitter, assume to be outliers

			com=sum(idx.*masked_spike)/sum(masked_spike);

			if isnan(com) || isempty(com)
				continue;
			end

			alignpoint=round(com);

			% if the alignpoint is too far from the center of the frame, dump the spike (likely outlier)

			if abs(alignpoint-frame_center)>jitter
				continue;
			end


		case 'min'

			% just take the min

			[val loc]=min(interp_window(:,1));

			alignpoint=loc;

			%abs(alignpoint-frame_center)

			if abs(alignpoint-frame_center)>jitter
				continue;
			end

		case 'max'

			% just take the max

			[val loc]=max(interp_window(:,1));

			alignpoint=loc;

			if abs(alignpoint-frame_center)>jitter
				continue;
			end


	end

	% new spike window

	new_spikewindow=interp_window(alignpoint-spike_window(1):alignpoint+spike_window(2),:);

	% get the spike time in the old sample space

	new_time=peak_time-frame(1)+(round(newframepoints(alignpoint))-1);

	if new_time<0 
		warning('ephysPipeline:spikedetect:spiketimerror',...
			'New spike time outside data vector');
	end

	NEWSPIKES.times(counter)=new_time;
	NEWSPIKES.windows(1:spike_window_length,counter,1:traces)=new_spikewindow;

	if isfield(SPIKES,'oldwindows') 
		new_spikewindow2=interp_window2(alignpoint-spike_window(1):alignpoint+spike_window(2),:);
		NEWSPIKES.oldwindows(1:spike_window_length,counter,1:traces)=new_spikewindow2;
	end

	counter=counter+1;


end

NEWSPIKES.times(counter:nspikes)=[];
NEWSPIKES.windows(:,counter:nspikes,:)=[];

if isfield(NEWSPIKES,'oldwindows')
	NEWSPIKES.oldwindows(:,counter:nspikes,:)=[];
end

counter=2;
while counter<=length(NEWSPIKES.times)
	dtime=NEWSPIKES.times(counter)-NEWSPIKES.times(counter-1);
	if dtime<censor*fs
		NEWSPIKES.times(counter)=[];
		NEWSPIKES.windows(:,counter,:)=[];
		if isfield(NEWSPIKES,'oldwindows')
			NEWSPIKES.oldwindows(:,counter,:)=[];
		end
	else
		counter=counter+1;
	end
end


