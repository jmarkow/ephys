function [SPIKES_PP SPIKES_PB]=ephys_spike_detect(DATA,THRESH,varargin)
%ephys_spike_detect.m performs spike detection on a vector with a pre-determined
%threshold
%
% [SPIKES_PP SPIKES_PB]=spike_detect(DATA,sr,traces,THRESH)
%
%
% DATA
% sample x trace matrix of voltage recrordings (threshold crossings detected on the first column, the rest are slaved)
%
% sr
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
window=[.0006 .0006]; % how large of a window to grab, seconds before and after spike
method='b';
visualize='y';
sr=25e3;
realign='y';
interpolate='y'; % do we want to interpolate for realignment (default 'y');
interpolate_fs=50e3; % what fs should we intepolate to? (50e3 has worked in my hands, consider going higher for low SNR)
align='com'; % you'll want to use COM here, others seem a bit unreliable
jitter=4; % how much jitter do we allow before tossing out a spike (in samples of original sr)?

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
		case 'sr'
			sr=varargin{i+1};
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
	if lower(method(1))=='p' || lower(method(1))=='a'
		SPIKES_PB.pos.trains=zeros(size(DATA(:,1)),'int8');
	end
	
	if lower(method(1))=='n' || lower(method(1))=='a'
		SPIKES_PB.neg.trains=zeros(size(DATA(:,1)),'int8');
	end

	if lower(method(1))=='b' || lower(method(1))=='a'
		SPIKES_PB.abs.trains=zeros(size(DATA(:,1)),'int8');
	end
end

% specify the spike window in terms of samples

spike_window=round(window.*sr);
spike_window_length=length(-spike_window(1):spike_window(2));
spike_window_center=floor(median([1:spike_window_length]));

if lower(interpolate(1))=='y'
	expansion=interpolate_fs/sr;
end

% collect the pos-going and neg-going spikes

if lower(method(1))=='p' || lower(method(1))=='a'

	% find our peaks
	[pos_times]=find(DATA(:,1)>THRESH)
	counter=1;

	SPIKES_PP.pos.times=[];
	SPIKES_PP.pos.values=[];
	SPIKES_PP.pos.windows=[];

	% if we can get a spike window around the peak, then collect for output
	counter=2;
	while counter<=length(pos_times)
		dtime=pos_times(counter)-pos_times(counter-1);
		if dtime<shadow*sr
			pos_times(counter)=[];
		else
			counter=counter+1;
		end
	end

	counter=1;
	for j=1:length(pos_times)

		if pos_times(j)-spike_window(1)>0 && pos_times(j)+spike_window(2)<length(DATA(:,1))


			% find the absolute minimum (or max) and use as the spike peak for alignment

			tmp_time=pos_times(j);
			tmp_window=DATA(tmp_time-spike_window(1):tmp_time+spike_window(2),1);

			% upsample the window with sinc interpolation

			timepoints=[1:length(tmp_window)]';
			samples=length(tmp_window);

			if lower(interpolate(1))=='y'
				newtimepoints=linspace(1,samples,expansion*samples)';

				% sinc interpolation

				%interp_window=sinc(newtimepoints(:,ones(size(timepoints)))-...
				%	timepoints(:,ones(size(newtimepoints)))')*tmp_window(:);

				interp_window=spline(timepoints,tmp_window,newtimepoints);

			else
				newtimepoints=timepoints;
				interp_window=tmp_window;
			end

			% need to add escape hatch with no interpolate

			% align by com (center of mass), min or max

			switch lower(align)

				case 'com'

					[val loc]=max(abs(interp_window));

					% take the samples around the min peak and compute the
					% COM

					% per logothetis and mahani, take continuous region > 50% of peak value

					compoints=find(abs(interp_window)>.5*val);
					
					% get the continuous region

					compoints=compoints(diff(compoints)==1);
					
					% any breakpoints between continuous regions, take the first region
					
					breakpoints=find(diff(compoints)>1);
					
					if ~isempty(breakpoints)
						compoints=compoints(1:breakpoints(1));
					end

					% toss out any points where alignment>jitter	
			
					com=sum(compoints.*abs(interp_window(compoints)))/sum(abs(interp_window(compoints)));

					if isnan(com)
						continue;
					end

					alignpoint=round(newtimepoints(round(com)));

					if abs(alignpoint-spike_window_center)>jitter
						continue;
					end
				
				case 'min'
					[val loc]=min(interp_window);
					alignpoint=round(newtimepoints(loc));
				case 'max'
					[val loc]=max(interp_window);
					alignpoint=round(newtimepoints(loc));

			end

			new_time=alignpoint;
			new_time=pos_times(j)-spike_window(1)+(new_time-1);

			if new_time-spike_window(1)>0 && new_time+spike_window(2)<length(DATA)

				SPIKES_PP.pos.times(counter)=new_time;
				SPIKES_PP.pos.values(counter)=DATA(new_time,1);
				
				for k=1:traces
					SPIKES_PP.pos.windows(:,counter,k)=DATA(new_time-spike_window(1):new_time+spike_window(2),k);
				end
				
				counter=counter+1;

			end

		end

	end

	counter=2;
	while counter<=length(SPIKES_PP.pos.times)
		dtime=SPIKES_PP.pos.times(counter)-SPIKES_PP.pos.times(counter-1);
		if dtime<shadow*sr
			SPIKES_PP.pos.times(counter)=[];
			SPIKES_PP.pos.values(counter)=[];
			SPIKES_PP.pos.windows(:,counter,:)=[];
		else
			counter=counter+1;
		end
	end

end

% same principle applies here...

if lower(method(1))=='n' || lower(method(1))=='a'

	[neg_times]=find(-DATA>THRESH);	
	counter=1;

	SPIKES_PP.neg.times=[];
	SPIKES_PP.neg.values=[];
	SPIKES_PP.neg.windows=[];

	counter=2;
	while counter<=length(neg_times)
		dtime=neg_times(counter)-neg_times(counter-1);
		if dtime<shadow*sr
			neg_times(counter)=[];
		else
			counter=counter+1;
		end
	end

	counter=1;
	for j=1:length(neg_times)

		if neg_times(j)-spike_window(1)>0 && neg_times(j)+spike_window(2)<length(DATA(:,1))


			% find the absolute minimum (or max) and use as the spike peak for alignment

			tmp_time=neg_times(j);
			tmp_window=DATA(tmp_time-spike_window(1):tmp_time+spike_window(2),1);

			% upsample the window with sinc interpolation

			timepoints=[1:length(tmp_window)]';
			samples=length(tmp_window);

			if lower(interpolate(1))=='y'
				newtimepoints=linspace(1,samples,expansion*samples)';

				% sinc interpolation

				%interp_window=sinc(newtimepoints(:,ones(size(timepoints)))-...
				%	timepoints(:,ones(size(newtimepoints)))')*tmp_window(:);

				interp_window=spline(timepoints,tmp_window,newtimepoints);

			else
				newtimepoints=timepoints;
				interp_window=tmp_window;
			end

			% need to add escape hatch with no interpolate

			% align by com (center of mass), min or max

			switch lower(align)

				case 'com'

					[val loc]=max(abs(interp_window));

					% take the samples around the min peak and compute the
					% COM

					% per logothetis and mahani, take continuous region > 50% of peak value

					compoints=find(abs(interp_window)>.5*val);
					
					% get the continuous region

					compoints=compoints(diff(compoints)==1);
					
					% any breakpoints between continuous regions, take the first region
					
					breakpoints=find(diff(compoints)>1);
					
					if ~isempty(breakpoints)
						compoints=compoints(1:breakpoints(1));
					end

					% toss out any points where alignment>jitter	
			
					com=sum(compoints.*abs(interp_window(compoints)))/sum(abs(interp_window(compoints)));
					
					if isnan(com)
						continue;
					end
					
					alignpoint=round(newtimepoints(round(com)));

					if abs(alignpoint-spike_window_center)>jitter
						continue;
					end

				case 'min'
					[val loc]=min(interp_window);
					alignpoint=round(newtimepoints(loc));
				case 'max'
					[val loc]=max(interp_window);
					alignpoint=round(newtimepoints(loc));

			end

			new_time=alignpoint;
			new_time=neg_times(j)-spike_window(1)+(new_time-1);

			if new_time-spike_window(1)>0 && new_time+spike_window(2)<length(DATA(:,1))

				SPIKES_PP.neg.times(counter)=new_time;
				SPIKES_PP.neg.values(counter)=DATA(new_time,1);

				SPIKES_PP.neg.windows(:,counter,:)=DATA(new_time-spike_window(1):new_time+spike_window(2),:);
				counter=counter+1;

			end

		end
	end

	counter=2;
	while counter<=length(SPIKES_PP.neg.times)
		dtime=SPIKES_PP.neg.times(counter)-SPIKES_PP.neg.times(counter-1);
		if dtime<shadow*sr
			SPIKES_PP.neg.times(counter)=[];
			SPIKES_PP.neg.values(counter)=[];
			SPIKES_PP.neg.windows(:,counter,:)=[];
		else
			counter=counter+1;
		end
	end


end

if lower(method(1))=='b' || lower(method(1))=='a'

	[abs_times]=find(abs(DATA)>THRESH);
	nspikes=length(abs_times);

	% censor period

	counter=2;
	while counter<=length(abs_times)
		dtime=abs_times(counter)-abs_times(counter-1);
		if dtime<shadow*sr
			abs_times(counter)=[];
		else
			counter=counter+1;
		end
	end

	SPIKES_PP.abs.times=zeros(1,nspikes);
	SPIKES_PP.abs.values=zeros(1,nspikes);
	SPIKES_PP.abs.windows=zeros(spike_window_length,nspikes,traces);

	counter=1;
	for j=1:length(abs_times)

		if abs_times(j)-spike_window(1)>0 && abs_times(j)+spike_window(2)<length(DATA(:,1))


			% find the absolute minimum (or max) and use as the spike peak for alignment

			tmp_time=abs_times(j);
			tmp_window=DATA(tmp_time-spike_window(1):tmp_time+spike_window(2),1);

			% upsample the window with sinc interpolation

			timepoints=[1:length(tmp_window)]';
			samples=length(tmp_window);

			if lower(interpolate(1))=='y'
				newtimepoints=linspace(1,samples,expansion*samples)';

				% sinc interpolation

				%interp_window=sinc(newtimepoints(:,ones(size(timepoints)))-...
				%	timepoints(:,ones(size(newtimepoints)))')*tmp_window(:);

				interp_window=spline(timepoints,tmp_window,newtimepoints);

			else
				newtimepoints=timepoints;
				interp_window=tmp_window;
			end

			% need to add escape hatch with no interpolate
			% align by com (center of mass), min or max

			switch lower(align)

				case 'com'

					[val loc]=max(abs(interp_window));

					% take the samples around the min peak and compute the
					% COM

					% per logothetis and mahani, take continuous region > 50% of peak value

					compoints=find(abs(interp_window)>.5*val);
					
					% get the continuous region

					compoints=compoints(diff(compoints)==1);
					
					% any breakpoints between continuous regions, take the first region
					
					breakpoints=find(diff(compoints)>1);
					
					if ~isempty(breakpoints)
						compoints=compoints(1:breakpoints(1));
					end

					% toss out any points where alignment>jitter	
			
					com=sum(compoints.*abs(interp_window(compoints)))/sum(abs(interp_window(compoints)));
					
					if isnan(com)
						continue;
					end

					alignpoint=round(newtimepoints(round(com)));

					if abs(alignpoint-spike_window_center)>jitter
						continue;
					end

				case 'min'
					[val loc]=min(interp_window);
					alignpoint=round(newtimepoints(loc));
				case 'max'
					[val loc]=max(interp_window);
					alignpoint=round(newtimepoints(loc));

			end

			new_time=alignpoint;
			new_time=abs_times(j)-spike_window(1)+(new_time-1);

			if new_time-spike_window(1)>0 && new_time+spike_window(2)<length(DATA(:,1))

				SPIKES_PP.abs.times(counter)=new_time;
				SPIKES_PP.abs.values(counter)=DATA(new_time,1);
				SPIKES_PP.abs.windows(1:spike_window_length,counter,1:traces)=...
					DATA(new_time-spike_window(1):new_time+spike_window(2),1:traces);

				counter=counter+1;

			end

		end
	end

	% how much of array is left unused...

	SPIKES_PP.abs.times(counter:nspikes)=[];
	SPIKES_PP.abs.values(counter:nspikes)=[];
	SPIKES_PP.abs.windows(:,counter:nspikes,:)=[];

	% make sure we haven't made any alignments that violate the censor period

	counter=2;
	while counter<=length(SPIKES_PP.abs.times)
		dtime=SPIKES_PP.abs.times(counter)-SPIKES_PP.abs.times(counter-1);
		if dtime<shadow*sr
			SPIKES_PP.abs.times(counter)=[];
			SPIKES_PP.abs.values(counter)=[];
			SPIKES_PP.abs.windows(:,counter,:)=[];
		else
			counter=counter+1;
		end
	end


end

% if the user wants binned point processes, give them to him/her!

if nargout>1
	% all of the spikes times are given as one, note that the bin size
	% here is ONE SAMPLE, hence we cannot have more than one spike in the bin

	if lower(method(1))=='p' || lower(method(1))=='a'
		SPIKES_PB.pos.trains(SPIKES_PP.pos.times)=1;
	end

	if lower(method(1))=='n' || lower(method(1))=='a'
		SPIKES_PB.neg.trains(SPIKES_PP.neg.times)=1;
	end

	if lower(method(1))=='b' || lower(method(1))=='a'
		SPIKES_PB.abs.trains(SPIKES_PP.abs.times)=1;
	end

end

% visualize the voltage trace, threshold(s) and spikes

if lower(visualize(1))=='y'

	nsamples=length(DATA);
	figure();
	plot([1:nsamples]./sr,DATA,'b');hold on
	ylabel({'Voltage (in V)';['Threshold (in V):  ' num2str(THRESH)]},'FontSize',13,'FontName','Helvetica');
	xlabel('T (in s)','FontSize',13,'FontName','Helvetica');
	plot([1:nsamples]./sr,ones(nsamples,1).*THRESH,'r');
	plot([1:nsamples]./sr,ones(nsamples,1).*-THRESH,'r');

	if lower(method(1))=='p' || lower(method(1))=='a'
		plot(SPIKES_PP.pos.times/sr,SPIKES_PP.pos.values,'r*','markersize',10);
	end

	if lower(method(1))=='n'  || lower(method(1))=='a'
		plot(SPIKES_PP.neg.times/sr,SPIKES_PP.neg.values,'k*','markersize',10);
	end

	if lower(method(1))=='b' || lower(method(1))=='a'
		plot(SPIKES_PP.abs.times/sr,SPIKES_PP.abs.values,'b*','markersize',10);
	end

	set(gca,'FontSize',11,'FontName','Helvetica')
	box off

end

