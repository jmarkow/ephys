function [SPIKECOUNT LFPWINS spikebins]=ephys_su_lfp_coherence_peaktrigger(LFPDATA,SPIKETIMES,varargin)
%computes IFR triggered LFPs
%
%	[LFPWINS_PEAK,LFPWINS_TROUGH,LFPWINS_RAND]=ephys_su_lfp_coherence_trigger(LFPCHANNEL,SUCHANNEL,SUCLUSTER,varargin)
%	
%	LFPCHANNEL
%	LFPCHANNEL to use
%
%	SUCHANNEL
%	channel with single units (already clustered and saved in 'sua')
%
%	SUCLUSTER
%	number of cluster with single unit
%
%	the following parameters may be passed as parameter/value pairs:
%
%		filedir
%		base directory with clustered single unit data in the folder "sua" and the lfp data in aggregated_data.mat
%
%		savedir
%		directory to save results (default: 'coherence')
%
%		lfp_winextract
%		length of the window about the peak and trough points to extract (in seconds, default: [.03 .03])
%		
%		peakedges
%		sets threshold for positive and negative-going threshold crossings of peakedges(1) and peakedges(2) in IFR (default: [300 300])
%
%		troughedges
%		sets threshold for negative and positive-going threshold crossings of peakedges(1) and peakedges(2) in IFR (default: [100 100])
%
%		fig_title
%		prefix used for graph filenames (default: 'noname')
%
%		debug
%		display detected peaks and troughs and pause for each trial
%
%		trial_min
%		minimum peaks/troughs detected for display (default: 20)
%
%		peak
%		how to detect peaks, use either the leading edge (e) or the mean of the positive and negative-going threshold crossings (default: 'm')
%
%		trough
%		how to detect troughs, use either the leading edge (e) or the mean of the positive and negative-going threshold crossings (default: 'm')
%		
%		medfilt_scale
%		median filter scale, in ms (default: 1.5)
%
%		lfp_fs
%		sampling rate of fields (default: 25e3)
%	
%
%

% TODO finish documentation
% TODO comment thoroughly
% TODO pare down to compute only IFR triggered fields


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETER COLLECTION  %%%%%%%%%%%%%%

if nargin<2
	error('ephysPipeline:ifrcoherence:notenoughparams','Need 2 arguments to continue, see documentation');
end

nparams=length(varargin);
filedir=pwd;
savedir=pwd;
lfp_winextract=[.1 .1]; % msec before and after peak or trough to grab LFP
peakedge=[1.5]; % take the leading edge of the peak
troughedge=[-2];
fig_title='noname';
debug=0;
freq_range=[10 50];
filt_order=4;
proc_fs=1e3;
medfilt_scale=.0015; % median filter scale (in s)
fs=25e3;
ifr_fs=25e3;
subtrials=[];

filt_type='bandpass';
xres=.001;
ifr_thresh=[];

if mod(nparams,2)>0
	error('ephysPipeline:argChk','Parameters must be specified as parameter/value pairs!');
end

for i=1:2:nparams
	switch lower(varargin{i})		
		case 'lfp_winextract'
			lfp_winextract=varargin{i+1};
		case 'debug'
			debug=varargin{i+1};
		case 'troughedge'
			troughedge=varargin{i+1};
		case 'freq_range'
			freq_range=varargin{i+1};
		case 'medfilt_scale'
			medfilt_scale=varargin{i+1};
		case 'fs'
			fs=varargin{i+1};
		case 'ifr_fs'
			ifr_fs=varargin{i+1};
		case 'proc_fs'
			proc_fs=varargin{i+1};
		case 'subtrials'
			subtrials=varargin{i+1};
		case 'xres'
			xres=varargin{i+1};
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SIGNAL CONDITIONING %%%%%%%%%%%%%%%%

downfact=fs/proc_fs;

if mod(downfact,1)>0
	error('ephysPipeline:spectcoherence:downsamplenotinteger','Need to downsample by integer');
end

% first median filter

if ~isempty(medfilt_scale)
	LFPDATA=medfilt1(double(LFPDATA),round(medfilt_scale*fs));
end

% anti-alias and downsample

[b,a]=butter(3,[300/(fs/2)],'low');
lfp_data=filtfilt(b,a,double(LFPDATA));
lfp_data=downsample(lfp_data,downfact);

clear EPHYS_DATA;

lfp_data=ephys_condition_signal(lfp_data,'l','freq_range',freq_range,'medfilt_scale',medfilt_scale,'medfilt',0,...
	'filt_order',filt_order,'filt_type',filt_type,'fs',proc_fs);
lfp_data=squeeze(lfp_data);

% need to account for subset of trials if used in single unit data 

[samples,ntrials]=size(lfp_data);

if isempty(subtrials)
	subtrials=1:ntrials;
end

lfp_data=lfp_data(:,subtrials);
[samples,ntrials]=size(lfp_data);
lfp_data=lfp_data';
lfp_time=[1:samples]./proc_fs;

indx=1:samples-1;
spikebins=-lfp_winextract(1):xres:lfp_winextract(2);

SPIKECOUNT=zeros(1,length(spikebins));
lfp_winextract=round(lfp_winextract*proc_fs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get the IFR sampling rate

% take IFR zscore from each trial, check for peaks and troughs
% extract LFPs triggered on "peaks" or "bursts" and "pauses"

% peak detection: check for zero crossings in the positive-going direction
% trough detection:  check for zero crossings in the negative-going direction

% for peaks and troughs, we can set the magnitude of change for the zero-crossings, i.e. x[i+1]>>x[i] | x[i+1]<<x[i]
% this will detect the **leading edge** of peaks and troughs, so interpret accordingly

% find the troughs and compute the spike correlogram centered at the trough
LFPWINS=[];

for i=1:ntrials

	currlfp=zscore(lfp_data(i,:));
	difflfp=diff(currlfp);
	currtimes=round(SPIKETIMES{i}*proc_fs);

	% only include burst spikes

	if ~isempty(ifr_thresh)
		currisi=[];
		tmp=[-inf currtimes inf];
		for j=2:length(tmp)-1
			nextspike=abs(tmp(j)-tmp(j+1))./proc_fs;
			prevspike=abs(tmp(j)-tmp(j-1))./proc_fs;
			currisi(j-1)=min([nextspike prevspike]);
		end
		currifr=1./currisi;
		currtimes=currtimes(currifr>ifr_thresh);
	end

	if isempty(currtimes), continue; end

	zerocross{i}=find(currlfp(indx)>troughedge & currlfp(indx+1)<troughedge)+1; 
	%randcross{i}=randsample(1:length(currlfp),length(poscross{i}));

	ncross(i)=length(zerocross);

	% grab LFP windows around the troughs

	for j=1:length(zerocross{i})	
		
		lfpcenter=zerocross{i}(j);

		% find the turning point, that's our in

		for k=lfpcenter:samples-1
			if difflfp(k)>=0
				lfpcenter=k;
				break;
			end
		end

		startidx=lfpcenter-lfp_winextract(1);
		stopidx=lfpcenter+lfp_winextract(2);

		zerocross{i}(j)=lfpcenter;

		if startidx>0 & stopidx<length(currlfp)
			tmp=ksdensity((currtimes-lfpcenter)./proc_fs,spikebins);
			if isempty(tmp), continue; end
			SPIKECOUNT=SPIKECOUNT+tmp;
			LFPWINS=[LFPWINS lfp_data(i,startidx:stopidx)'];
		end
		
	end

	% grab LFP windows around the peaks
	% show the results if debugging

	if debug
		lfpfig=figure();
		plot(lfp_time,currlfp);
		hold on;
		if ~isempty(zerocross{i})
			plot(lfp_time(zerocross{i}),currlfp(zerocross{i}),'g*');
		end
	
		pause();
		close([lfpfig]);
	end

end
