function [spikeamp,spikephase]=ephys_su_lfp_spikeprob(LFPDATA,SPIKETIMES,SPIKETRIALS,varargin)
%computes SPIKE triggered LFP STATS
%
%	[LFPWINS_PEAK,LFPWINS_TROUGH,LFPWINS_RAND]=ephys_su_lfp_coherence_ifr(LFPCHANNEL,SUCHANNEL,SUCLUSTER,varargin)
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

lfp_winextract=[.3 .3]; % msec before and after peak or trough to grab LFP
fig_title='noname';
debug=0;
freq_range=[18 25];
filt_order=5;
proc_fs=1e3;
medfilt_scale=.0015; % median filter scale (in s)
fs=25e3;
ifr_fs=25e3;
subtrials=[];
smoothing=[]; % if ~isempty, sgolay using this span (quartic)
boundary=[.4 .4]; % do not take any spikes 300 msec from the edge
ifr_cutoff=[];

if mod(nparams,2)>0
	error('ephysPipeline:argChk','Parameters must be specified as parameter/value pairs!');
end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'filedir'
			filedir=varargin{i+1};
		case 'savedir'
			savedir=varargin{i+1};
		case 'lfp_winextract'
			lfp_winextract=varargin{i+1};
		case 'debug'
			debug=varargin{i+1};
		case 'peakedges'
			peakedges=varargin{i+1};
		case 'troughedges'
			troughedges=varargin{i+1};
		case 'fig_title'
			fig_title=varargin{i+1};
		case 'freq_range'
			freq_range=varargin{i+1};
		case 'peak'
			peak=varargin{i+1};
		case 'trough'
			trough=varargin{i+1};
		case 'medfilt_scale'
			medfilt_scale=varargin{i+1};
		case 'fs'
			fs=varargin{i+1};
		case 'proc_fs'
			proc_fs=varargin{i+1};
		case 'subtrials'
			subtrials=varargin{i+1};
		case 'ifr_cutoff'
			ifr_cutoff=varargin{i+1};
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
lfpdata=filtfilt(b,a,double(LFPDATA));
lfpdata=downsample(lfpdata,downfact);

clear EPHYS_DATA;

[samples,ntrials]=size(lfpdata);

if isempty(subtrials)
	subtrials=1:ntrials;
end

lfpdata=lfpdata(:,subtrials);
[samples,ntrials]=size(lfpdata);
%lfpdata=lfpdata';
lfp_time=[1:samples]./proc_fs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% since we'll have so many more trials here, we may need to 
% set a trial limit

lfp_winextract=round(lfp_winextract.*proc_fs);
winlength=length([-lfp_winextract(1):lfp_winextract(2)]);

% need to account for subset of trials if used in single unit data 

%lfpdata=lfpdata(subtrials,:);
% vector for peak/trough detection

% count the number of extractions in each case to pre-allocate (major speedup)
% bank of filters at the relevant time-scales

% set up filters...

trials=unique(SPIKETRIALS);
spikephase=zeros(1,length(SPIKETRIALS));
spikeamp=zeros(size(spikephase));

[b,a]=ellip(4,.1,40,[freq_range]/(proc_fs/2));
lfpdata=filtfilt(b,a,lfpdata);

% take the hilbert transform

lfpdata=hilbert(lfpdata);

counter=0;

boundary(2)=(samples./proc_fs)-boundary(2)

for i=1:length(trials)

	currtrial=trials(i);
	
	% convert spike times to seconds

	currtimes=SPIKETIMES(SPIKETRIALS==currtrial)./fs;
	
	todel=(currtimes<=boundary(1))|currtimes>=boundary(2);

	% if there is an IFR cutoff, set all spikes < IFR cutoff to NaN

	if ~isempty(ifr_cutoff)

		currifr=zeros(1,length(currtimes));
		tmp=[-inf currtimes inf];
		for j=2:length(tmp)-1
			prev_spike=abs(tmp(j)-tmp(j-1));
			next_spike=abs(tmp(j+1)-tmp(j));
			currifr(j-1)=1./min([prev_spike next_spike]);
		end

		currifr
		todel=todel|(currifr<ifr_cutoff);

	end

	% convert to LFP samples

	currtimes=round(currtimes*proc_fs)

	% filter LFP and grab phase and power

	currlfp=lfpdata(:,currtrial);
	lfpsamples=currlfp(currtimes);

	lfpphases=angle(lfpsamples);
	lfpamps=abs(lfpsamples);

	lfpphases(todel)=NaN;
	lfpamps(todel)=NaN;

	spikephase(counter+1:counter+length(lfpsamples))=lfpphases;
	spikeamp(counter+1:counter+length(lfpsamples))=lfpamps;

	% if the time is outside the allowed range set to NaN


	% take the current spike times
	% grab the instantaneous phase at each spike along with power

	counter=counter+length(lfpsamples);

end

% now pre-allocate wave matrices

