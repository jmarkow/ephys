function [LFPWINS_PEAK LFPWINS_TROUGH LFPWINS_RAND TIMEPOINTS]=ephys_su_lfp_coherence_ifr(LFPDATA,IFRDATA,varargin)
%computes IFR triggered LFPs
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
peakedges=[50 50]; % take the leading edge of the peak
troughedges=[50 50];
fig_title='noname';
randreps=100;
singletrialplots=30; % trials chosen at random to plot with fields, spikes, and IFR
debug=0;
freq_range=[3 90];
filt_order=5;
proc_fs=1e3;
medfilt_scale=.0015; % median filter scale (in s)
fs=25e3;
hist_facecolor='none';
hist_edgecolor='m';
ifr_fs=25e3;
subtrials=[];
filt_type='bandpass';
centroid=1;
smoothing=[]; % if ~isempty, sgolay using this span (quartic)

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
		case 'ifr_fs'
			ifr_fs=varargin{i+1};
		case 'proc_fs'
			proc_fs=varargin{i+1};
		case 'hist_facecolor'
			hist_facecolor=varargin{i+1};
		case 'hist_edgecolor'
			hist_edgecolor=varargin{i+1};
		case 'subtrials'
			subtrials=varargin{i+1};
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

%lfp_data=ephys_condition_signal(lfp_data,'l','freq_range',freq_range,'medfilt_scale',medfilt_scale,'medfilt',0,...
%	'filt_order',filt_order,'filt_type',filt_type,'fs',proc_fs);
%lfp_data=squeeze(lfp_data);


% need to account for subset of trials if used in single unit data 

[samples,ntrials]=size(lfp_data);

if isempty(subtrials)
	subtrials=1:ntrials;
end

lfp_data=lfp_data(:,subtrials);
[samples,ntrials]=size(lfp_data);
lfp_data=lfp_data';
lfp_time=[1:samples]./proc_fs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get the IFR sampling rate

lfp_winextract=round(lfp_winextract.*proc_fs);
winlength=length([-lfp_winextract(1):lfp_winextract(2)]);
TIMEPOINTS=[-lfp_winextract(1):lfp_winextract(2)]./proc_fs;
% take IFR zscore from each trial, check for peaks and troughs
% extract LFPs triggered on "peaks" or "bursts" and "pauses"

% peak detection: check for zero crossings in the positive-going direction
% trough detection:  check for zero crossings in the negative-going direction

% for peaks and troughs, we can set the magnitude of change for the zero-crossings, i.e. x[i+1]>>x[i] | x[i+1]<<x[i]
% this will detect the **leading edge** of peaks and troughs, so interpret accordingly


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PREALLOCATION FOR WINDOW EXTRACTION %


[ifrtrials,ifrsamples]=size(IFRDATA)
ifr_time=[1:ifrsamples]./ifr_fs;

% need to account for subset of trials if used in single unit data 

%lfp_data=lfp_data(subtrials,:);
% vector for peak/trough detection

indx=[1:ifrsamples-1];
LFPWINS_TROUGH.waveforms=[];
LFPWINS_PEAK.waveforms=[];
LFPWINS_RAND.waveforms=[];

% count the number of extractions in each case to pre-allocate (major speedup)

zerocount=0;
poscount=0;
randcount=0;

% random shift, constant across trials for 

randshift=randi(ifrsamples,1,1);

for i=1:ntrials

	if ~isempty(smoothing)
		currifr=smooth(IFRDATA(i,:),smoothing,'sgolay',4);
	else
		currifr=IFRDATA(i,:);
	end

	currlfp=lfp_data(i,:);

	zerocross{i}=find(currifr(indx)>troughedges(1) & currifr(indx+1)<troughedges(2))+1; 
	poscross{i}=find(currifr(indx)<peakedges(1) & currifr(indx+1)>peakedges(2))+1;
	
	% wrap the random shift if necessary

	randcross{i}=mod(poscross{i}+randshift,ifrsamples);
	%randcross{i}=randsample(1:length(currifr),length(poscross{i}));

	ncross(i)=length(zerocross);

	% grab LFP windows around the troughs

	for j=1:length(zerocross{i})	
		
		lfpcenter=round((zerocross{i}(j)/ifr_fs)*proc_fs);

		% attempt to find the end of the trough

		startidx=lfpcenter-lfp_winextract(1);
		stopidx=lfpcenter+lfp_winextract(2);

		if startidx>0 & stopidx<length(currlfp)
			zerocount=zerocount+1;
		end

	end

	% grab LFP windows around the peaks

	for j=1:length(poscross{i})

		% take all contiguous bursting points

		if centroid

			edge=length(currifr);
			
			for k=poscross{i}(j)+1:length(currifr)
				if currifr(k)<peakedges(2)
					edge=k-1;
					break;
				end
			end

			poscross{i}(j)=round(mean([poscross{i}(j):edge]));
			
		end

		lfpcenter=round((poscross{i}(j)/ifr_fs)*proc_fs);
		startidx=lfpcenter-lfp_winextract(1);
		stopidx=lfpcenter+lfp_winextract(2);

		if startidx>0 & stopidx<length(currlfp)
			poscount=poscount+1;
		end

	end

	for j=1:length(randcross{i})

		lfpcenter=round((randcross{i}(j)/ifr_fs)*proc_fs);
		startidx=lfpcenter-lfp_winextract(1);
		stopidx=lfpcenter+lfp_winextract(2);

		if startidx>0 & stopidx<length(currlfp)
			randcount=randcount+1;
		end

	end

	% show the results if debugging

	if debug
		lfpfig=figure();
		plotyy(ifr_time,currifr,lfp_time,currlfp);
		hold on;
		if ~isempty(zerocross{i})
			plot(ifr_time(zerocross{i}),currifr(zerocross{i}),'g*');
		end

		if ~isempty(poscross{i})
			plot(ifr_time(poscross{i}),currifr(poscross{i}),'r*');
		end
		if ~isempty(randcross{i})
			plot(ifr_time(randcross{i}),currifr(randcross{i}),'m*');
		end
		pause();
		close([lfpfig]);
	end

end

% now pre-allocate wave matrices

LFPWINS_TROUGH.waveforms=zeros(winlength,zerocount,'single');
LFPWINS_PEAK.waveforms=zeros(winlength,poscount,'single');

zerocount=1;
poscount=1;
randcount=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WINDOW EXTRACTION %%%%%%%%%%%%%%%%%%


for i=1:ntrials

	currlfp=lfp_data(i,:);
	
	% find zero crossings where IFR > 3*Std and next point < 0
	% grab LFP windows around the troughs

	for j=1:length(zerocross{i})

		lfpcenter=round((zerocross{i}(j)/ifr_fs)*proc_fs);

		startidx=lfpcenter-lfp_winextract(1);
		stopidx=lfpcenter+lfp_winextract(2);

		if startidx>0 & stopidx<length(currlfp)
			LFPWINS_TROUGH.waveforms(1:winlength,zerocount)=currlfp(startidx:stopidx)';
			zerocount=zerocount+1;
		end

	end

	% grab LFP windows around the peaks

	for j=1:length(poscross{i})

		lfpcenter=round((poscross{i}(j)/ifr_fs)*proc_fs);

		startidx=lfpcenter-lfp_winextract(1);
		stopidx=lfpcenter+lfp_winextract(2);

		if startidx>0 & stopidx<length(currlfp)
			LFPWINS_PEAK.waveforms(1:winlength,poscount)=currlfp(startidx:stopidx)';
			poscount=poscount+1;
		end

	end

	for j=1:length(randcross{i})

		lfpcenter=round((randcross{i}(j)/ifr_fs)*proc_fs);

		startidx=lfpcenter-lfp_winextract(1);
		stopidx=lfpcenter+lfp_winextract(2);

		if startidx>0 & stopidx<length(currlfp)
			LFPWINS_RAND.waveforms(1:winlength,randcount)=currlfp(startidx:stopidx)';
			randcount=randcount+1;
		end

	end

end

