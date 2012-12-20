function [abscoh,t,f]=ephys_su_lfp_coherence_gaborhist(LFPDATA,SPIKETIMES,HISTOGRAM,varargin)
% THIS FUNCTION IS CURRENTLY INCOMPLETE
%computes coherograms between LFPs and single units
%
%	[abscoh,t,f]=ephys_su_lfp_coherence_tf(lfpdata,spiketimes,HISTOGRAM,varargin)
%	
%	lfpdata
%	matrix with lfp data (samples x trials)
%
%	spiketimes
%	cell array with spike times (in seconds)
%
%	HISTOGRAM
%	time frequency histogram structure (result of ephys_visual_histogram or loading histogram.mat)
%
%	the following may be passed as parameter/value pairs:
%
%		savedir
%		directory to save results in (default: 'coherence')
%
%		filedir
%		directory that contains the single-unit data and LFP subdirectory (default: 'pwd')
%
%		colors
%		colormap for coherence graph (default: 'jet')
%		
%		n
%		spectrogram window length (default: 6.25e3 samples)
%		
%		nfft
%		spectrogram nfft (default: 10e3 samples)
%
%		overlap
%		spectrogram overlap (default: 6e3 samples)
%		
%		min_f
%		minimum frequency to display (default: 10)
%		
%		max_f
%		maximum frequency to display (default: 100)
%
%		w
%		time-bandwidth product (default: 3)
%
%		ntapers
%		number of tapers to use (default: 2*w-1)
%
%		freq_range
%		frequency range for LFP filtering (default: 300)
%
%		lfp_fs
%		lfp sampling rate (default: 25e3)
%
%		proc_fs
%		processing sampling rate (default: 1e3, downsamples signal after filtering)
%
%		trial_range
%		only compute the coherence for these trials (leave blank for all, default: [])
%
%		medfilt_scale
%		median filter timescale in ms (default: 1.5)
%
%		nphasebins
%		number of phase bins for amplitude-weighted frequency/phase 2D histogram (default: 50)
%
%		lfp_channel
%		channel used for LFP processing
%
%		spike_channel
%		channel used for spike processing
%
%		spike_cluster
%		cluster used for spike processing
%		
%
% see also ephys_visual_histogra.m,ephys_su_lfp_coherence_spect.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETER COLLECTION %%%%%%%%%%%%%%%%%

if nargin<3
	error('ephysPipeline:tfcoherence:notenoughparams','Need four arguments to continue, see documentation...');
end

nparams=length(varargin);

savedir=pwd;
filedir=pwd;
fig_title='noname';
colors='jet';
phasecolors='hsv';
debug=0;
nfft=512;
n=512;
overlap=511;
min_f=1;
max_f=100;
w=2;
alpha=[.001 .01 .05];
smoothplot=0;
proc_fs=500; % processing frequency, will downsample lfp
window_sig=.15;

ntapers=[];
beta=23/20; % parameter for z-transforming coherence
freq_range=[5 125]; % let's just refilter
filt_order=7;
lfp_fs=25e3; % default Intan sampling rate
trial_range=[];
medfilt_scale=1.5; % median filter scale (in ms)
nphasebins=50;

clim=[];

spike_channel='null';
spike_cluster='null';
lfp_channel='null';

subtrials=[];


angles=-pi/4:pi/16:pi/4; % for contour image
wsigma=.05:.02:.15; %timescales in milliseconds
padding=[];

if mod(nparams,2)>0
	error('ephysPipeline:argChk','Parameters must be specified as parameter/value pairs!');
end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'filedir'
			filedir=varargin{i+1};
		case 'savedir'
			savedir=varargin{i+1};
		case 'fig_title'
			fig_title=varargin{i+1};
		case 'nfft'
			nfft=varargin{i+1};
		case 'n'
			n=varargin{i+1};
		case 'min_f'
			min_f=varargin{i+1};
		case 'max_f'
			max_f=varargin{i+1};
		case 'freq_range'
			freq_range=varargin{i+1};
		case 'nphasebins'
			nphasebins=varargin{i+1};
		case 'overlap'
			overlap=varargin{i+1};
		case 'subtrials'
			subtrials=varargin{i+1};
		case 'proc_fs'
			proc_fs=varargin{i+1};
		case 'spike_channel'
			spike_channel=varargin{i+1};
		case 'spike_cluster'
			spike_cluster=varargin{i+1};
		case 'lfp_channel'
			lfp_channel=varargin{i+1};
		case 'method'
			method=varargin{i+1};
		case 'debug'
			debug=varargin{i+1};
		case 'wsigma'
			wsigma=varargin{i+1};
		case 'angles'
			angles=varargin{i+1};
		case 'padding'
			padding=varargin{i+1};

	end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(nfft)
	nfft=max([n 2^nextpow2(n)]);
else
	nfft=2^nextpow2(nfft);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DATA COLLECTION %%%%%%%%%%%%%%%%%%%%


[nsamples,ntrials,nchannels]=size(LFPDATA);

if ~isempty(padding)
	padsamples=padding*lfp_fs;
	pad=zeros(padsamples,ntrials,nchannels);
	LFPDATA=[pad;LFPDATA;pad];
end


ntimescales=length(wsigma);

% where to grab the files from?
% filter the LFP data

if length(SPIKETIMES)~=size(LFPDATA,2)
	error('ephysPipeline:spectcoherence:unequaltrials','Unequal number of trials');
end

% downsample factor

downfact=lfp_fs/proc_fs;

if mod(downfact,1)>0
	error('ephysPipeline:spectcoherence:downsamplenotinteger','Need to downsample by integer');
end

[b,a]=butter(2,[200/(lfp_fs/2)],'low');
lfp_data=downsample(filtfilt(b,a,double(LFPDATA)),downfact);

lfp_data=ephys_condition_signal(lfp_data,'l','freq_range',freq_range,'medfilt_scale',medfilt_scale,'medfilt',0,...
	'fs',proc_fs,'filt_order',filt_order);
lfp_data=squeeze(lfp_data);

[nsamples,ntrials]=size(lfp_data)
[path,name,ext]=fileparts(savedir);

savedir=fullfile(savedir,'phasehist',...
	[ 'tf (lfpch' num2str(lfp_channel) '_such' num2str(spike_channel) '_cl' num2str(spike_cluster) ' contour)'] );

if ~exist(savedir,'dir')
	mkdir(savedir);
end

dt=1/proc_fs;

[nsamples,ntrials]=size(lfp_data);
lfp_data=lfp_data';

% pre-allocate the binned spike matrix, specify at a high enough
% fs so that we don't add multiple spikes to a given bin, sampling
% rate of LFP should be more than sufficient (25e3 for Intan)

binspike_data=zeros(ntrials,nsamples,'double');
binedges=[1:nsamples]./proc_fs;

% number of rows and columns in the spectrogram

[t,f,startidx,stopidx]=getspecgram_dim(nsamples,n,overlap,nfft,proc_fs,min_f,max_f);

if ~isempty(padding)
	t=t-padding;
end

rows=length(f);
columns=length(t);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% AVERAGE SPECTRA %%%%%%%%%%%%%%%%%%%%



% for each trial grab the phase and frequency of a contour and use to generate a new phase histogram

phase_edges=linspace(0,2*pi,nphasebins);
phase_edges(1)=-inf;
fvec=startidx:stopidx;
amp_weighted_histogram=zeros(length(fvec),length(phase_edges));
cutoff=0; % to be determined

consensus_ave=zeros(rows,columns);
gabor_ave=zeros(rows,columns);

t_samples=t*proc_fs;

for i=1:ntrials

	disp([num2str(i)]);	
	consensus_tmp=zeros(rows,columns);
	gabor_tmp=zeros(rows,columns);

	currdata=lfp_data(i,:);
	spect=[];

	for j=1:ntimescales

		% get the stft and reassignment

		[stft dx]=chirp_stft(currdata,'fs',proc_fs,'n',n,...
			'overlap',overlap,'nfft',nfft,'wsigma',wsigma(j));

		% build the complex contour image

		%consensus=contour_consensus(stft,dx,angles);

		% accumulate

		%consensus_tmp=consensus_tmp+consensus./ntimescales;

		gabor_tmp=gabor_tmp+stft./ntimescales;

	end

	%consensus_ave=consensus_ave+consensus_tmp./ntrials;
	gabor_ave=gabor_ave+gabor_tmp./ntrials;

	% where are the spikes (timestamps are in seconds)?

	currtimes=unique(round(SPIKETIMES{i}*proc_fs));
	newtimes=[];
	
	for j=1:length(currtimes)
		newtimes=[ newtimes find(currtimes(j)==t_samples) ];
	end

	mag=abs(gabor_tmp);
	phase=angle(gabor_tmp);

	for j=1:2:rows
		phase(j,:)=phase(j,:)+pi;
	end

	phase=mod(phase,2*pi);

	% get likelihood and phase

	spikemags=mag(fvec,newtimes);
	spikephases=phase(fvec,newtimes);

	for j=1:length(newtimes)

		[density,phasebins]=histc(spikephases(:,j),phase_edges);

		for k=1:length(fvec)
			amp_weighted_histogram(k,phasebins(k))=...
				amp_weighted_histogram(k,phasebins(k))+spikemags(k,j);
		end
	end


end

figure();imagesc(amp_weighted_histogram(5:end,:))
