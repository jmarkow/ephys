function [features lfp_data binspike_data]=ephys_su_lfp_coherence_tf(MICDATA,LFPDATA,SPIKETIMES,TRIALS,varargin)
%THIS FUNCTION IS INCOMPLETE USE AT YOUR OWN RISK!
%computes coherograms between LFPs and acoustic features
%
%	[abscoh,t,f]=ephys_su_lfp_audio(LFPCHANNEL,SUCHANNEL,SUCLUSTER,HISTOGRAM,varargin)
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
%		trial_range
%		only compute the coherence for these trials (leave blank for all, default: [])
%
%		medfilt_scale
%		median filter timescale in ms (default: 1.5)
%
%		nphasebins
%		number of phase bins for amplitude-weighted frequency/phase 2D histogram (default: 50)
%
% see also ephys_visual_histogra.m,ephys_su_lfp_coherence_spect.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETER COLLECTION %%%%%%%%%%%%%%%%%

if nargin<4
	error('ephysPipeline:tfcoherence:notenoughparams','Need four arguments to continue, see documentation...');
end

nparams=length(varargin);

savedir=pwd;
filedir=pwd;
fig_title='noname';
colors='jet';
debug=0;

nfft=1024;
n=250;
overlap=240;
lfp_fs=25e3;

n_audio=409;
overlap_audio=384;
audio_fs=25e3;
nfft_audio=1024;

randomizations=1e3;
p_cutoff=.05;

min_f=20;
max_f=100;
w=2;
alpha=[.001 .01 .05];
smoothplot=0;

ntapers=[];
freq_range=[300]; % let's just refilter
resample_fs=1e3; % default Intan sampling rate
trial_range=[];
medfilt_scale=1.5; % median filter scale (in ms)
nphasebins=50;
null=''; % can be inpoiss (for inhomogeneous Poisson based on IFR)
	 % poiss (homogeneous based on mean FR)
	 % wn (for white noise LFP, keep same spikes)
nullrate=[]; % for inpoiss should be a vector of firing rates computed from average IFR
	     % for poiss a scalar with average firing rate
clim=[];

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
		case 'freq_range'
			freq_range=varargin{i+1};
		case 'trial_range'
			trial_range=varargin{i+1};
		case 'nphasebins'
			nphasebins=varargin{i+1};
		case 'null'
			null=varargin{i+1};
		case 'nullrate'
			nullrate=varargin{i+1};
		case 'clim'
			clim=varargin{i+1};
		case 'overlap'
			overlap=varargin{i+1};

	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(nfft)
	nfft=max([n 2^nextpow2(n)]);
else
	nfft=2^nextpow2(nfft);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DATA COLLECTION %%%%%%%%%%%%%%%%%%%%

% first get the audio features

[T features]=compute_sap_features(MICDATA,'fs',audio_fs,'n',n_audio','overlap',overlap_audio,'nfft',nfft_audio);

% all features need to be resampled to match the audio features

audio_fs=1./(T(2)-T(1));

% get distance between window centers, new sampling rate

lfp_data=ephys_condition_signal(LFPDATA,'l','freq_range',freq_range,'medfilt_scale',medfilt_scale);
lfp_data=squeeze(lfp_data);
lfp_data=lfp_data(:,TRIALS);
spike_data=SPIKETIMES(TRIALS);

[p,q]=rat(lfp_fs/audio_fs);
startidx=floor((n_audio)/2); % the pad for the spectrogram

lfp_data=resample(lfp_data(startidx:end-startidx,:),q,p);

[nsamples,ntrials]=size(lfp_data)
[nsamples,ntrials]=size(features.am)

% use null spikes or white noise fields if the user desires

lfp_data=lfp_data';

% number of rows and columns in the spectrogram

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% AVERAGE SPECTRA %%%%%%%%%%%%%%%%%%%%

% bin spikes at the same fs as the fields

featurenames=fieldnames(features);

binspike_data=zeros(ntrials,nsamples,'double');

% number of rows and columns in the spectrogram

[t,f,startidx,stopidx]=getspecgram_dim(nsamples,n,overlap,nfft,lfp_fs,min_f,max_f);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% AVERAGE SPECTRA %%%%%%%%%%%%%%%%%%%%

% bin spikes at the same fs as the fields

disp('Binning spikes...');

for i=1:ntrials

	spike_locs=round(spike_data{i}*audio_fs);
	spike_locs(spike_locs>nsamples)=[];

	if length(unique(spike_locs))<length(spike_locs)
		warning('ephysPipeline:spectcoherence:toomanyspikesperbin',...
			'Multiple spikes in a single bin, try increasing fs');
	end

	binspike_data(i,spike_locs)=1;

end

% testing with just covariance first

%for i=1:length(featurenames)
%	[coh,t,f,lfp_spect_mean,spike_spect_mean,stats]=ephys_tfcoherence(lfp_data,zscore(features.(featurenames{i}))','fs',audio_fs,...
%		'nfft',nfft,'n',n,'overlap',overlap,'ntapers',ntapers,'alpha',alpha,'w',w);
%	figure();imagesc(t,f(1:50),abs(coh(1:50,:)));colorbar()
%	title(featurenames{i})
%	axis xy
%end

% start with simple cross-covariance, set significance boundaries with phase
% scrambled randomization

sigma=.01;

kernedges=[-3*sigma:1/audio_fs:3*sigma];

if exist('normpdf')>0
	kernel=normpdf(kernedges,0,sigma);
else
	kernel=(1/(sigma*sqrt(2*pi)))*exp((-(kernedges-0).^2)./(2*sigma^2));
end

[b,a]=butter(3,[100 200]/(audio_fs/2),'bandpass');
randomizations=10e3;

for i=1:length(featurenames)

	% assign threshold through Monte Carlo, randomize phase

	counter=1;
	c=zeros(ntrials,size(lfp_data,2)*2-1);
	crand=zeros(randomizations,size(lfp_data,2)*2-1);

	for j=1:ntrials

		%testsig=conv(binspike_data(j,:),kernel,'same');
		
		testsig=lfp_data(j,:);
		
		%testsig=filtfilt(b,a,lfp_data(j,:));

		[c(j,:),lags]=xcov(features.(featurenames{i})(:,j),testsig,'coeff');

		% phase scramble
	
	end

	parfor j=1:randomizations

		% select trial at random

		currtrial=randi(ntrials,1);
		testsig=lfp_data(currtrial,:);

		fftsig=fft(testsig);
		sigamp=abs(fftsig);
		sigtheta=angle(fftsig);
		scrtheta=angle(fft(rand(size(testsig))));

		scrsig=real(ifft(sigamp.*exp(1i.*scrtheta)));

		crand(j,:)=xcov(features.(featurenames{i})(:,currtrial),scrsig,'coeff');

	end

	% to start just phase scramble the lfps or spikes given the 
	% feature vector, take max and min at each lag
	% to define significance

	sem=std(c)./sqrt(ntrials);

	figure();plot(lags,mean(c),'m-','linewidth',1.25)
	hold on
	plot(lags,mean(c)+1.96*sem,'m--');
	plot(lags,mean(c)-1.96*sem,'m--');
	plot(lags,prctile(crand,5),'k--','color',[.3 .3 .3]);
	plot(lags,prctile(crand,95),'k--','color',[.3 .3 .3]);
	title(featurenames{i})

	pause();

end
