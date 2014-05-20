function [abscoh,t,f]=ephys_su_lfp_coherence_gaborhist(LFPDATA,SPIKETIMES,varargin)
% THIS FUNCTION IS CURRENTLY INCOMPLETE
%computes coherograms between LFPs and single units
%
%	[abscoh,t,f]=ephys_su_lfp_coherence_gaborhist(lfpdata,spiketimes,HISTOGRAM,varargin)
%	
%	lfpdata
%	matrix with lfp data (samples x trials)
%
%	spiketimes
%	cell array with spike times (in seconds)
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

if nargin<2
	error('ephysPipeline:tfcoherence:notenoughparams','Need 2 arguments to continue, see documentation...');
end

nparams=length(varargin);

savedir=pwd;

colors='hot';

debug=0;
nfft=1024;
n=512;
overlap=511;
min_f=1;
max_f=100;
proc_fs=500; % processing frequency, will downsample lfp

freq_range=[20 120]; % let's just refilter
filt_order=7;
filt_type='bandpass';
lfp_fs=25e3; % default Intan sampling rate
trial_range=[];
medfilt_scale=1.5; % median filter scale (in ms)
nphasebins=13;
method='contour';

spike_channel='null';
spike_cluster='null';
lfp_channel='null';

subtrials=[];

angles=-pi/4:pi/16:pi/4; % for contour image
%wsigma=.05:.02:.15; %timescales in milliseconds
wsigma=.15;
padding=.5;
debug=0;
ifr_thresh=[];

if mod(nparams,2)>0
	error('ephysPipeline:argChk','Parameters must be specified as parameter/value pairs!');
end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'savedir'
			savedir=varargin{i+1};
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
		case 'ifr_thresh'
			ifr_thresh=varargin{i+1};
		case 'binary'
			binary=varargin{i+1};

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

for i=1:nchannels
	lfp_data(:,:,i)=filtfilt(b,a,double(LFPDATA(:,:,i)));
end

lfp_data=downsample(lfp_data,downfact);

lfp_data=ephys_condition_signal(lfp_data,'l','freq_range',freq_range,'medfilt_scale',medfilt_scale,'medfilt',0,...
	'fs',proc_fs,'filt_order',filt_order,'filt_type',filt_type);
lfp_data=squeeze(lfp_data);

[nsamples,ntrials]=size(lfp_data);
[path,name,ext]=fileparts(savedir);


savedir=fullfile(savedir,'phasehist',...
	[ 'tf (lfpch' num2str(lfp_channel) '_such' num2str(spike_channel) '_cl' num2str(spike_cluster) ' ' method ')'] );

savefilename=[ name '_tfcoherence_lfpch'...
       	num2str(lfp_channel) '_such' num2str(spike_channel) '_cl ' num2str(spike_cluster)];

if ~exist(savedir,'dir')
	mkdir(savedir);
end

[nsamples,ntrials,nchannels]=size(lfp_data);
%lfp_data=lfp_data';

% pre-allocate the binned spike matrix, specify at a high enough
% fs so that we don't add multiple spikes to a given bin, sampling
% rate of LFP is overkill (25e3 for Intan)

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

%phase_edges=linspace(0,2*pi,nphasebins);
phase_edges=0:2*pi/nphasebins:2*pi;
%phase_edges(1)=-inf;

% fvec defines the frequency range we're interested in

fvec=startidx:stopidx;
amp_weighted_histogram=zeros(length(fvec),length(phase_edges));
cutoff=0; % to be determined

consensus_ave=zeros(rows,columns);
gabor_ave=zeros(rows,columns);

t_samples=round(t*proc_fs);
size(lfp_data)

for i=1:ntrials

	disp([ 'Trial ' num2str(i)]);	

	consensus_tmp=zeros(rows,columns);
	gabor_tmp=zeros(rows,columns);

	currdata=squeeze(lfp_data(:,i,:));
	spect=[];

	parfor j=1:ntimescales

	       % get the stft and reassignment

	       consensus_channel_tmp=zeros(rows,columns);
	       gabor_channel_tmp=zeros(rows,columns);

	       stft=[];
	       consensus=[];

	       for k=1:nchannels

		       [stft dx]=chirp_stft(currdata(:,k),'fs',proc_fs,'n',n,...
			       'overlap',overlap,'nfft',nfft,'wsigma',wsigma(j));

		       % build the complex contour image, combine across channels

		       if lower(method(1))=='c'
			       consensus=contour_consensus(stft,dx,angles);
		       else
			       consensus=zeros(size(stft));
		       end

		       consensus_channel_tmp=consensus_channel_tmp+consensus./nchannels;
		       gabor_channel_tmp=gabor_channel_tmp+stft./nchannels;

	       end

	       % accumulate

		consensus_tmp=consensus_tmp+consensus_channel_tmp./ntimescales;
		gabor_tmp=gabor_tmp+gabor_channel_tmp./ntimescales;

	end

	if debug

		figure(1)
		subplot(2,1,1);imagesc(abs(consensus_tmp(fvec,:)));axis xy;colorbar();
		subplot(2,1,2);imagesc(angle(consensus_tmp(fvec,:)));axis xy;colorbar();

		figure(2);
		subplot(2,1,1);imagesc(abs(gabor_tmp(fvec,:)));axis xy;colorbar();
		subplot(2,1,2);imagesc(mod(unwrap(angle(gabor_tmp(fvec,:)),[],2),2*pi));axis xy;colorbar();

		pause();

	end

	consensus_ave=consensus_ave+consensus_tmp./ntrials;
	gabor_ave=gabor_ave+gabor_tmp./ntrials;

	%gabor_ave=spectrogram(currdata,n,overlap,nfft);
	% where are the spikes (timestamps are in seconds)?

	currtimes=round(SPIKETIMES{i}*proc_fs);

	if ~isempty(ifr_thresh)

		currifr=zeros(1,length(currtimes));

		tmptimes=[-inf currtimes inf]

		counter=1;
		for j=2:length(tmptimes)-1
			prevspike=abs(tmptimes(j)-tmptimes(j-1));
			nextspike=abs(tmptimes(j)-tmptimes(j+1));
			currifr(counter)=1./(min([prevspike nextspike])./proc_fs);
			counter=counter+1;
		end

		currifr(currifr==-inf|currifr==inf)=0;
		currifr
		currtimes=currtimes(currifr>ifr_thresh);

	end
	
	% get the column of the spectrogram where each spike occurs
	
	newtimes=[];

	%tmp=t_samples(t_samples>0);
	
	for j=1:length(currtimes)
		newtimes=[ newtimes find(currtimes(j)==t_samples) ];
	end

	%newtimes=setdiff(tmp,currtimes);

	% mag and phase, offset by pi

	if lower(method(1))=='c'
		mag=abs(consensus_tmp);
		phase=angle(consensus_tmp);
	else
		mag=abs(gabor_tmp);
		phase=angle(gabor_tmp);
	end

	for j=1:2:rows
		phase(j,:)=phase(j,:)+pi;
	end

	phase=mod(phase+pi,2*pi);

	% get the magnitude and phase at the spike times

	spikemags=mag(fvec,newtimes);

	if binary
		spikemags=ones(size(spikemags));
	end

	spikephases=phase(fvec,newtimes);

	% for the weighting amp for now

	for j=1:length(newtimes)

		% check phases throughout the entire frequency range

		[density,phasebins]=histc(spikephases(:,j),phase_edges);

		% add to the amp weighted histogram using the amplitude of the 
		% time-frequency image (Gabor or consensus)

		for k=1:length(fvec)
			amp_weighted_histogram(k,phasebins(k))=...
				amp_weighted_histogram(k,phasebins(k))+spikemags(k,j);
		end

	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%

fig=figure('Visible','off');

imagesc(phase_edges(1:end-1),f(fvec),amp_weighted_histogram(:,1:end-1));
axis xy
set(gca,'XTick',[ phase_edges(1) phase_edges(ceil(length(phase_edges)/2)) phase_edges(end-1) ],...
	'XTickLabel',{'0' 'p' '2p'});
set(gca,'FontName','Symbol','FontSize',20);

xlabel('Phase','FontName','Helvetica','FontSize',14);
ylabel('Frequency (Hz)','FontName','Helvetica','FontSize',14);
colormap(hot);

box off;

multi_fig_save(fig,savedir,[savefilename '_ampweightedphase' ],'eps,png');
close([fig]);

save(fullfile(savedir,[savefilename '.mat']),'amp_weighted_histogram','t','fvec','f',...
	'phase_edges','gabor_ave','consensus_ave')
