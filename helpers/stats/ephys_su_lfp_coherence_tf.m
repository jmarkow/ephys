function [abscoh,t,f]=ephys_su_lfp_coherence_tf(LFPDATA,SPIKETIMES,HISTOGRAM,varargin)
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
nfft=1024;
n=512;
overlap=506;
min_f=1;
max_f=100;
w=2;
alpha=[.001 .01 .05];
smoothplot=0;
proc_fs=1.25e3; % processing frequency, will downsample lfp

ntapers=[];
beta=23/20; % parameter for z-transforming coherence
freq_range=[20 125]; % let's just refilter
filt_order=6;
lfp_fs=25e3; % default Intan sampling rate
trial_range=[];
medfilt_scale=1.5; % median filter scale (in ms)
nphasebins=50;
null=''; % can be inpoiss (for inhomogeneous Poisson based on IFR)
	 % poiss (homogeneous based on mean FR)
	 % wn (for white noise LFP, keep same spikes)
nullrate=[]; % for inpoiss should be a vector of firing rates computed from average IFR
	     % for poiss a scalar with average firing rate
clim=[];

spike_channel='null';
spike_cluster='null';
lfp_channel='null';

subtrials=[];

method='slepian';

angles=-pi/4:pi/16:pi/4; % for contour image
ts=50:20:150; %timescales in milliseconds
contour_thresh=0;

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

	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(nfft)
	nfft=max([n 2^nextpow2(n)]);
else
	nfft=2^nextpow2(nfft);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DATA COLLECTION %%%%%%%%%%%%%%%%%%%%

resolution=w*1/(n/proc_fs);
disp(['Resolution:  ' num2str(resolution)  ' Hz']);
disp(['NFFT:  ' num2str(nfft)]);

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

[nsamples,ntrials]=size(lfp_data);
[path,name,ext]=fileparts(savedir);

switch lower(method(1))
	case 's'
		savedir=fullfile(savedir,'coherence',...
			[ 'tf (lfpch' num2str(lfp_channel) '_such' num2str(spike_channel) '_cl' num2str(spike_cluster) ' slepian)'] );
	case 'c'
		savedir=fullfile(savedir,'coherence',...
			[ 'tf (lfpch' num2str(lfp_channel) '_such' num2str(spike_channel) '_cl' num2str(spike_cluster) ' contour)'] );
end


if ~exist(savedir,'dir')
	mkdir(savedir);
end

dt=1/proc_fs;

% use null spikes or white noise fields if the user desires

switch lower(null)
	
	case 'inpoiss'

		% rate-varying Poisson process

		disp('Inhomogeneous poisson randomization');

		if length(nullrate)~=nsamples
			error('For null option inpoiss nullrate must be the same length as nsamples in lfp_data');
		end

		for i=1:ntrials
			SPIKETIMES{i}=find(poissrnd(nullrate.*dt))./proc_fs;
		end

		name=[name '_inpoiss'];

	case 'poiss'

		% homogeneous Poisson process

		disp('Homogeneous poisson randomization');
	
		if length(nullrate)>1
			error('For null option poiss nullrate must be a scalar');
		end

		for i=1:ntrials
			SPIKETIMES{i}=find(poissrnd(nullrate(ones(1,nsamples)).*dt))./proc_fs;
		end

		name=[name '_poiss_fr' num2str(nullrate)];

	case 'wn'

		disp('LFP white noise randomization');

		peak_value=max(abs(lfp_data(:)));
		lfp_data=randn(nsamples,ntrials).*peak_value;

		name=[name '_lfpwn'];

	case 'poiss+wn'

		% homogeneous Poisson process

		disp('Poiss+wn randomization');
	
		if length(nullrate)>1
			error('For null option poiss nullrate must be a scalar');
		end

		for i=1:ntrials
			SPIKETIMES{i}=find(poissrnd(nullrate(ones(1,nsamples)).*dt))./proc_fs;
		end

		peak_value=max(abs(lfp_data(:)));
		lfp_data=randn(nsamples,ntrials).*peak_value;

		name=[name '_lfpwn+poiss_fr' num2str(nullrate)];

	otherwise

end

savefilename=[ name '_tfcoherence_lfpch'...
       	num2str(lfp_channel) '_such' num2str(spike_channel) '_cl ' num2str(spike_cluster)];

% need to account for subset of trials if used in single unit data 

if ~isempty(trial_range)
	disp(['Truncating trials to ' num2str([trial_range(1) trial_range(end)])]);
	lfp_data=lfp_data(:,trial_range);
	SPIKETIMES=SPIKETIMES(trial_range);
end

[nsamples,ntrials]=size(lfp_data);
lfp_data=lfp_data';

% pre-allocate the binned spike matrix, specify at a high enough
% fs so that we don't add multiple spikes to a given bin, sampling
% rate of LFP should be more than sufficient (25e3 for Intan)

binspike_data=zeros(ntrials,nsamples,'double');
binedges=[1:nsamples]./proc_fs;

% number of rows and columns in the spectrogram

[t,f,startidx,stopidx]=getspecgram_dim(nsamples,n,overlap,nfft,proc_fs,min_f,max_f);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% AVERAGE SPECTRA %%%%%%%%%%%%%%%%%%%%

% bin spikes at the same fs as the fields

disp('Binning spikes...');

for i=1:ntrials
	binspike_data(i,:)=histc(SPIKETIMES{i},binedges);
end

if any(binspike_data>1)
	warning('ephysPipeline:spectcoherence:toomanyspikesperbin',...
		'Multiple spikes in a single bin, try increasing fs');
end


disp('Computing spectra and coherence');

switch lower(method(1))

	case 's' % slepian


		[coh,t,f,lfp_spect_mean,spike_spect_mean,stats,forcingfunction]=ephys_tfcoherence(lfp_data,binspike_data,'fs',proc_fs,...
			'nfft',nfft,'n',n,'overlap',overlap,'ntapers',ntapers,'alpha',alpha,'w',w);

		fid=fopen(fullfile(savedir,'alpha_log.log'),'w');

		for i=1:length(stats.null)	
			fprintf(fid,'For alpha %g coherence cutoff:\t%g\t%g (Bonferroni corr.)\n',...
				stats.alpha(i),stats.null(i),stats.null_boncorrected(i));
		end

		fprintf(fid,'Variance:\t%g',stats.var);
		fclose(fid);

		fig_title=['LFPCH' num2str(lfp_channel) ' SUCH' num2str(spike_channel)...
       			' SUCL' num2str(spike_cluster) ' NTAP' num2str(stats.ntapers) ' RES ' num2str(resolution) ' Hz' ' NTRIALS' num2str(ntrials)];

	case 'c' % contour

		[coh,t,f,lfp_spect_mean,spike_spect_mean,stats]=ephys_tfcoherence_contour(lfp_data,binspike_data,'fs',proc_fs,...
			'nfft',nfft,'n',n,'overlap',overlap,'angles',angles,'contour_thresh',contour_thresh,'ts',ts);

		fig_title=['LFPCH' num2str(lfp_channel) ' SUCH' num2str(spike_channel) ' SUCL' num2str(spike_cluster) ...
			' CONTOURTHRESH ' num2str(contour_thresh)];

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTTING CODE %%%%%%%%%%%%%%%%%%%%%%

% coherence is the absolute value of the cross spectrum over the power of the
% the two spectra

abscoh=abs(coh);

% phase, could plot over the abs value as quivers per Pesaran et al.

phasecoh=mod(unwrap(angle(coh)),2*pi);

% z-transformed coherence

% we can indicate significant points using the asymptotic values per Jarvis and Mitra (2001)
% spike spectrogram

spike_fig=figure('Visible','off','position',[0 0 round(600*nsamples/proc_fs) 800]);

spike_spect.image=spike_spect_mean./max(spike_spect_mean(:));
spike_spect.t=t;
spike_spect.f=f;

time_frequency_raster(HISTOGRAM,spike_spect,'tf_min_f',min_f,'tf_max_f',max_f,'scale','log',...
	'fig_num',spike_fig,'scalelabel','dB','fig_title',fig_title,'tfimage_colors',colors,'smoothplot',smoothplot);

set(spike_fig,'PaperPositionMode','auto');
multi_fig_save(spike_fig,savedir,[savefilename '_spikespectrogram' ],'eps,png');

close([spike_fig]);

% lfp spectrogram

lfp_fig=figure('Visible','off','position',[0 0 round(600*nsamples/proc_fs) 800]);

lfp_spect.image=lfp_spect_mean./max(lfp_spect_mean(:));
lfp_spect.t=t;
lfp_spect.f=f;

time_frequency_raster(HISTOGRAM,lfp_spect,'tf_min_f',min_f,'tf_max_f',max_f,'scale','log',...
	'fig_num',lfp_fig,'scalelabel','dB','fig_title',fig_title,'tfimage_colors',colors,'smoothplot',smoothplot);

set(lfp_fig,'PaperPositionMode','auto');
multi_fig_save(lfp_fig,savedir,[savefilename '_lfpspectrogram' ],'eps,png');

close([lfp_fig]);

% coherogram

spect_fig=figure('Visible','off','position',[0 0 round(600*nsamples/proc_fs) 800]);

plotabscoh.image=abscoh;
plotabscoh.t=t;
plotabscoh.f=f;

time_frequency_raster(HISTOGRAM,plotabscoh,'tf_min_f',min_f,'tf_max_f',max_f,'scale','linear',...
	'fig_num',spect_fig,'scalelabel','Coherence','fig_title',fig_title,'tfimage_colors',colors,'tf_clim',clim,'smoothplot',smoothplot);

set(spect_fig,'PaperPositionMode','auto');
multi_fig_save(spect_fig,savedir,savefilename,'eps,png','renderer','painters');
close([spect_fig])

% also show phase angle

phase_fig=figure('Visible','off','position',[0 0 round(600*nsamples/proc_fs) 800]);

plotabscoh.image=phasecoh;

% plot phase with HSV

time_frequency_raster(HISTOGRAM,plotabscoh,'tf_min_f',min_f,'tf_max_f',max_f,'scale','linear',...
	'fig_num',phase_fig,'scalelabel','Phase','fig_title',fig_title,'tfimage_colors',phasecolors,'smoothplot',smoothplot);

set(phase_fig,'PaperPositionMode','auto');
multi_fig_save(phase_fig,savedir,[savefilename '_phase' ],'eps,png','renderer','painters');
close([phase_fig]);

% bin the phase angle at each frequency bin, maybe ten bins to start

phase_edges=linspace(0,2*pi,nphasebins);
startidx=max([find(f<=min_f)]);

if isempty(startidx)
	startidx=1;
end

%
%
stopidx=min([find(f>=max_f)]);

% now bin at the frequencies we're interested in, after removing all spurious (coherence <p=.05)

fvec=startidx:stopidx;
amp_weighted_histogram=zeros(length(fvec),length(phase_edges));
%cutoff=stats.null(stats.alpha==.05)+stats.var;
cutoff=0;

for i=1:length(fvec)

	[density,phasebins]=histc(phasecoh(fvec(i),:),phase_edges);

	% take the amplitude for each bin

	for j=1:length(phase_edges)
		
		% where are points that are in this bin?

		currpoints=find(phasebins==j);
		tmp=abscoh(fvec(i),currpoints);
		tmp(tmp<cutoff)=0;

		amp_weighted_histogram(i,j)=sum(tmp);

	end

end

amp_fig=figure('Visible','off','position',[0 0 700 500]);
imagesc(phase_edges(1:end-1),f(fvec),amp_weighted_histogram(:,1:end-1));
axis xy
colormap(hot)
ylabel('Fs (Hz)');
xlabel('Phase');
box off
prettify_axis(amp_fig,'ticklength',[.025 .025],'linewidth',3);
prettify_axislabels(amp_fig,'fontsize',18,'font','helvetica');
axis tight;
set(gca,'XTick',[ phase_edges(1) phase_edges(ceil(length(phase_edges)/2)) phase_edges(end-1) ],...
	'XTickLabel',{'0' 'p' '2p'});
set(gca,'FontName','Symbol','FontSize',20);

pos=get(gca,'pos');
set(gca,'pos',[pos(1) pos(2) pos(3)*.95 pos(4)]);
pos=get(gca,'pos');
hc2=colorbar('location','eastoutside','position',[pos(1)+pos(3)+.01 pos(2)+pos(4)/2 .02 pos(4)/2]);
set(hc2,'linewidth',2,'FontSize',12,'FontName','Helvetica');
ylabel(hc2,'sum(abs(coherence))','FontSize',15,'FontName','Helvetica');

set(amp_fig,'PaperPositionMode','auto');
multi_fig_save(amp_fig,savedir,[savefilename '_cohweightedphase' ],'eps,png');
close([amp_fig]);

% plot the tapers as well

if lower(method(1))=='s'
	
	proc_fs
	tapers_fig=figure('Visible','off','position',[0 0 800 600]);
	plot([1:n]./proc_fs,stats.tapers,'linewidth',1.25);
	title(['Tapers, n ' num2str(stats.ntapers) ' w ' num2str(w)])
	xlabel('Time (s)');
	prettify_axislabels(tapers_fig,'fontsize',15,'font','Helvetica');
	prettify_axis(tapers_fig,'ticklength',[.02 .02],'fontsize',14,'font','helvetica','linewidth',3);
	box off
	axis tight

	multi_fig_save(tapers_fig,savedir,[savefilename '_tapers' ],'eps,png');
	close([tapers_fig]);
end

save(fullfile(savedir,[savefilename '.mat']),...
	'coh','t','f','forcingfunction','fvec','phase_edges','w','n','nfft','overlap',...
	'lfp_data','binspike_data','amp_weighted_histogram','stats');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

