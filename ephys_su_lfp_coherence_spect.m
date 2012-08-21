function [abscoh,m_coherence_null,m_coherence_err,freqs]=ephys_su_lfp_coherence_spect(LFPCHANNEL,SUCHANNEL,SUCLUSTER,varargin)
%computes coherency spectra between fields and single units
%
%	ephys_su_lfp_coherence_spect(LFPCHANNEL,SUCHANNEL,SUCLUSTER,HISTOGRAM,varargin)
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
%		nfft
%		spectrum nfft (default: 2^nextpow2(length(signal)))
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
%		alpha
%		confidence line that coherence is significant (default .001)
%
%		ntapers
%		number of tapers to use (default: 2*w-1)
%
%		beta
%		parameter for z-transform of coherence (default: 1.5)
%
%		freq_range
%		frequency range for LFP filtering (default: 300)
%
%		lfp_fs
%		lfp sampling rate (default: 25e3)
%
%		medfilt_scale
%		median filter scale in ms (default: 1.5)
%
% see also ephys_su_lfp_coherence_ifr.m, ephys_su_lfp_coherence_tf.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETER COLLECTION %%%%%%%%%%%%%%%%%

if nargin<3
	error('ephysPipeline:spectcoherence:notenoughparams','Need 3 arguments to continue, see documentation');
end

nparams=length(varargin);

filedir=pwd;
savedir=pwd;
nfft=[];
min_f=10;
max_f=100;
w=3;

% about 200 trials seems to be our limit here...

ntapers=[];
beta=1.5; % parameter for z-transforming coherence
freq_range=[300]; % let's just refilter
lfp_fs=25e3; % default Intan sampling rate
trial_range=[];
medfilt_scale=1.5; % median filter scale (in ms)
alpha=.001; % alpha for null hypothesis line 

if mod(nparams,2)>0
	error('ephysPipeline:argChk','Parameters must be specified as parameter/value pairs!');
end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'savedir'
			savedir=varargin{i+1};
		case 'filedir'
			filedir=varargin{i+1};
		case 'dir'
			dir=varargin{i+1};
		case 'lfp_winextract'
			lfp_winextract=varargin{i+1};
		case 'debug'
			debug=varargin{i+1};
		case 'peakedges'
			peakedges=varargin{i+1};
		case 'troughedges'
			troughedges=varargin{i+1};
		case 'nfft'
			nfft=varargin{i+1};
		case 'freq_range'
			freq_range=varargin{i+1};
		case 'trial_range'
			trial_range=varargin{i+1};
		case 'medfilt_scale'
			medfilt_scale=varargin{i+1};
		case 'alpha'
			alpha=varargin{i+1};

	end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DATA COLLECTION %%%%%%%%%%%%%%%%%%%%


sua_mat=fullfile(filedir,'sua',['sua_channels ' num2str(SUCHANNEL) '.mat']);

load(fullfile(filedir,'aggregated_data.mat'),'CHANNELS','EPHYS_DATA'); % get the channel map and LFPs
load(sua_mat,'smooth_spikes','clust_spike_vec','subtrials'); % smooth spikes

if isempty(find(LFPCHANNEL==CHANNELS))
	error('ephysPipeline:spectcoherence:lfpchanneldne','LFP channel %g does not exist',LFPCHANNEL);
end

lfp_data=ephys_denoise_signal(EPHYS_DATA,CHANNELS,LFPCHANNEL);
clear EPHYS_DATA;
lfp_data=ephys_condition_signal(lfp_data,'l','freq_range',freq_range,'medfilt_scale',medfilt_scale);
lfp_data=squeeze(lfp_data);
spike_data=clust_spike_vec{1}{SUCLUSTER}; 

lfp_data=lfp_data(:,subtrials);
spike_data=clust_spike_vec{1}{SUCLUSTER}; 

[samples,trials]=size(lfp_data);
lfp_data=lfp_data';
% bin the spikes 

binspike_data=zeros(trials,samples,'double');


if isempty(nfft)
	nfft=max([samples 2^nextpow2(samples)]);
else
	nfft=2^nextpow2(nfft);
end

if isempty(ntapers)
    ntapers=2*(w)-1;
end

resolution=w*1/(samples/lfp_fs);
disp(['Resolution:  ' num2str(resolution)  ' Hz']);
disp(['NFFT:  ' num2str(nfft)]);

% take the average spectra to compute our cross spectrum

% normalize spectrum by samples and sampling rate

freqs=lfp_fs/2*linspace(0,1,nfft/2+1);

% smooth with multi-taper

[tapers,lambda]=dpss(samples,w,ntapers);

% pre allocate cell arrays to store taper and trial estimates

cross_spect_mean=zeros(1,nfft);
lfp_spect_mean=zeros(1,nfft);
spike_spect_mean=zeros(1,nfft);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% AVERAGE SPECTRA %%%%%%%%%%%%%%%%%%%%



for i=1:trials

	spike_locs=round(spike_data{i}*lfp_fs);

	if length(unique(spike_locs))<length(spike_locs)
		warning('ephysPipeline:spectcoherence:toomanyspikesperbin',...
			'Multiple spikes in a single bin, try increasing fs');
	end

	binspike_data(i,spike_locs)=1;

end


parfor i=1:trials

	% demean the data and compute the cross spectrum with the fields

	currlfp=lfp_data(i,:)-mean(lfp_data(i,:));
	currspikes=binspike_data(i,:)-mean(binspike_data(i,:));

	% multi-taper estimate
	% index by trial and taper, then we can estimate jackknife error bars

	
	cross_spect_mean_tmp=zeros(1,nfft);
	lfp_spect_mean_tmp=zeros(1,nfft);
	spike_spect_mean_tmp=zeros(1,nfft);

	for j=1:ntapers

		spectlfp=fft(currlfp.*tapers(:,j)',nfft);
		spectspikes=fft(currspikes.*tapers(:,j)',nfft);

		cross_power=(spectlfp.*conj(spectspikes))./(samples*lfp_fs);

		lfp_power=(abs(spectlfp).^2)./(samples*lfp_fs);
		spike_power=(abs(spectspikes).^2)./(samples*lfp_fs);

		% average across tapers, store across trials and just jackknife across trials

		cross_spect_mean_tmp=cross_spect_mean_tmp+cross_power./(ntapers);
		lfp_spect_mean_tmp=lfp_spect_mean_tmp+lfp_power./(ntapers);
		spike_spect_mean_tmp=spike_spect_mean_tmp+spike_power./(ntapers);
	end

	cross_spect_mean=cross_spect_mean+cross_spect_mean_tmp./(trials);
	lfp_spect_mean=lfp_spect_mean+lfp_spect_mean_tmp./(trials);
	spike_spect_mean=spike_spect_mean+spike_spect_mean_tmp./(trials);


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% STATS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



dof=ntapers*trials;

coh=cross_spect_mean./sqrt(lfp_spect_mean.*spike_spect_mean);
abscoh=abs(coh);

% asymptotic errorbars

abscoh=abscoh(1:nfft/2+1);

% take the range that we're interested in

startidx=max([find(freqs<=min_f)]);

if isempty(startidx)
	startidx=1;
end

stopidx=min([find(freqs>=max_f)]);

if isempty(stopidx)
	stopidx=length(freqs);
end

abscoh=abscoh(startidx:stopidx);
freqs=freqs(startidx:stopidx);

m_coherence_err=1.96./(sqrt(dof));
m_coherence_null=sqrt(1-alpha^(1/(dof/2-1)));

upperconf=abscoh+m_coherence_err;
lowerconf=abscoh-m_coherence_err;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTTING CODE %%%%%%%%%%%%%%%%%%%%%%

% save plots, data and parameters

[path,name,ext]=fileparts(savedir);

savedir=fullfile(savedir,'coherence',[ 'spect (ch' num2str(SUCHANNEL) '_cl' num2str(SUCLUSTER) ')']);

if ~exist(savedir,'dir')
	mkdir(savedir);
end

savefilename=[ name '_coherence_lfpch'...
       	num2str(LFPCHANNEL) '_such' num2str(SUCHANNEL) '_cl ' num2str(SUCLUSTER)];

% construct patch vectors

xdata=[freqs fliplr(freqs)];
ydata=[lowerconf fliplr(upperconf)];

spect_fig=figure('Visible','off','Position',[0 0 800 600]);
hold on;
patch(xdata,ydata,1,'facecolor',[.52 .81 .93],'edgecolor','none');
h(1)=plot(freqs,abscoh,'-','color',[0 0 1],'linewidth',1.25);
h(2)=plot(freqs,ones(size(abscoh)).*m_coherence_null,'k--','linewidth',1.25);
set(gca,'TickDir','out','TickLength',[.02 .02],'layer','top');
axis tight
box off

l=legend(h,'Coh.',['p=' num2str(alpha)]);
legend boxoff;

xlabel('Fs (Hz)');
ylabel('Coherence');
title({['LFP Ch. ' num2str(LFPCHANNEL) ' SU Ch. ' num2str(SUCHANNEL) 'cl' num2str(SUCLUSTER) ...
	' NTapers ' num2str(ntapers)' ' trials ' num2str(trials) ' res. ' num2str(resolution) ' Hz ']})
prettify_axislabels(spect_fig);

multi_fig_save(spect_fig,savedir,savefilename,'eps,png');

save(fullfile(savedir,[ savefilename '.mat']),'m_coherence_err','m_coherence_null','alpha','coh','abscoh','freqs','nfft','samples','ntapers','w');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

