function [LFPWINS_PEAK,LFPWINS_TROUGH,LFPWINS_RAND]=ephys_su_lfp_coherence_ifr(LFPCHANNEL,SUCHANNEL,SUCLUSTER,varargin)
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

if nargin<3
	error('ephysPipeline:ifrcoherence:notenoughparams','Need 3 arguments to continue, see documentation');
end

nparams=length(varargin);
filedir=pwd;
savedir=pwd;
lfp_winextract=[.03 .03]; % msec before and after peak or trough to grab LFP
peakedges=[200 200]; % take the leading edge of the peak
troughedges=[100 100];
fig_title='noname';
randreps=100;
singletrialplots=30; % trials chosen at random to plot with fields, spikes, and IFR
debug=0;
trial_min=20;
freq_range=[30 60];
options=statset('UseParallel','Always');
peak='e';
trough='m';
medfilt_scale=1.5; % median filter scale (in ms)
lfp_fs=25e3;
hist_facecolor='none';
hist_edgecolor='m';

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
		case 'lfp_fs'
			lfp_fs=varargin{i+1};
		case 'hist_facecolor'
			hist_facecolor=varargin{i+1};
		case 'hist_edgecolor'
			hist_edgecolor=varargin{i+1};
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DATA COLLECTION %%%%%%%%%%%%%%%%%%%%

% where to grab the files from?

if isempty(savedir)
	name='noname';
else
	[path,name,ext]=fileparts(savedir);
end
savename=[ name '_lfpch_' num2str(LFPCHANNEL) '_such' num2str(SUCHANNEL) '_clust' num2str(SUCLUSTER) '_freqs' num2str(freq_range)];
load(fullfile(pwd,'../../aggregated_data.mat'),'CHANNELS','EPHYS_DATA'); % get the channel map and LFPs

% first let's get the smooth spike traces and IFR (use IFR on a trial by trial basis)

sua_mat=fullfile(pwd,'sua',['sua_channels ' num2str(SUCHANNEL) '.mat']);
load(sua_mat,'smooth_spikes','IFR','TIME','clust_spike_vec','subtrials'); % smooth spikes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SIGNAL CONDITIONING %%%%%%%%%%%%%%%%

if isempty(find(LFPCHANNEL==CHANNELS))
	error('ephysPipeline:ifrcoherence:lfpchanneldne','LFP channel %g does not exist',LFPCHANNEL);
end

lfp_data=ephys_denoise_signal(EPHYS_DATA,CHANNELS,LFPCHANNEL);

clear EPHYS_DATA;
lfp_data=ephys_condition_signal(lfp_data,'l','freq_range',freq_range,'medfilt_scale',medfilt_scale);
lfp_data=squeeze(lfp_data);

% need to account for subset of trials if used in single unit data 

lfp_data=lfp_data(:,subtrials);
[samples,trials]=size(lfp_data);
lfp_data=lfp_data';
lfp_time=[1:samples]./lfp_fs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get the IFR sampling rate

ifr_fs=1./(TIME(2)-TIME(1));

lfp_winextract=round(lfp_winextract.*lfp_fs);
winlength=length([-lfp_winextract(1):lfp_winextract(2)]);

% take IFR zscore from each trial, check for peaks and troughs
% extract LFPs triggered on "peaks" or "bursts" and "pauses"

% peak detection: check for zero crossings in the positive-going direction
% trough detection:  check for zero crossings in the negative-going direction

% for peaks and troughs, we can set the magnitude of change for the zero-crossings, i.e. x[i+1]>>x[i] | x[i+1]<<x[i]
% this will detect the **leading edge** of peaks and troughs, so interpret accordingly


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PREALLOCATION FOR WINDOW EXTRACTION %

[ntrials,samples]=size(IFR{1}{SUCLUSTER});
ifr=zeros(ntrials,samples);
ifr_time=[1:samples]./ifr_fs;

% ifr and lfp amplitudes in a trial x sample matrix

ifr_data=IFR{1}{SUCLUSTER};

% need to account for subset of trials if used in single unit data 

%lfp_data=lfp_data(subtrials,:);

spike_data=clust_spike_vec{1}{SUCLUSTER};

% vector for peak/trough detection

indx=[1:samples-1];
LFPWINS_TROUGH.waveforms=[];
LFPWINS_PEAK.waveforms=[];
LFPWINS_RAND.waveforms=[];

% count the number of extractions in each case to pre-allocate (major speedup)

zerocount=0;
poscount=0;
randcount=0;

for i=1:ntrials

	%normifr(i,:)=zscore(ifr_data(i,:));
	normifr(i,:)=ifr_data(i,:);
	normlfp(i,:)=zscore(lfp_data(i,:));
	phaselfp(i,:)=angle(hilbert(normlfp(i,:)));

	% unwrap the phase

	phaselfp(i,:)=mod(unwrap(phaselfp(i,:)),2*pi);

	currifr=normifr(i,:);
	currlfp=normlfp(i,:);

	zerocross{i}=find(currifr(indx)>troughedges(1) & currifr(indx+1)<troughedges(2))+1; 
	poscross{i}=find(currifr(indx)<peakedges(1) & currifr(indx+1)>peakedges(2))+1;

	ncross(i)=length(zerocross);

	% grab LFP windows around the troughs

	for j=1:length(zerocross{i})

		if lower(trough(1))=='m'
			trough_end=min(find(currifr(zerocross{i}(j)+1:end)>troughedges(2)))+zerocross{i}(j);

			if ~isempty(trough_end)
				zerocross{i}(j)=round((zerocross{i}(j)+trough_end)/2);
			else
				zerocross{i}(j)=[];
				continue;
			end
		end
		
		lfpcenter=round((zerocross{i}(j)/ifr_fs)*lfp_fs);

		% attempt to find the end of the trough

		startidx=lfpcenter-lfp_winextract(1);
		stopidx=lfpcenter+lfp_winextract(2);

		if startidx>0 & stopidx<length(currlfp)
			zerocount=zerocount+1;
		end

	end

	% grab LFP windows around the peaks

	for j=1:length(poscross{i})

		if lower(peak(1))=='m'

			peak_end=min(find(currifr(poscross{i}(j)+1:end)<peakedges(2)))+poscross{i}(j);

			if ~isempty(peak_end)
				poscross{i}(j)=round((poscross{i}(j)+peak_end)/2);
			else
				poscross{i}(j)=[];
				continue;
			end
		end

		lfpcenter=round((poscross{i}(j)/ifr_fs)*lfp_fs);
		startidx=lfpcenter-lfp_winextract(1);
		stopidx=lfpcenter+lfp_winextract(2);

		if startidx>0 & stopidx<length(currlfp)
			poscount=poscount+1;
		end

	end

	population=[1:samples];
	randcross{i}=population(randsample(length(population),max([length(poscross{i}) length(zerocross{i})])));

	for j=1:length(randcross{i})

		lfpcenter=round((randcross{i}(j)/ifr_fs)*lfp_fs);

		startidx=lfpcenter-lfp_winextract(1);
		stopidx=lfpcenter+lfp_winextract(2);

		if startidx>0 & stopidx<length(currlfp)
			randcount=randcount+1;
		end

	end

	% show the results if debugging

	if debug
		lfpfig=figure();
		plot(lfp_time,normifr(i,:));
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
LFPWINS_RAND.waveforms=zeros(winlength,randcount,'single');
LFPWINS_TROUGH.phaseangle=zeros(zerocount,1,'single');
LFPWINS_PEAK.phaseangle=zeros(poscount,1,'single');
LFPWINS_RAND.phaseangle=zeros(randcount,1,'single');

zerocount=1;
poscount=1;
randcount=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WINDOW EXTRACTION %%%%%%%%%%%%%%%%%%


for i=1:ntrials

	currifr=normifr(i,:);
	currlfp=normlfp(i,:);
	
	% find zero crossings where IFR > 3*Std and next point < 0
	

	% grab LFP windows around the troughs

	for j=1:length(zerocross{i})

		lfpcenter=round((zerocross{i}(j)/ifr_fs)*lfp_fs);

		startidx=lfpcenter-lfp_winextract(1);
		stopidx=lfpcenter+lfp_winextract(2);

		if startidx>0 & stopidx<length(currlfp)
			LFPWINS_TROUGH.waveforms(1:winlength,zerocount)=currlfp(startidx:stopidx)';
			LFPWINS_TROUGH.phaseangle(zerocount)=phaselfp(i,lfpcenter);
			zerocount=zerocount+1;
		end

	end

	% grab LFP windows around the peaks

	for j=1:length(poscross{i})

		lfpcenter=round((poscross{i}(j)/ifr_fs)*lfp_fs);

		startidx=lfpcenter-lfp_winextract(1);
		stopidx=lfpcenter+lfp_winextract(2);

		if startidx>0 & stopidx<length(currlfp)
			LFPWINS_PEAK.waveforms(1:winlength,poscount)=currlfp(startidx:stopidx)';
			LFPWINS_PEAK.phaseangle(poscount)=phaselfp(i,lfpcenter);
			poscount=poscount+1;
		end

	end

	for j=1:length(randcross{i})

		lfpcenter=round((randcross{i}(j)/ifr_fs)*lfp_fs);

		startidx=lfpcenter-lfp_winextract(1);
		stopidx=lfpcenter+lfp_winextract(2);

		if startidx>0 & stopidx<length(currlfp)
			LFPWINS_RAND.waveforms(1:winlength,poscount)=currlfp(startidx:stopidx)';
			LFPWINS_RAND.phaseangle(randcount)=phaselfp(i,lfpcenter);
			randcount=randcount+1;
		end

	end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RANDOM POINT WINDOW EXTRACTION %%%%%


% take random positions in each trial for significance evaluation

disp('Computing random repetitions')
mean_rand=zeros(randreps,winlength);

parfor i=1:randreps

	randcross={};
	currlfp=[];

	counter=0;
	for j=1:ntrials
		
		currlfp=normlfp(j,:);
		randcross{j}=randsample(samples,ncross(j)); % random samples
		
		for k=1:length(randcross{j})

			lfpcenter=round((randcross{j}(k)/ifr_fs)*lfp_fs);

			startidx=lfpcenter-lfp_winextract(1);
			stopidx=lfpcenter+lfp_winextract(2);

			if startidx>0 & stopidx<length(currlfp)
				counter=counter+1;
				%rand_waveforms=[rand_waveforms currlfp(startidx:stopidx)'];
			end

		end

	end

	rand_waveforms=zeros(winlength,counter);

	counter=1;
	for j=1:ntrials

		currlfp=normlfp(j,:);
		
		for k=1:length(randcross{j})

			lfpcenter=round((randcross{j}(k)/ifr_fs)*lfp_fs);

			startidx=lfpcenter-lfp_winextract(1);
			stopidx=lfpcenter+lfp_winextract(2);

			if startidx>0 & stopidx<length(currlfp)
				rand_waveforms(1:winlength,counter)=currlfp(startidx:stopidx);
				counter=counter+1;
				%rand_waveforms=[rand_waveforms currlfp(startidx:stopidx)'];
			end

		end

	end

	mean_rand(i,:)=mean(rand_waveforms,2);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[lfptrials,lfpsample]=size(lfp_data);

% plot the mean spike triggered lfp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTTING CODE %%%%%%%%%%%%%%%%%%%%%%

timevec=[-lfp_winextract(1):lfp_winextract(2)]./lfp_fs*1e3;
mean_lfp_peak=mean(LFPWINS_PEAK.waveforms,2);

stalfp_fig=figure('Position',[0 0 800 600],'Visible','off');

if isempty(savedir)
	set(stalfp_fig,'Visible','on');
else
	
	savedir=fullfile(savedir,'coherence',['stalfp (ch' num2str(SUCHANNEL) '_cl' num2str(SUCLUSTER) ')']);

	if ~exist(savedir,'dir')
		mkdir(savedir);
	end

end

% 99 % confidence interval

conf_upper=ones(size(timevec)).*prctile(mean_rand(:),97.5);
conf_lower=ones(size(timevec)).*prctile(mean_rand(:),2.5);

fill([timevec fliplr(timevec)],[conf_upper fliplr(conf_lower)],'k','FaceColor',[.8 .8 .8],'EdgeColor','none');

set(gca,'Layer','Top')
h=[];
hold on

[samples,trials]=size(LFPWINS_PEAK.waveforms);

if trials>trial_min

	h(end+1)=plot(timevec,mean_lfp_peak,'m-','linewidth',2);

	% jackknife the standard error, sqrt(n-1/n*SIGMA[mean_j - mean_sample]^2)
	
	LFPWINS_PEAK.sem=[std(LFPWINS_PEAK.waveforms')./sqrt(trials)]';
	plot(timevec,mean_lfp_peak-1.96*LFPWINS_PEAK.sem,'m--');
	plot(timevec,mean_lfp_peak+1.96*LFPWINS_PEAK.sem,'m--');
	axis tight;

end

% label with N, need a statistical test here...

mean_lfp_trough=mean(LFPWINS_TROUGH.waveforms,2);
[samples,trials]=size(LFPWINS_TROUGH.waveforms);

if trials>trial_min

	h(end+1)=plot(timevec,mean_lfp_trough,'g-','linewidth',2);

	LFPWINS_TROUGH.sem=[std(LFPWINS_TROUGH.waveforms')./sqrt(trials)]';

	plot(timevec,mean_lfp_trough-1.96*LFPWINS_TROUGH.sem,'g--');
	plot(timevec,mean_lfp_trough+1.96*LFPWINS_TROUGH.sem,'g--');
	axis tight;

end

%h(3)=plot(timevec,ones(size(timevec)).*prctile(mean_rand(:),97.5),'--','color',[.7 .7 .7],'linewidth',2);
%plot(timevec,ones(size(timevec)).*prctile(mean_rand(:),2.5),'--','color',[.7 .7 .7],'linewidth',2);
set(gca,'TickDir','out','FontName','Helvetica','FontSize',18,'linewidth',1.25,'TickLength',[.025 .025])
box off
ylabel('STA-LFP (zscore)','FontName','Helvetica','FontSize',20);
xlabel('Time (ms)','FontName','Helvetica','FontSize',20);
title(['LFP (channel ' num2str(LFPCHANNEL) ') SU (channel ' num2str(SUCHANNEL) ' cluster ' num2str(SUCLUSTER) ')']);
% include channel, bird ID

L=legend(h,'FR Peak','FR Trough');
set(L,'FontName','Helvetica','FontSize',15,'Location','NorthWest')
legend boxoff

if ~isempty(savedir)
	multi_fig_save(stalfp_fig,savedir,[ savename '_spect' ],'eps,png');
	close([stalfp_fig]);
end

% include statistical tests and polar histograms of phase values!, that should be it...
% now compute polar histograms of the phase angles at the peaks and troughs

if length(LFPWINS_PEAK.phaseangle)>0
	[rtest_p.peak]=circ_rtest(LFPWINS_PEAK.phaseangle);
else
	rtest_p.peak=NaN;
end

if length(LFPWINS_TROUGH.phaseangle)>0
	[rtest_p.trough]=circ_rtest(LFPWINS_TROUGH.phaseangle);
else
	rtest_p.trough=NaN;
end

if length(LFPWINS_RAND.phaseangle)>0
	[rtest_p.rand]=circ_rtest(LFPWINS_RAND.phaseangle);
else
	rtest_p.trough=NaN;
end


stats_fig=figure('Position',[0 0 800 600],'Visible','off');
if isempty(savedir)
	set(stats_fig,'Visible','on');
end

subplot(1,3,1);
pretty_polar(LFPWINS_PEAK.phaseangle,12,'x_label',{'FR Peak';...
	['Rayleigh test: p ' sprintf('%.3f',rtest_p.peak)]},'fignum',stats_fig,...
	'edgecolor','m','facecolor',[.9333 .51 .9333]);
subplot(1,3,2);
pretty_polar(LFPWINS_TROUGH.phaseangle,12,'x_label',{'FR Trough',...
	['Rayleigh test: p ' sprintf('%.3f',rtest_p.trough)]},'fignum',stats_fig,'labels',{'','','',''},...
	'edgecolor',[0 .545 .275],'facecolor','g');
subplot(1,3,3);
pretty_polar(LFPWINS_RAND.phaseangle,12,'x_label',{'Random',...
	['Rayleigh test: p ' sprintf('%.3f',rtest_p.rand)]},'fignum',stats_fig,'labels',{'','','',''},...
	'edgecolor',[.2 .2 .2],'facecolor',[.6 .6 .6]);

if ~isempty(savedir)
	multi_fig_save(stats_fig,savedir,...
		[ savename '_circstats' ],'eps,png');
	close([stats_fig]);
	save(fullfile(savedir,[savename '.mat']),'LFPWINS_PEAK','LFPWINS_TROUGH','LFPWINS_RAND',...
		'rtest_p');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

