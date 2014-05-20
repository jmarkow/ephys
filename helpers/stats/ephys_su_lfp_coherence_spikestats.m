function ephys_su_lfp_coherence_spikestats(LFPDATA,SPIKETIMES,varargin)
%generates plots with aligned fields, spikes and IFR
%
%
%
%
%
%

if nargin<2
	error('ephysPipeline:spikestats:notenoughparams','Need 3 arguments to continue, see documentation');
end

nparams=length(varargin);

filedir=pwd;
savedir=pwd;
lfp_fs=25e3;
fig_title=[];
medfilt_scale=1.5; % median filter scale (in ms)
freq_range=[10 40];
filt_order=5;
singletrials=1:10; % default to 1-20
time_range=[]; % if defined only visualize data in this subregion

spike_channel='null';
spike_cluster='null';
lfp_channel='null';
proc_fs=1.25e3; % processing frequency, will downsample lfp

if mod(nparams,2)>0
	error('ephysPipeline:argChk','Parameters must be specified as parameter/value pairs!');
end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'filedir'
			filedir=varargin{i+1};
		case 'savedir'
			savedir=varargin{i+1};
		case 'lfp_fs'
			lfp_fs=varargin{i+1};
		case 'freq_range'
			freq_range=varargin{i+1};
		case 'filt_order'
			filt_order=varargin{i+1};
		case 'median_scale'
			median_scale=varargin{i+1};
		case 'singletrials'
			singletrials=varargin{i+1};
		case 'medfilt_scale'
			medfilt_scale=varargin{i+1};
		case 'time_range'
			time_range=varargin{i+1};
		case 'spike_channel'
			spike_channel=varargin{i+1};
		case 'spike_cluster'
			spike_cluster=varargin{i+1};
		case 'lfp_channel'
			lfp_channel=varargin{i+1};
	end
end

downfact=lfp_fs/proc_fs;

if mod(downfact,1)>0
	error('ephysPipeline:spectcoherence:downsamplenotinteger','Need to downsample by integer');
end

if proc_fs<500
	warning('ephysPipeline:aliasing',...
		'Anti-aliasing filter set to 200 Hz, consider increasing the sampling rate');
end

[path,name,ext]=fileparts(filedir);
savename=[ name '_lfpch_' num2str(lfp_channel) '_such' num2str(spike_channel) '_clust' num2str(spike_cluster) '_freqs' num2str(freq_range)];

downfact=lfp_fs/proc_fs;

if mod(downfact,1)>0
	error('ephysPipeline:spectcoherence:downsamplenotinteger','Need to downsample by integer');
end

[b,a]=butter(2,[200/(lfp_fs/2)],'low');
lfp_data=downsample(filtfilt(b,a,double(LFPDATA)),downfact);

lfp_data=ephys_condition_signal(lfp_data,'l','freq_range',freq_range,'medfilt_scale',medfilt_scale,'medfilt',1,...
	'fs',proc_fs,'filt_order',filt_order);
lfp_data=squeeze(lfp_data);

lfp_fs=proc_fs;
[samples,trials]=size(lfp_data);
lfp_data=lfp_data';

ifr_fs=lfp_fs;

% get the ifr

ifr_data=zeros(size(lfp_data));

for i=1:length(SPIKETIMES)
	ifr_data(i,:)=ephys_ifr(SPIKETIMES{i}*lfp_fs,samples,lfp_fs);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SIGNAL CONDITIONING %%%%%%%%%%%%%%%%


if isempty(fig_title)
	fig_title=[' LFPCH ' num2str(lfp_channel) ' SUCH' num2str(spike_channel) ...
		' SUCL' num2str(spike_cluster) ' FILT ' num2str(freq_range) ' NTRIALS' num2str(trials)];
end


colorline=figure('visible','off');

x=[1:samples]./ifr_fs;
x=[x;x];

y=zscore(mean(lfp_data));
y=[y;y];

z=zeros(size(x));

col=mean(ifr_data);
col=[col;col];

surface(x,y,z,col,'facecol','no','edgecol','interp','linew',2);

grays=colormap('gray');
grays(1:10,:)=[];

colormap(1-grays);

ylabel('Mean LFP (STD)');
xlabel('Time (s)');
box off

axis tight;
prettify_axis(gca,'LineWidth',2,'TickLength',[.025 .025],'FontName','Helvetica','FontSize',20);
prettify_axislabels(gca,'FontName','Helvetica','FontSize',20);
savedir=fullfile(savedir,'coherence',[ 'spikestats (ch' num2str(spike_channel) '_cl' num2str(spike_cluster) ')' ]);

if ~exist(savedir,'dir')
	mkdir(savedir);
end

multi_fig_save(colorline,savedir,[ savename '_colorline' ],'eps,png');

close([colorline]);

% single trials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SINGLE TRIALS %%%%%%%%%%%%%%%%%%%%%%

trialvec=1:trials;
trialperm=trialvec(randperm(trials));

% add an options to isolate a given time region and color code according to phase if the user wants


singletrialsavedir=fullfile(savedir,'singletrials');

if ~exist(singletrialsavedir,'dir')
	mkdir(singletrialsavedir);
end

% remove any trials if we exceed the number of actual trials in the data

singletrials(singletrials>trials)=[];

for i=singletrials
	
	% collect what we need for plotting

	currlfp=lfp_data(i,:);
	currifr=ifr_data(i,:);
	currspike=SPIKETIMES{i}*ifr_fs;	
	currphase=mod(unwrap(angle(hilbert(lfp_data))),2*pi);

	singletrialfig=figure('Position',[0 0 800 600],'Visible','off');
	ephys_visual_spike_lfp(currlfp,currspike,currifr,'fignum',singletrialfig,...
		'lfp_fs',lfp_fs,'ifr_fs',ifr_fs,'spikevis','ifr','time_range',time_range);
	set(singletrialfig,'paperpositionmode','auto');
	multi_fig_save(singletrialfig,singletrialsavedir,[ savename '_IFRsingletrial' num2str(i) ],'eps,png');
	close([singletrialfig]);

	singletrialfig=figure('Position',[0 0 800 600],'Visible','off');
	ephys_visual_spike_lfp(currlfp,currspike,currifr,'fignum',singletrialfig,...
		'lfp_fs',lfp_fs,'ifr_fs',ifr_fs,'spikevis','ticks','time_range',time_range);
	set(singletrialfig,'paperpositionmode','auto');
	multi_fig_save(singletrialfig,singletrialsavedir,[ savename '_SPIKESsingletrial' num2str(i) ],'eps,png');
	close([singletrialfig]);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MEAN PLOTTING %%%%%%%%%%%%%%%%%%%%%%

lfp_mean=mean(lfp_data);
lfp_sem=std(lfp_data)./sqrt(trials);

ymin=floor(min(lfp_mean-lfp_sem));
ymax=ceil(max(lfp_mean+lfp_sem));

mean_fig=figure('Position',[0 0 800 600],'Visible','off');
if isempty(savedir)
	set(mean_fig,'Visible','on');
end

lowerconf=lfp_mean-1.96*lfp_sem;
upperconf=lfp_mean+1.96*lfp_sem;
timevec_lfp=[1:length(lfp_mean)]./lfp_fs;

xdata=[timevec_lfp fliplr(timevec_lfp)];
ydata=[lowerconf fliplr(upperconf)];

ax(1)=gca;
%patch(xdata,ydata,1,'facecolor',[1 .71 .76],'edgecolor','k');
plot(timevec_lfp,lfp_mean,'g--','linewidth',1.05);
hold on;
set(gca,'YTick',[ymin ymax],'TickDir','out','TickLength',[.02 .02],'layer','top','ycolor','g',...
	'linewidth',2,'fontsize',12,'fontname','helvetica');
ylim([ymin-eps ymax+eps]);
title([fig_title],'fontname','helvetica','fontsize',15)
box off
lab_axis=ylabel('LFP Amp. ($\mu$V)','fontsize',15,'fontname','helvetica','interpreter','latex');

xlabel('Time (s)','FontSize',13,'FontName','Helvetica');

ifr_mean=mean(ifr_data);
ifr_sem=std(ifr_data)./sqrt(trials);

timevec_ifr=[1:length(ifr_mean)]./ifr_fs;

ymin=floor(min(ifr_mean-ifr_sem));
ymax=ceil(max(ifr_mean+ifr_sem));

if ymin>=ymax
	ymin=0;
	ymax=100;
end

lowerconf=ifr_mean-1.96*ifr_sem;
upperconf=ifr_mean+1.96*ifr_sem;

xdata=[timevec_ifr fliplr(timevec_ifr)];
ydata=[lowerconf fliplr(upperconf)];

ax(2)=axes('Yaxislocation','right','xtick',[],'color','none');
%patch(xdata,ydata,1,'facecolor',[1 .71 .76],'edgecolor','none');
hold on;
plot(timevec_ifr,ifr_mean,'m-','linewidth',1);
set(ax(2),'YTick',[ymin ymax],'TickDir','out','TickLength',[.02 .02],'layer','top','ycolor','m',...
	'linewidth',2,'fontsize',12,'fontname','helvetica');
ylim([ymin-eps ymax+eps]);
ylabel('IFR (impulses/sec)','fontsize',15,'fontname','helvetica')
box off

%
linkaxes(ax,'x');
xlim([timevec_ifr(1) timevec_ifr(end)]);

if ~isempty(savedir)
	set(mean_fig,'paperpositionmode','auto');
	multi_fig_save(mean_fig,savedir,[ savename '_mean' ],'eps,png');
	save(fullfile(savedir,[savename '.mat']),'ifr_mean','lfp_mean','ifr_sem','lfp_sem',...
		'lfp_data','ifr_data','SPIKETIMES');
	close([mean_fig]);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

