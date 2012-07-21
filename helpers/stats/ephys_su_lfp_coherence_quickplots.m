function ephys_su_lfp_coherence_quickplots(LFPCHANNEL,SUCHANNEL,SUCLUSTER,varargin)
%generates plots with aligned fields, spikes and IFR
%
%
%
%
%
%

nparams=length(varargin);

filedir=pwd;
savedir=pwd;
lfp_fs=25e3;
fig_title=[];
medfilt_scale=1.5; % median filter scale (in ms)
freq_range=[25 100];
singletrialplots=20;

if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs');
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
		case 'median_scale'
			median_scale=varargin{i+1};
		case 'singletrialplots'
			singletrialplots=varargin{i+1};
		case 'medfilt_scale'
			medfilt_scale=varargin{i+1};
	end
end

[path,name,ext]=fileparts(filedir);
savename=[ name '_lfpch_' num2str(LFPCHANNEL) '_such' num2str(SUCHANNEL) '_clust' num2str(SUCLUSTER) '_freqs' num2str(freq_range)];



load(fullfile(filedir,'aggregated_data.mat'),'CHANNELS','EPHYS_DATA'); % get the channel map and LFPs

% first let's get the smooth spike traces and IFR (use IFR on a trial by trial basis)

sua_mat=fullfile(filedir,'sua',['sua_channels ' num2str(SUCHANNEL) '.mat']);
load(sua_mat,'smooth_spikes','IFR','TIME','clust_spike_vec','subtrials'); % smooth spikes

% filter the lfp_data

lfp_data=ephys_denoise_signal(EPHYS_DATA,CHANNELS,LFPCHANNEL);
clear EPHYS_DATA;
lfp_data=ephys_condition_signal(lfp_data,'l','freq_range',freq_range,'medfilt_scale',medfilt_scale);
lfp_data=squeeze(lfp_data);

ifr_data=IFR{1}{SUCLUSTER};
spike_data=clust_spike_vec{1}{SUCLUSTER};
ifr_fs=1./(TIME(2)-TIME(1));

lfp_data=lfp_data(:,subtrials);
[samples,trials]=size(lfp_data);
lfp_data=lfp_data';

if isempty(fig_title)
	fig_title=[' LFPCH ' num2str(LFPCHANNEL) ' SUCH' num2str(SUCHANNEL) ...
		' SUCL' num2str(SUCLUSTER) ' FILT ' num2str(freq_range) ' NTRIALS' num2str(trials)];
end

trialvec=1:trials;
trialperm=trialvec(randperm(trials));

if ~isempty(savedir)

	savedir=fullfile(savedir,'coherence','quickplots');
	if ~exist(savedir,'dir')
		mkdir(savedir);
	end

	if singletrialplots>trials
		singletrialplots=trials;
	end

	for i=1:singletrialplots
		singletrialfig=figure('Position',[0 0 800 600],'Visible','off');
		currtrial=trialperm(i);
		currlfp=lfp_data(currtrial,:);
		currifr=ifr_data(currtrial,:);
		currspike=spike_data{currtrial}*ifr_fs;	
		ephys_visual_spike_lfp(currlfp,currspike,currifr,'fignum',singletrialfig,'lfp_fs',lfp_fs,'ifr_fs',ifr_fs);
		set(singletrialfig,'paperpositionmode','auto');
		multi_fig_save(singletrialfig,savedir,[ savename '_singletrial' num2str(currtrial) ],'eps,png');
		close([singletrialfig]);
	end
end

lfp_mean=mean(lfp_data);
lfp_sem=std(lfp_data)./sqrt(trials);

ymin=floor(min(lfp_mean-lfp_sem));
ymax=ceil(max(lfp_mean+lfp_sem));

mean_fig=figure('Position',[0 0 800 600],'Visible','off');
if isempty(savedir)
	set(mean_fig,'Visible','on');
end

lowerconf=lfp_mean-lfp_sem;
upperconf=lfp_mean+lfp_sem;
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

lowerconf=ifr_mean-ifr_sem;
upperconf=ifr_mean+ifr_sem;

xdata=[timevec_ifr fliplr(timevec_ifr)];
ydata=[lowerconf fliplr(upperconf)];

ax(2)=axes('Yaxislocation','right','xtick',[],'color','none');
%patch(xdata,ydata,1,'facecolor',[1 .71 .76],'edgecolor','none');
hold on;
plot(timevec_ifr,ifr_mean,'m-','linewidth',1.25);
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
	save(fullfile(savedir,[savename '.mat']),'ifr_mean','lfp_mean','ifr_sem','lfp_sem');
	close([mean_fig]);
end

