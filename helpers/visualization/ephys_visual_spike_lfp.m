function FIGNUM=ephys_visual_spike_lfp(LFP,SPIKETIMES,IFR,varargin)
%visualize field and aligned spikes with IFR
%
%
%
%
%
%
%

nparams=length(varargin);

if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs');
end

fignum=[];
fig_title=[];
spike_fs=25e3;
ifr_fs=25e3;
lfp_fs=12.5e3;

for i=1:2:nparams
	switch lower(varargin{i})
		case 'fignum'
			fignum=varargin{i+1};
		case 'fig_title'
			fig_title=varargin{i+1};
		case 'lfp_fs'
			lfp_fs=varargin{i+1};
		case 'ifr_fs'
			ifs_fs=varargin{i+1};
		case 'spike_fs'
			spike_fs=varargin{i+1};
	end
end

if isempty(fignum)
	fignum=figure();
	FIGNUM=fignum;
else
	FIGNUM=fignum;
end


% define the time vector over the length of the LFP/IFR

timevec_lfp=[1:length(LFP)]./lfp_fs;
timevec_ifr=[1:length(IFR)]./ifr_fs;
% create the spike vector

spikes=SPIKETIMES(:)'./spike_fs;

spikevec=[spikes;spikes];
trialvec=[ones(size(spikes));zeros(size(spikes))];

ymin=floor(min(LFP));
ymax=ceil(max(LFP));

ax(1)=subaxis(6,1,1,1,1,3,'margin',.15,'spacingvert',0);
plot(timevec_lfp,LFP,'b-','linewidth',1.15);
set(gca,'XTick',[],'YTick',[ymin ymax],'TickDir','out','TickLength',[.02 .02],'xcolor',get(fignum,'color'));
ylim([ymin-eps ymax+eps]);
ylabel('LFP Amp. ($\mu$V)','fontsize',13,'fontname','helvetica','interpreter','latex')
box off

ax(2)=subaxis(6,1,1,4,1,1,'margin',.15,'spacingvert',0);
plot(spikevec,trialvec,'k-');
axis off

ymin=floor(min(IFR));
ymax=ceil(max(IFR));

if ymin>=ymax
	ymin=0;
	ymax=100;
end

ax(3)=subaxis(6,1,1,5,1,1,'margin',.15,'spacingvert',.05);
stairs(timevec_ifr,IFR,'r-');
set(gca,'YTick',[ymin ymax],'TickDir','out','TickLength',[.02 .02],'FontSize',13,'FontName','Helvetica','layer','top');
ylim([ymin-eps ymax+eps]);
ylabel('IFR (impulses/sec)','fontsize',13,'fontname','helvetica');
box off
xlabel('Time (s)','FontSize',13,'FontName','Helvetica')

linkaxes(ax,'x');
xlim([timevec_ifr(1) timevec_ifr(end)]);
