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
lfp_fs=25e3;
spikevis='ticks'; % can be ticks, IFR or both
time_range=[];
phase=[]; % needs to be defined for each sample of the LFP, will also color spikes
          % and IFR
phase_colors='hsv'; % can defined a colormap defined as value x [rgb] matrix

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
		case 'phase'
			phase=varargin{i+1};
		case 'phase_colors'
			phase_colors=varargin{i+1};
		case 'spike_fs'
			spike_fs=varargin{i+1};
		case 'spikevis'
			spikevis=varargin{i+1};
		case 'time_range'
			time_range=varargin{i+1};
	end
end

spikevis=lower(spikevis);

% TODO phase coding

if isempty(fignum)
	fignum=figure();
	FIGNUM=fignum;
else
	FIGNUM=fignum;
end

% define the time vector over the length of the LFP/IFR

timevec_lfp=[1:length(LFP)]./lfp_fs;

% create the spike vector

ymin=floor(min(LFP));
ymax=ceil(max(LFP));

if ymin>=ymax
	ymax=ymin+1;
end

% set plot number according to user selection

switch lower(spikevis)
	case 'ticks'
		nplots=5;
	case 'ifr'
		nplots=5;
	case 'both'
		nplots=6;
	otherwise
		error('Did not understand spikevis argument (options are ticks, ifr ,or both)');
end

% use the first three subplots positions for the field

ax=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LFP PLOTTING%%%%%%%%%%%%%%%%%%%%%%%%

ax(end+1)=subaxis(nplots,1,1,1,1,3,'margin',.15,'spacingvert',0);
plot(timevec_lfp,LFP,'b-','linewidth',1.15);
set(gca,'XTick',[],'YTick',[ymin ymax],'TickDir','out','TickLength',[.02 .02],'xcolor',get(fignum,'color'));
ylim([ymin-eps ymax+eps]);
ylabel('LFP Amp. ($\mu$V)','fontsize',13,'fontname','helvetica','interpreter','latex')
box off
graphcount=4;

% spike ticks, might consider color coding according to phase

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SPIKE PLOTTING%%%%%%%%%%%%%%%%%%%%%%

if strcmp(spikevis,'ticks') || strcmp(spikevis,'both')

	spikes=SPIKETIMES(:)'./spike_fs;

	% plot spikes by simply drawing lines between 0 and 1 at the spike times

	spikevec=[spikes;spikes];
	trialvec=[ones(size(spikes));zeros(size(spikes))];
	ax(end+1)=subaxis(nplots,1,1,graphcount,1,1,'margin',.15,'spacingvert',0);
	plot(spikevec,trialvec,'k-');
	axis off
	graphcount=graphcount+1;

end

% IFR plotting, could also color code according to phase

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% IFR PLOTTING %%%%%%%%%%%%%%%%%%%%%%%

if strcmp(spikevis,'ifr') || strcmp(spikevis,'both')

	timevec_ifr=[1:length(IFR)]./ifr_fs;
	
	ymin=floor(min(IFR));
	ymax=ceil(max(IFR));

	if ymin>=ymax
		ymin=0;
		ymax=100;
	end

	ax(end+1)=subaxis(nplots,1,1,graphcount,1,1,'margin',.15,'spacingvert',.05);
	stairs(timevec_ifr,IFR,'r-','linewidth',1.5);
	set(gca,'YTick',[ymin ymax],'TickDir','out','TickLength',[.02 .02],'FontSize',13,'FontName','Helvetica','layer','top');
	ylim([ymin-eps ymax+eps]);
	ylabel('IFR (impulses/sec)','fontsize',13,'fontname','helvetica');
	box off
	xlabel('Time (s)','FontSize',13,'FontName','Helvetica')
	graphcount=graphcount+1;

end

linkaxes(ax,'x');

if length(time_range)~=2
	xlim([timevec_lfp(1) timevec_lfp(end)]);
else
	xlim([time_range(1) time_range(2)]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

