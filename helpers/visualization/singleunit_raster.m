function singleunit_raster(TIME,HISTOGRAM,CLUSTER,PROC_DATA,varargin)
%
%
%
%
%


nparams=length(varargin);

if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs');
end

savedir=pwd;
savefilename_raster='raster';
ntrials=[];
spikeheight=.3;
channels=[];
figtitle='';
min_f=1;
max_f=10e3;
hist_colors='jet';
singletrials=5;

for i=1:2:nparams
	switch lower(varargin{i})
		case 'spikeheight'
			spikeheight=varargin{i+1};
		case 'savedir'
			savedir=varargin{i+1};
		case 'savefilename_raster'
			savefilename_raster=varargin{i+1};
		case 'ntrials'
			ntrials=varargin{i+1};
		case 'channels'
			channels=varargin{i+1};
		case 'figtitle'
			figtitle=varargin{i+1};
		case 'min_f'
			min_f=varargin{i+1};
		case 'max_f'
			max_f=varargin{i+1};
		case 'hist_colors'
			hist_colors=varargin{i+1};
		case 'singletrials'
			singletrials=varargin{i+1};
	end
end

% how many trials to show here?


nplots=length(CLUSTER.windows);

if isempty(ntrials)

	trialmax=-inf;

	for i=1:nplots
		tmpmax=max(CLUSTER.trials{i});
		if tmpmax>trialmax
			trialmax=tmpmax;
		end
	end

	ntrials=trialmax;
end

% form plot vectors/matrices

for j=1:nplots

	spike_vec{j}=[];
	trial_vec{j}=[];

	for k=1:ntrials

		clusterspikes=CLUSTER.times{j}(CLUSTER.trials{j}==k);
		clusterspikes=clusterspikes./CLUSTER.parameters.fs;

		tmp=repmat(clusterspikes,2,1);
		spike_vec{j}=[spike_vec{j} tmp];

		tmp=[k+spikeheight;k-spikeheight];
		tmp=repmat(tmp,1,length(clusterspikes));
		trial_vec{j}=[trial_vec{j} tmp];

		% for storage resort the spike vector by cluster

		clustspiketimes{j}{k}=clusterspikes;

	end
end

if ~isempty(HISTOGRAM)
	startidx=max([find(HISTOGRAM.f<=min_f)]);

	if isempty(startidx)
		startidx=1;
	end

	stopidx=min([find(HISTOGRAM.f>=max_f)]);

	if isempty(stopidx)
		stopidx=length(HISTOGRAM.f);
	end
end

for j=1:nplots

	raster_fig=figure('visible','off','Units','Pixels','Position',[0 0 600 800]);

	if ~isempty(HISTOGRAM)
		ax(1)=subaxis(6,1,1,1,1,2,'margin',.1,'spacingvert',0);
		imagesc(HISTOGRAM.t,HISTOGRAM.f(startidx:stopidx),HISTOGRAM.imask(startidx:stopidx,:));
		colormap(hist_colors);
		freezeColors;
		set(gca,'ydir','normal','tickdir','out','xtick',[],'ytick',[min_f max_f],'linewidth',1.5,'ticklength',[.025 .025],...
			'FontSize',11,'FontName','Helvetica');
		ylim([min_f max_f]);
		box off;
		title({[ figtitle ];[ 'Channel ' num2str(channels)]},'FontSize',18,'FontName','Helvetica');

		ax(2)=subaxis(6,1,1,3,1,1,'spacingvert',0,'margin',0.1,'paddingbottom',.025);
		plot(TIME,HISTOGRAM.mean_osc,'-k');
		ylabel('Osc.','FontSize',13,'FontName','Helvetica');
		axis tight;
		set(gca,'tickdir','out','xtick',[],'ytick',[]);

		ax(3)=subaxis(6,1,1,4,1,3,'spacingvert',0.025,'margin',0.1,'paddingbottom',0);

		plot(spike_vec{j},trial_vec{j},'-','color','k');	
		axis([0 TIME(end) 0 ntrials+1]);
		xlabel('Time (in s)','FontSize',13,'FontName','Helvetica');
		ylabel('Trial','FontSize',13,'FontName','Helvetica');
		box off
		set(gca,'tickdir','out','linewidth',1.5,'ticklength',[.025 .025],'FontSize',11,'FontName','Helvetica','ydir','rev');

		linkaxes(ax,'x');

	else
		plot(spike_vec{j},trial_vec{j},'-','color','k');
		axis([0 TIME(end) 0 ntrials+1]);

		xlabel('Time (in s)','FontSize',13,'FontName','Helvetica');
		ylabel('Trial','FontSize',13,'FontName','Helvetica');
		box off
		set(gca,'tickdir','out','linewidth',1.5,'ticklength',[.025 .025],'FontSize',11,'FontName','Helvetica','ydir','rev');

	end

	if ~isempty(figtitle)
		name=figtitle;
	end

	set(raster_fig,'PaperPositionMode','auto');
	multi_fig_save(raster_fig,savedir,...
		[ savefilename_raster num2str(j) ],'eps,png');
	close([raster_fig]);

	% generate single trial plots, simply the spikes marked along with the filtered traces

	randpopulation=randperm(ntrials);

	if singletrials<=ntrials
		randtrials=randpopulation(1:singletrials);
	else
		randtrials=randpopulation(1:ntrials);
		singletrials=ntrials;
	end

	singletrialdir=fullfile(savedir,'singletrials',['ch' num2str(channels) ],[ 'clust' num2str(j)]);

	if ~isempty(CLUSTER.windows)

		if ~exist(fullfile(singletrialdir),'dir');
			mkdir(singletrialdir);
		end

		disp('Plotting single trials');

		for k=1:singletrials

			singletrialfig=figure('Visible','off');
			idx=randtrials(k);
			currspikes=CLUSTER.times{j}(CLUSTER.trials{j}==idx);
			plot(TIME,PROC_DATA(:,idx,1),'k-');
			hold on
			plot(currspikes./CLUSTER.parameters.fs,PROC_DATA(currspikes,idx,1),'r*');
			set(gca,'TickDir','out','TickLength',[.02 .02],'FontName','Helvetica','FontSize',11);
			xlabel('Time (s)','FontName','Helvetica','FontSize',13,'interpreter','latex');
			ylabel('Voltage ($\mu$V)','FontName','Helvetica','FontSize',13,'interpreter','latex');
			box off
			axis tight
			multi_fig_save(singletrialfig,singletrialdir,['trial' num2str(idx)],'png','res',150);

			close([singletrialfig]);

		end
	end
end
