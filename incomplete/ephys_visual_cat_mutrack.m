function [smooth_trace,stats]=ephys_ratetracker_plot(FR,varargin)
% generates a ratetracking plot given concatenated data
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETER COLLECTION %%%%%%%%%%%%%%%%%

% TODO:  remove multi-channel support

nparams=length(varargin);

if mod(nparams,2)>0
	error('ephysPipeline:argChk','Parameters must be specified as parameter/value pairs!');
end

savedir='ephys_cat';
figtitle='';

subtrials=[];
channels=[];
bound=.5;
% remove eps generation, too slow here...
fbstart='';
fr_smoothing=50;
fs=500;
bin_trials=200;
clip=1;

for i=1:2:nparams
	switch lower(varargin{i})
		case 'fr_smoothing'
			fr_smoothing=varargin{i+1};
		case 'bound'
			bound=varargin{i+1};
		case 'fbstart'
			fbstart=varargin{i+1};
		case 'fs'
			fs=varargin{i+1};
		case 'savedir'
			savedir=varargin{i+1};
		case 'bin_trials'
			bin_trials=varargin{i+1};
		case 'clip'
			clip=varargin{i+1};

	end
end

if ~isempty(savedir)
	savedir=fullfile(savedir,'mua');
end

fbstart_datevec=datevec((FR(1).datenums(1)));
fbstart_day=fbstart_datevec(3);

% convert date string into hours since midnight for display

time_elapsed=zeros(1,length(FR(1).datenums));

for i=1:length(time_elapsed)
	tmp=datevec(FR(1).datenums(i));
	time_elapsed(i)=etime(tmp,[fbstart_datevec(1:3) 0 0 0]);
end

time_elapsed=time_elapsed./(3600);

% rate plot for each channel as a function of TOD

nchannels=length(FR);
channel_labels=cat(1,FR(:).channel);

nplots=nchannels*2;
ax=[];

% generate smooth traces, useful for correlations, etc.

smooth_trace.timescale=[];
smooth_trace.x=[];
smooth_trace.y.raw=[];
smooth_trace.y.norm=[];

stats.distmat.raw=[];
stats.distmat.norm=[];

panel_fig=figure('visible','off');

for i=1:nchannels

	disp(['Processing channel ' num2str(channel_labels(i)) ]);

	% feed in the time-series and datenumbers	
	
	bound_smp=round(bound*FR(i).fs);	

	% plot average firing rate (consider other summary statistics)
	%
	
	fr_mu=mean(FR(i).time_series(bound_smp:end-bound_smp,:));
	ax(end+1)=subplot(2,nchannels,i);

	fr_smoothed.y.raw=smooth(time_elapsed,fr_mu,fr_smoothing,'lowess');
	fr_smoothed.x=time_elapsed;
	fr_smoothed.timescale=fr_smoothing;
	fr_smoothed.y.norm=(fr_smoothed.y.raw-min(fr_smoothed.y.raw))./(max(fr_smoothed.y.raw)-min(fr_smoothed.y.raw));

	plot(time_elapsed,fr_mu,'ro','markerfacecolor','r','markersize',4);
	hold on;
	plot(time_elapsed,fr_smoothed.y.raw,'k-','linewidth',2);
	
	set(gca,'TickDir','in','TickLength',[.02 .02],'FontSize',9,'FontName','Helvetica');

	ylim([prctile(fr_mu,clip) prctile(fr_mu,100-clip)]);

	%tmp=ephys_ratetracker_plot(fr_mu,FR(i).datenums,'fig_num',panel_fig);

	% could also bin and show as heatmap across the array...
	%
	% need to analyze song duration, LFP, etc. etc.

	title(['CH ' num2str(channel_labels(i))]);
	
	if i==1
		ylabel('FR');
	end

	ax(end+1)=subplot(2,nchannels,i+nchannels);

	rms_smooth=smooth(time_elapsed,FR(i).threshold,'lowess');

	plot(time_elapsed,FR(i).threshold,'ro','markerfacecolor','r','markersize',4);
	hold on;
	plot(time_elapsed,rms_smooth,'k-','linewidth',2);
	
	set(gca,'TickDir','in','TickLength',[.02 .02],'FontSize',9,'FontName','Helvetica');

	ylim([prctile(FR(i).rms.song,clip) prctile(FR(i).rms.song,100-clip)]);

	%ephys_ratetracker_plot(FR(i).rms.song,FR(i).datenums,'fig_num',panel_fig);
	
	smooth_trace.timescale=[smooth_trace.timescale fr_smoothed.timescale(:)];
	smooth_trace.x=[smooth_trace.x fr_smoothed.x(:)];
	smooth_trace.y.raw=[smooth_trace.y.raw fr_smoothed.y.raw(:)];
	smooth_trace.y.norm=[smooth_trace.y.norm fr_smoothed.y.norm(:)];

	if i==1
		ylabel('RMS');
		%xlabel('Time (days)');
	end

	% set to an appropriate size, save
	
	set(panel_fig,'position',[ 0 600 1700 200],'PaperPositionMode','auto');
	
end

linkaxes(ax,'x');
xlim([min(time_elapsed)-1 max(time_elapsed)+1]);

tightfig;

hist_fig=figure('visible','off');

ncols=4;
nrows=ceil(nchannels/ncols);

for i=1:nchannels

	bound_smp=round(bound*FR(i).fs);
	fr_mu=mean(FR(i).time_series(bound_smp:end-bound_smp,:));

	% take the first and last 50 trials, plot as histograms

	ntrials=length(FR(i).datenums);

	% first two hours for pre
	%
	
	% last two hours for post
	%

	startpt=min(find(time_elapsed>1));
	stoppt=max(find(time_elapsed<(max(time_elapsed)-1)));

	if startpt-1<bin_trials
		startpt=min(bin_trials,ntrials);
	end

	if ntrials-stoppt<bin_trials
		stoppt=ntrials-min(bin_trials,ntrials-1);
	end

	if isempty(startpt) | isempty(stoppt)
		startpt=min(bin_trials,ntrials);
		stoppt=ntrials-min(bin_trials,ntrials-1);
	end

	predata=fr_mu(1:startpt);
	postdata=fr_mu(stoppt:end);

	fr_min=prctile([predata(:);postdata(:)],clip)-10;
	fr_max=prctile([predata(:);postdata(:)],100-clip)+10;

	bins=linspace(fr_min,fr_max,15);

	[n_pre]=histc(predata,bins);
	[n_post]=histc(postdata,bins);

	n_pre=n_pre./sum(n_pre);
	n_post=n_post./sum(n_post);

	[x,y_pre]=ephys_stairplot(n_pre,bins);
	[~,y_post]=ephys_stairplot(n_post,bins);

	subplot(nrows,ncols,i);
	stairs(x,y_pre,'b-','linewidth',1.5);hold on;
	stairs(x,y_post,'r-','linewidth',1.5);hold on;
	box off;
	axis tight;
	xlims=xlim();
	xlims=round(xlims./10).*10;
	xlim(xlims);

	ylims=ylim();
	ylims=round(ylims.*10)./10;
	ylims=[0 ylims(2)+.1];
	ylim(ylims);

	set(gca,'XTick',xlims,'YTick',ylims,'FontSize',9);
	title(['CH ' num2str(channel_labels(i))]);

	if i==1
		xlabel('Firing rate');
		ylabel('P');
	end

	set(hist_fig,'position',[0 600 100*ncols 100*nrows],'PaperPositionMode','auto');
end
%tightfig;

stats.distmat.trend.raw=1-squareform(pdist(smooth_trace.y.raw','spearman'));
stats.distmat.trend.norm=1-squareform(pdist(smooth_trace.y.norm','spearman'));

stats.distmat.ts.raw=zeros(nchannels,nchannels);
stats.distmat.ts.norm=zeros(nchannels,nchannels);

disp('Computing TS correlation matrix');
[nblanks formatstring]=fb_progressbar(100);
fprintf(1,['Progress:  ' blanks(nblanks)]);

counter=1;
for i=1:nchannels
	for j=1:nchannels

		fprintf(1,formatstring,round((counter/(nchannels^2))*100));
		counter=counter+1;

		bound_smp=round(bound*FR(i).fs);
		fr1=FR(i).time_series(bound_smp:end-bound_smp,:);
		fr2=FR(j).time_series(bound_smp:end-bound_smp,:);

		rho=corr(fr1,fr2,'type','spearman');
		mu_rho_norm=mean(diag(rho));
		stats.distmat.ts.norm(i,j)=mu_rho_norm;

	end
end

fprintf(1,'\n');

% get trial-averaged correlation coefficient between trials (do they covary across trials?)


%stats.distmat.ts.raw=1-squareform(pdist(
%stats.distmat.ts.norm=
% compute correlation of average timeseries

corr_fig=figure('visible','off');

ch_labels=cat(1,FR(:).channel);

imagesc(ch_labels,ch_labels,stats.distmat.trend.raw);
title('FR trend correlation by channel (Spearman)');
set(gca,'FontSize',10,'FontName','Helvetica',...
	'TickLength',[0 0],'xtick',ch_labels(:),'ytick',ch_labels(:));
cb=colorbar('south');
axis square;

set(cb,'position',[.15 .01 .7 .03]);
set(cb,'TickDir','in','TickLength',[.01 .01],'FontSize',8)

set(corr_fig,'position',[0 600 300 300],'PaperPositionMode','auto');

corr_fig2=figure('visible','off');

ch_labels=cat(1,FR(:).channel);

imagesc(ch_labels,ch_labels,stats.distmat.ts.norm);
title('FR ts correlation by channel (Spearman)');
set(gca,'FontSize',10,'FontName','Helvetica',...
	'TickLength',[0 0],'xtick',ch_labels(:),'ytick',ch_labels(:));
cb=colorbar('south');
axis square;

set(cb,'position',[.15 .01 .7 .03]);
set(cb,'TickDir','in','TickLength',[.01 .01],'FontSize',8)

set(corr_fig,'position',[0 600 300 300],'PaperPositionMode','auto');

if ~isempty(savedir)

	mkdir(savedir);

	multi_fig_save(panel_fig,savedir,'mua_panels','eps,png,fig','res',100);
	close([panel_fig]);
	multi_fig_save(corr_fig,savedir,'mua_corr_trend','eps,png,fig','res',100);
	close([corr_fig]);
	multi_fig_save(corr_fig2,savedir,'mua_corr_ts','eps,png,fig','res',100);
	close([corr_fig2]);
	multi_fig_save(hist_fig,savedir,'mua_hist','eps,png,fig','res',100);
	close([hist_fig]);

	save(fullfile(savedir,'cat_mua.mat'),'smooth_trace','stats');

end




%
