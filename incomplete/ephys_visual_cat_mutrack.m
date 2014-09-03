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

	end
end

if ~isempty(savedir)
	savedir=fullfile(savedir,'mua');
end

% rate plot for each channel as a function of TOD

nchannels=length(FR);
channel_labels=cat(1,FR(:).channel);

% for each channel, generate plot
%
%

% each channel is a column

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
	ephys_ratetracker_plot(fr_mu,FR(i).datenums,'fig_num',panel_fig);

	% could also bin and show as heatmap across the array...
	%
	% need to analyze song duration, LFP, etc. etc.

	title(['CH ' num2str(channel_labels(i))]);
	set(gca,'XTick',[]);	
	
	if i==1
		ylabel('FR');
	end

	ax(end+1)=subplot(2,nchannels,i+nchannels);
	tmp=ephys_ratetracker_plot(FR(i).rms.song,FR(i).datenums,'fig_num',panel_fig);
	
	smooth_trace.timescale=[smooth_trace.timescale tmp.timescale(:)];
	smooth_trace.x=[smooth_trace.x tmp.x(:)];
	smooth_trace.y.raw=[smooth_trace.y.raw tmp.y.raw(:)];
	smooth_trace.y.norm=[smooth_trace.y.norm tmp.y.norm(:)];

	if i==1
		ylabel('RMS');
		xlabel('Time (days)');
	end

	% set to an appropriate size, save
	
	set(panel_fig,'position',[ 0 600 1700 200],'PaperPositionMode','auto');
	
end

linkaxes(ax,'x');
tightfig;

hist_fig=figure('visible','off');

ncols=4;
nrows=ceil(nchannels/ncols);

for i=1:nchannels

	bound_smp=round(bound*FR(i).fs);
	fr_mu=mean(FR(i).time_series(bound_smp:end-bound_smp,:));

	% take the first and last 50 trials, plot as histograms

	ntrials=length(FR(i).datenums);
	
	predata=fr_mu(1:min(50,ntrials));
	postdata=fr_mu(end-min(50,ntrials-1):end);

	fr_min=prctile([predata(:);postdata(:)],0)-10;
	fr_max=prctile([predata(:);postdata(:)],100)+10;

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

% generate correlation matrices...
% should be equivalent, good to compute both either way

stats.distmat.raw=1-squareform(pdist(smooth_trace.y.raw','correlation'));
stats.distmat.norm=1-squareform(pdist(smooth_trace.y.norm','correlation'));

corr_fig=figure('visible','off');

ch_labels=cat(1,FR(:).channel);

imagesc(ch_labels,ch_labels,stats.distmat.raw);
title('FR time-course correlation by channel (Pearson)');
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
	multi_fig_save(corr_fig,savedir,'mua_corr','eps,png,fig','res',100);
	close([corr_fig]);
	multi_fig_save(hist_fig,savedir,'mua_hist','eps,png,fig','res',100);
	close([hist_fig]);

	save(fullfile(savedir,'cat_mua.mat'),'smooth_trace','stats');

end




%
