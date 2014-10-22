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
bound=.6;
% remove eps generation, too slow here...
fbstart='';
fr_smoothing=200;
fs=500;
bin_trials=50;
clip=5;

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TIME SERIES PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

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

	fr_smoothed.y.raw=smooth(time_elapsed,fr_mu,fr_smoothing,'rlowess');
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

	bound_smp=round(bound*FR(i).rms.fs);
	fr_mu=mean(FR(i).rms.time_series(bound_smp:end-bound_smp,:));

	rms_smooth=smooth(time_elapsed,fr_mu,fr_smoothing,'rlowess');

	plot(time_elapsed,fr_mu,'ro','markerfacecolor','r','markersize',4);
	hold on;
	plot(time_elapsed,rms_smooth,'k-','linewidth',2);
	
	set(gca,'TickDir','in','TickLength',[.02 .02],'FontSize',9,'FontName','Helvetica');

	ylim([prctile(fr_mu,clip) prctile(fr_mu,100-clip)]);

	%ephys_ratetracker_plot(FR(i).rms.song,FR(i).datenums,'fig_num',panel_fig);
	
	smooth_trace.timescale=[smooth_trace.timescale fr_smoothing];
	smooth_trace.x=[smooth_trace.x time_elapsed];
	smooth_trace.y.raw=[smooth_trace.y.raw rms_smooth];
	smooth_trace.y.norm=[smooth_trace.y.norm (rms_smooth-min(rms_smooth))./(max(rms_smooth)-min(rms_smooth))];

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% HISTOGRAM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

hist_fig=figure('visible','off');

ncols=4;
nrows=ceil(nchannels/ncols);
ntrials=length(FR(1).datenums);
if bin_trials>ntrials/2
	bin_trials=floor(ntrials/2);
end

for i=1:nchannels

	bound_smp=round(bound*FR(i).rms.fs);
	fr_mu=mean(FR(i).rms.time_series(bound_smp:end-bound_smp,:));

	% take the first and last 50 trials, plot as histograms


	% first two hours for pre
	%
	
	% last two hours for post
	%

	startpt=bin_trials;
	stoppt=ntrials-bin_trials;

	predata=fr_mu(1:startpt);
	postdata=fr_mu(stoppt:end);

	fr_min=floor(prctile([predata(:);postdata(:)],clip));
	fr_max=ceil(prctile([predata(:);postdata(:)],100-clip));

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
	%xlims=xlim();
	%xlims=round(xlims./10).*10;
	xlims=[fr_min fr_max];
	xlim([xlims]);

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% HISTOGRAM (RMS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

hist_fig2=figure('visible','off');

ncols=4;
nrows=ceil(nchannels/ncols);

for i=1:nchannels

	bound_smp=round(bound*FR(i).rms.fs);
	fr_mu=mean(FR(i).rms.time_series(bound_smp:end-bound_smp,:));

	% take the first and last 50 trials, plot as histograms

	ntrials=length(FR(i).datenums);

	% first two hours for pre
	%
	
	% last two hours for post
	%

	startpt=bin_trials;
	stoppt=ntrials-bin_trials;

	predata=fr_mu(1:startpt);
	postdata=fr_mu(stoppt:end);

	fr_min=floor(prctile([predata(:);postdata(:)],clip));
	fr_max=ceil(prctile([predata(:);postdata(:)],100-clip));

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
	%xlims=round(xlims);
	xlims(1)=floor(xlims(1));
	xlims(2)=ceil(xlims(2));
	xlim(xlims);

	ylims=ylim();
	ylims=round(ylims.*10)./10;
	ylims=[0 ylims(2)+.1];
	ylim(ylims);

	set(gca,'XTick',xlims,'YTick',ylims,'FontSize',9);
	title(['CH ' num2str(channel_labels(i))]);

	if i==1
		xlabel('RMS');
		ylabel('P');
	end

	set(hist_fig2,'position',[0 600 100*ncols 100*nrows],'PaperPositionMode','auto');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% STATS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%tightfig;

stats.distmat.trend.raw=1-squareform(pdist(smooth_trace.y.raw','correlation'));
stats.distmat.trend.norm=1-squareform(pdist(smooth_trace.y.norm','correlation'));

stats.distmat.ts.raw=zeros(nchannels,nchannels);
stats.distmat.ts.norm=zeros(nchannels,nchannels);

disp('Computing TS correlation matrix');
[nblanks formatstring]=fb_progressbar(100);
fprintf(1,['Progress:  ' blanks(nblanks)]);

% TODO: look at correlation between RMS change (pre-post) and time-series correlation

counter=1;
for i=1:nchannels
	for j=1:nchannels

		fprintf(1,formatstring,round((counter/(nchannels^2))*100));
		counter=counter+1;

		bound_smp=round(bound*FR(i).rms.fs);
		fr1=mean(zscore(FR(i).rms.time_series(bound_smp:end-bound_smp,:)),2);
		fr2=mean(zscore(FR(j).rms.time_series(bound_smp:end-bound_smp,:)),2);	

		rho=corrcoef(fr1,fr2);
		%mu_rho_norm=mean(diag(rho));
		%stats.distmat.ts.norm(i,j)=mu_rho_norm;
		stats.distmat.ts.norm(i,j)=rho(2,1);

	end
end

fprintf(1,'\n');

stats.prepost.change=zeros(1,nchannels);
stats.channels=channel_labels;

for i=1:nchannels
	
	fr_mu=mean(FR(i).rms.time_series(bound_smp:end-bound_smp,:));

	ntrials=length(FR(i).datenums);

	startpt=bin_trials;
	stoppt=ntrials-bin_trials;

	predata=fr_mu(1:startpt);
	postdata=fr_mu(stoppt:end);


	f=fit(time_elapsed(:),fr_mu(:),'poly1','robust','LAR'); % robust linear regression on the data
	f_coeffs=coeffvalues(f);

	stats.prepost.change(i)=(mean(postdata)-mean(predata))./sqrt(std(postdata)*std(predata));
	stats.prepost.change(i)=f_coeffs(1);


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CORRELATION FIGURES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% get trial-averaged correlation coefficient between trials (do they covary across trials?)

% compute correlation of average timeseries

corr_fig=figure('visible','off');

ch_labels=cat(1,FR(:).channel);

imagesc(ch_labels,ch_labels,stats.distmat.trend.raw);
title('FR trend correlation by channel (Pearson)');
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
title('FR ts correlation by channel (Pearson)');
set(gca,'FontSize',10,'FontName','Helvetica',...
	'TickLength',[0 0],'xtick',ch_labels(:),'ytick',ch_labels(:));
cb=colorbar('south');
axis square;

set(cb,'position',[.15 .01 .7 .03]);
set(cb,'TickDir','in','TickLength',[.01 .01],'FontSize',8)

set(corr_fig2,'position',[0 600 300 300],'PaperPositionMode','auto');

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
	multi_fig_save(hist_fig2,savedir,'mua_hist_rms','eps,png,fig','res',100);
	close([hist_fig2]);


	save(fullfile(savedir,'cat_mua.mat'),'smooth_trace','stats');

end




%
