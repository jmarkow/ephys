function ephys_pipeline_mutracking_plot_trials(DATA,varargin)
%
%
%
%
%
%

% creates a scatter plot of multi-unit activity against time with a smooth average

% take the data, square and take the average

intan_fs=25e3;
freq_range=[300 6e3];
plotcolor=[.7 .7 .7];
windowsize=15;
mads=20;
grpid=[];
grpcolors=[0 0 0;...
	0 0 1;...
	1 0 0];
startid=1; % which datenum indicates feedback start?
grplabels={'No FB','Escape down','Escape up'};
nparams=length(varargin);

if mod(nparams,2)>0
	error('ephysPipeline:argChk','Parameters must be specified as parameter/value pairs!');
end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'grpid'
			grpid=varargin{i+1};
		case 'startid'
			startid=varargin{i+1};
	end
end


% bandpass filter to recover the multi-unit data

[b,a]=ellip(5,.2,40,freq_range./(intan_fs/2),'bandpass');
filtdata=filtfilt(b,a,double(DATA)).^2;

% square and take the average (any other measures to consider here?  maybe std instead?)

plotdata=sqrt(mean(filtdata));
%plotdata=mean(abs(hilbert(filtdata)));
% consider making the x axis time and not trials
% clean up data in sliding window, remove anything >=3MADS+MEDIAN

if isempty(grpid)
	grpid=ones(1,length(plotdata));
end


todel=[];
for i=1:length(plotdata)-windowsize

	winmed=median(plotdata(i:i+windowsize-1));
	winmad=mad(plotdata(i:i+windowsize-1));
	todel=[todel find((plotdata>winmed+mads*winmad)|(plotdata<winmed-mads*winmad))];
end


plotdata(todel)=[];
grpid(todel)=[];

ylimits=[ min(plotdata)-20 max(plotdata)+20 ];

fig=figure();

%get number of groups from grpid

grps=unique(grpid);

% set the x vector relative to the startid

if any(grps==1)

	basemean=mean(plotdata(grpid==1));
	basevar=std(plotdata(grpid==1));

	upperconf=basemean+basevar;
	lowerconf=basemean-basevar;

	xvec=[1 length(plotdata)];
	xdata=[xvec fliplr(xvec)];
	ydata=[lowerconf.*ones(size(xvec)) upperconf.*ones(size(xvec))];

	patch(xdata,ydata,1,'facecolor',[.75 .75 .75],'edgecolor','none');
	hold on;
	plot(xvec,[basemean.*ones(size(xvec))],'k--','linewidth',2)
end


for i=1:length(grps)
	ax(i)=plot(find(grpid==grps(i)),plotdata(grpid==grps(i)),...
		'.','color',grpcolors(grps(i),:),'markersize',10);
	hold on;
	%plot(smooth(plotdata,40),'color',[0 0 0],'linewidth',1.5);
end

plot(smooth(plotdata,20),'color',[0 0 0],'linewidth',2);
box off;
axis tight;
ylabel('MU activity','FontSize',20,'FontName','Helvetica');
xlabel('Trial','FontSize',20,'FontName','Helvetica');
set(gca,'layer','top','linewidth',2,'ticklength',[.025 .025],...
	'tickdir','out','fontsize',20,'fontname','helvetica');

L=legend(ax,grplabels(grps));
legend boxoff;
set(L,'FontSize',20,'FontName','Helvetica','Location','NorthWest');
ch=get(L,'children');
set(ch(1:3:end),'markersize',30);
