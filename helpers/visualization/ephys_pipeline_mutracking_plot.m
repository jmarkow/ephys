function ephys_pipeline_mutracking_plot(DATA,DATENUMS,varargin)
%
%
%
%
%
%

% creates a scatter plot of multi-unit activity against time with a smooth average

% take the data, square and take the average

intan_fs=25e3;
freq_range=[300 2e3];
plotcolor=[.7 .7 .7];
windowsize=15;
mads=30;
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

if isempty(grpid)
	grpid=ones(size(DATA,2));
end

% bandpass filter to recover the multi-unit data

[b,a]=ellip(5,.2,40,freq_range./(intan_fs/2),'bandpass');
filtdata=filtfilt(b,a,double(DATA));

% square and take the average (any other measures to consider here?  maybe std instead?)

plotdata=sqrt(mean(filtdata.^2));
%plotdata=mean(filtdata.^2);
%plotdata=mean(abs(hilbert(filtdata)));
% consider making the x axis time and not trials
% clean up data in sliding window, remove anything >=3MADS+MEDIAN

todel=[];
for i=1:length(plotdata)-windowsize

	winmed=median(plotdata(i:i+windowsize-1))
	winmad=mad(plotdata(i:i+windowsize-1))

	todel=[todel find((plotdata>winmed+mads*winmad)|(plotdata<winmed-mads*winmad))];
end

plotdata(todel)=[];
grpid(todel)=[];
DATENUMS(todel)=[];

ylimits=[ min(plotdata)-20 max(plotdata)+20 ];

fig=figure();

%get number of groups from grpid

grps=unique(grpid)

% shade the baseline (if it has been provided)

startdate=datestr(DATENUMS(startid))
times=zeros(size(DATENUMS));

for i=1:length(times)
	times(i)=etime(datevec(DATENUMS(i)),datevec(DATENUMS(startid)))/(60*60*24);
end

if any(grps==1)

	basemean=mean(plotdata(grpid==1));
	basevar=std(plotdata(grpid==1));

	upperconf=basemean+basevar;
	lowerconf=basemean-basevar;

	xvec=[min(times) max(times)];
	xdata=[xvec fliplr(xvec)];
	ydata=[lowerconf.*ones(size(xvec)) upperconf.*ones(size(xvec))];

	patch(xdata,ydata,1,'facecolor',[.75 .75 .75],'edgecolor','none');
	hold on;
	plot(xvec,[basemean.*ones(size(xvec))],'k--','linewidth',2)
end

% set the x vector relative to the startid

for i=1:length(grps)
	ax(i)=plot(times(grpid==grps(i)),plotdata(grpid==grps(i)),...
		'.','color',grpcolors(grps(i),:),'markersize',10);
	hold on;
	%plot(times(grpid==grps(i)),smooth(plotdata(grpid==grps(i)),20),'color',[0 0 0],'linewidth',1.5);
end

box off;
axis tight;
ylabel('MU activity','FontSize',20,'FontName','Helvetica');
xlabel('Days since feedback on','FontSize',20,'FontName','Helvetica');
set(gca,'layer','top','linewidth',2,'ticklength',[.025 .025],...
	'tickdir','out','fontsize',20,'fontname','helvetica');


L=legend(ax,grplabels(grps));
legend boxoff;
set(L,'FontSize',20,'FontName','Helvetica','Location','NorthWest');
%ch=get(L,'children');
%set(ch(1:3:end),'markersize',30);


