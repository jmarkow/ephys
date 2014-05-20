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

fig=figure();

%get number of groups from grpid

grps=unique(grpid)

% shade the baseline (if it has been provided)

startdate=datestr(DATENUMS(startid))
times=zeros(size(DATENUMS));

for i=1:length(times)
	times(i)=etime(datevec(DATENUMS(i)),datevec(DATENUMS(startid)))/(60*60*24);
end

% set the x vector relative to the startid

for i=1:length(grps)
	ax(i)=plot(times(grpid==grps(i)),DATA(grpid==grps(i)),...
		'.','color',grpcolors(grps(i),:),'markersize',10);
	hold on;
	%plot(times(grpid==grps(i)),smooth(DATA(grpid==grps(i)),20),'color',[0 0 0],'linewidth',1.5);
end

box off;
axis tight;
ylabel('MU activity','FontSize',20,'FontName','Helvetica');
xlabel('Days','FontSize',20,'FontName','Helvetica');
set(gca,'layer','top','linewidth',2,'ticklength',[.025 .025],...
	'tickdir','out','fontsize',20,'fontname','helvetica');


%ch=get(L,'children');
%set(ch(1:3:end),'markersize',30);


