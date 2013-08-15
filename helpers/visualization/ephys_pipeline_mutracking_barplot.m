function [groupdata,fig,stats]=ephys_pipeline_mutracking_barplot(DATA,varargin)
%
%
%
%
%
%

% creates a scatter plot of multi-unit activity against time with a smooth average

% take the data, square and take the average

stats=[];
intan_fs=25e3;
plotcolor=[.7 .7 .7];
windowsize=15;
mads=30;
grpid=[];
grpcolors=[0 0 0;...
	0 0 1;...
	1 0 0];
startid=1; % which datenum indicates feedback start?

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
		case 'edges'
			edges=varargin{i+1};
	end
end

% bandpass filter to recover the multi-unit data


% for each group get the first and last 100-200 trials, histogram and run simple statistics (ranksum maybe)

groups=unique(grpid);

for i=1:length(groups)
	groupdata{i}=[];
end

for i=2:length(edges)

	currgrp=edges(i-1):edges(i);
	grpidx=grpid(edges(i-1));

	% take the first and last 100, or percentage

	grpidx

	len=length(edges(i-1):edges(i));

	control=DATA(edges(i-1):edges(i-1)+len/10);
	obs=DATA(edges(i)-len/10:edges(i));
	
	mucontrol=mean(control);
	stdcontrol=std(control);

	control=(control-mucontrol)./stdcontrol;
	obs=(obs-mucontrol)./stdcontrol;

	[p,h]=ranksum(control,obs)

	groupdata{1}=[groupdata{1} control];
	groupdata{grpidx}=[groupdata{grpidx} obs];
end

bins=[-8:.5:8];
fig=figure();
for i=1:length(groupdata)
	
	[N,BINVEC]=pretty_histogram(groupdata{i},bins)
	N=N./sum(N);
	subplot(1,2,1);
	stairs(BINVEC,N,'color',grpcolors(i,:),'linewidth',2);
	hold on;

end
box off;
set(gca,'layer','top','linewidth',1.5,'ticklength',[.025 .025],'tickdir','out',...
	'fontsize',20,'fontname','helvetica');
%%%%% bar chart

% plot means and 3*sem

barheights=zeros(1,length(groupdata));
barerr=zeros(1,length(groupdata));

for i=1:length(groupdata)
	barheights(i)=mean(groupdata{i});
	barerr(i)=2*std(groupdata{i})./sqrt(length(groupdata{i}));
end


for i=1:length(groupdata)
	subplot(1,2,2);
	x=i-1+.75;
	h=bar(x,barheights(i),1);
	baseh=get(h,'BaseLine')
	set(baseh,'linewidth',2)
	uistack(baseh,'top')
	hold on;
	plot([x;x],[barheights(i)-barerr(i);barheights(i)+barerr(i)],'k-','linewidth',1.5)
	set(h,'FaceColor','none','EdgeColor',grpcolors(i,:),'linewidth',3);
	set(gca,'xcolor',get(fig,'color'))
	
end
xlim([0 length(groupdata)+.5]);

box off;
set(gca,'layer','top','linewidth',1.5,'ticklength',[.025 .025],...
	'tickdir','out','fontsize',20,'fontname','helvetica','xtick',[]);
ylabel('Power');



