function [H,P]=pairedsample_plot(SAMPLE1,SAMPLE2,varargin)
%generates a paired sample plot overlaid on a simple boxplot, also
%uses MATLAB's ttest implementation to compute a paired sample t-test
%
%
%

nparams=length(varargin);

tail='right';
showbox=1;
boxwidth=.25;
labels{1}='Pre-treatment';
labels{2}='Post-treatment';

if mod(nparams,2)>0
	error('ephysPipeline:argChk','Parameters must be specified as parameter/value pairs!');
end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'tail'
			tail=varargin{i+1};
		case 'showbox'
			showbox=varargin{i+1};
		case 'boxwidth'
			boxwidth=varargin{i+1};
		case 'labels'
			labels=varargin{i+1};
	end
end

if ~isvector(SAMPLE1) || ~isvector(SAMPLE2)
	error('ephysPipeline:pairedsamplesplot:inputnotvector','Both inputs must be vectors.');
end

SAMPLE1=SAMPLE1(:);
SAMPLE2=SAMPLE2(:);

if length(SAMPLE1)~=length(SAMPLE2)
	error('ephysPipeline:pairedsampleplot:unequalsamples','Sample sizes must be the same.');
end

% set up spectrogram parameters

pairedfig=figure('Visible','off');

% draw box indicating either 25-75th percentile or 95% confidence (1.96*SEM)
% draw line indicating the median

if showbox

	extent= [ 1-boxwidth/2 1+boxwidth/2 ];

	% first plot the box

	lowerconf=prctile(SAMPLE1,25).*ones(1,2);
	upperconf=prctile(SAMPLE1,75).*ones(1,2);

	x= [ extent fliplr(extent) ];
	y= [ lowerconf upperconf ];

	patch(x,y,1,'edgecolor','none','facecolor',[ 0 .7461 1]);
	hold on;

	% show the median

	plot(extent,median(SAMPLE1).*ones(1,2),'-','linewidth',3,'color',[1 .6445 0])

	%%%%% sample 2

	extent= [ 2-boxwidth/2 2+boxwidth/2 ];

	% first plot the box

	lowerconf=prctile(SAMPLE2,25).*ones(1,2);
	upperconf=prctile(SAMPLE2,75).*ones(1,2);

	x= [ extent fliplr(extent) ];
	y= [ lowerconf upperconf ];

	patch(x,y,1,'edgecolor','none','facecolor',[ 0 .7461 1]);
	hold on;

	% show the median

	plot(extent,median(SAMPLE2).*ones(1,2),'-','linewidth',3,'color',[1 .6445 0])

end

% plot the data points connected by lines


plot([ SAMPLE1 SAMPLE2 ]','m-o','color','m','linewidth',1.05,'markersize',10,...
	'markerfacecolor','m','markeredgecolor','m')

maxval=max([SAMPLE1(:);SAMPLE2(:)]);

axis([ .5 2.5 0 maxval+eps]);

set(gca,'XTick',[1 2],'XTickLabel',labels);
ylabel('Impedance ($\Omega$)','interpreter','latex','FontSize',15,'FontName','Helvetica');
prettify_axis(gca,'FontSize',15,'FontName','Helvetica','ticklength',[.03 .03]);
prettify_axislabels(gca,'FontSize',15,'FontName','Helvetica');

[H,P]=ttest(SAMPLE1,SAMPLE2,.01,tail);

set(pairedfig,'Visible','on')


