function ephys_evalisi(ISI,varargin)
%
%
%
%
%
%

fs=25e3; % default Intan fs
isipoints=[0:.01:25];
colors={'m','g','y','b'};

nparams=length(varargin);

if mod(nparams,2)>0
	error('ephysPipeline:argChk','Parameters must be specified as parameter/value pairs!');
end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'fs'
			fs=varargin{i+1};
		case 'isipoints'
			isipoints=varargin{i+1};
	end
end

if ~iscell(ISI)
	tmp=ISI;
	clear ISI;
	ISI{1}=tmp;
end


for i=1:length(ISI)
	[density,xi]=ksdensity((ISI{i}/fs)*1e3,isipoints,'support','positive');
	plot(xi,density,'-','color',colors{i},'linewidth',3);
	hold on;
end

hline=findobj(gca,'type','line');
set(hline,'clipping','off');
box off
%set(h,'FaceColor',[.7 .7 .7],'EdgeColor','k','LineWidth',1.5);
xlabel('ISI (ms)');
ylabel('Probability density');

ylimits(1)=min(density);
ylimits(2)=ceil(max(density)*100)/100;

if ylimits(1)<ylimits(2)
	set(gca,'YLim',ylimits,'YTick',ylimits);
end

prettify_axis(gca,'FontSize',12,'FontName','Helvetica');
prettify_axislabels(gca,'FontSize',15,'FontName','Helvetica');
xlim([xi(1) xi(end)]);


