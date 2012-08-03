function fig_num=multi_unit_raster(HISTOGRAM,RASTER,varargin)
%creates a multi-unit raster
%
%
%

nparams=length(varargin);

if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs');
end



fs=25e3;
min_f=1;
max_f=10e3;
hist_colors='jet';
raster_colors='hot';
fig_num=[];

for i=1:2:nparams
	switch lower(varargin{i})
		case 'fs'
			fs=varargin{i+1};
		case 'min_f'
			min_f=varargin{i+1};
		case 'max_f'
			max_f=varargin{i+1};
		case 'raster_colors'
			raster_colors=varargin{i+1};
		case 'hist_colors'
			hist_colors=varargin{i+1};
		case 'fig_title'
			fig_title=varargin{i+1};
		case 'fig_num'
			fig_num=varargin{i+1};
		case 'yaxis'
			yaxis=vararign{i+1};
	end
end

if isempty(fig_num)
	fig_num=figure();
end

startidx=max([find(HISTOGRAM.f<=min_f)]);

if isempty(startidx)
	startidx=1;
end

stopidx=min([find(HISTOGRAM.f>=max_f)]);

if isempty(stopidx)
	stopidx=length(HISTOGRAM.f);
end

time=[1:length(HISTOGRAM.mean_osc)]./fs;

ax(1)=subaxis(6,1,1,1,1,2,'margin',.1,'spacingvert',0);
imagesc(HISTOGRAM.t,HISTOGRAM.f(startidx:stopidx),HISTOGRAM.imask(startidx:stopidx,:));
colormap(hist_colors);
freezeColors;
set(gca,'ydir','normal','tickdir','out','xtick',[],'ytick',[min_f max_f],'linewidth',1.5,'ticklength',[.025 .025],...
	'FontSize',16,'FontName','Helvetica');
ylim([min_f max_f]);
%ylabel('Hz','FontSize',13,'FontName','Helvetica');
box off;
title([fig_title],'FontSize',18,'FontName','Helvetica');

ax(2)=subaxis(6,1,1,3,1,1,'spacingvert',0,'margin',0.1,'paddingbottom',.025);
plot(time,HISTOGRAM.mean_osc,'-k');
ylabel('Osc.','FontSize',20,'FontName','Helvetica');
axis tight;
set(gca,'tickdir','out','xtick',[],'ytick',[]);

% any trials to exclude?

ax(3)=subaxis(6,1,1,4,1,3,'spacingvert',0.025,'margin',0.1,'paddingbottom',0);
imagesc(RASTER.t,RASTER.trials,RASTER.image);
axis tight;
colormap(raster_colors);
freezeColors;
%plot(t,data(:,channelvis(i))-mean(data,2));
xlabel('Time (in s)','FontSize',13,'FontName','Helvetica');
ylabel('Trial','FontSize',13,'FontName','Helvetica');
box off
set(gca,'tickdir','out','linewidth',1.5,'ticklength',[.025 .025],'FontSize',11,'FontName','Helvetica');

linkaxes(ax,'x');



