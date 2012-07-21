function fig_num=time_frequency_raster(HISTOGRAM,TFIMAGE,varargin)
%time frequency raster plots
%
%
%
%

nparams=length(varargin);

if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs');
end

SR=25e3;
hist_min_f=1e3;
hist_max_f=10e3;
tf_min_f=1;
tf_max_f=100;
hist_colors='jet';
tfimage_colors='jet';
fig_num=[];
fig_title=[];
scale='log';
scalelabel=[];
colorbarsize=.02; % normalized units, height of the colorbar

for i=1:2:nparams
	switch lower(varargin{i})
		case 'sr'
			SR=varargin{i+1};
		case 'tf_min_f'
			tf_min_f=varargin{i+1};
		case 'tf_max_f'
			tf_max_f=varargin{i+1};
		case 'hist_min_f'
			hist_min_f=varargin{i+1};
		case 'hist_max_f'
			hist_max_f=varargin{i+1};
		case 'tfimage_colors'
			tfimage_colors=varargin{i+1};
		case 'hist_colors'
			hist_colors=varargin{i+1};
		case 'fig_title'
			fig_title=varargin{i+1};
		case 'fig_num'
			fig_num=varargin{i+1};
		case 'yaxis'
			yaxis=varargin{i+1};
		case 'scale'
			scale=varargin{i+1};
		case 'scalelabel'
			scalelabel=varargin{i+1};
	end
end

if isempty(scalelabel)
	switch lower(scale)
		case 'log'
			scalelabel='dB';
		case 'linear'
			scalelabel='Power';
	end
end

if isempty(fig_num)
	fig_num=figure();
end

hist_startidx=max([find(HISTOGRAM.f<=hist_min_f)]);

if isempty(hist_startidx)
	hist_startidx=1;
end

hist_stopidx=min([find(HISTOGRAM.f>=hist_max_f)]);

tf_startidx=max([find(TFIMAGE.f<=tf_min_f)]);

if isempty(tf_startidx)
	tf_startidx=1;
end

tf_stopidx=min([find(TFIMAGE.f>=tf_max_f)]);

time=[1:length(HISTOGRAM.mean_osc)]./SR;

ax(1)=subaxis(4,1,1,1,1,2,'margin',.12,'spacingvert',0,'paddingbottom',0);
imagesc(HISTOGRAM.t,HISTOGRAM.f(hist_startidx:hist_stopidx),HISTOGRAM.imask(hist_startidx:hist_stopidx,:));
axis tight;
colormap(hist_colors);
freezeColors;
set(gca,'ydir','normal','tickdir','out','xtick',[],'ytick',[hist_min_f hist_max_f],'linewidth',2,'ticklength',[.025 .025],...
	'FontSize',13,'FontName','Helvetica','xcolor',get(fig_num,'color'));
ylim([hist_min_f hist_max_f]);
%ylabel('Hz','FontSize',13,'FontName','Helvetica');
box off;
title([fig_title],'FontSize',18,'FontName','Helvetica');

ax(2)=subaxis(4,1,1,3,1,2,'spacingvert',0.05,'margin',0.12,'paddingbottom',0);

switch lower(scale)
	case 'log'
		disp('Log scale');
		imagesc(TFIMAGE.t,TFIMAGE.f(tf_startidx:tf_stopidx),...
			20*log10(TFIMAGE.image(tf_startidx:tf_stopidx,:)));
	case 'linear'
		disp('Linear scale');
		imagesc(TFIMAGE.t,TFIMAGE.f(tf_startidx:tf_stopidx),...
			TFIMAGE.image(tf_startidx:tf_stopidx,:));
	otherwise
		error('Did not understand the scale parameters (must be log or linear)');
end
		
axis tight;
colormap(tfimage_colors);
freezeColors;
xlabel('Time (in s)','FontSize',15,'FontName','Helvetica');
ylabel('Fs (Hz)','FontSize',15,'FontName','Helvetica');
box off;
set(gca,'ydir','normal','tickdir','out','linewidth',2,'ticklength',[.025 .025],'FontSize',13,'FontName','Helvetica');

% put a colorbar at the bottom

pos=get(ax(2),'pos');
set(ax(2),'pos',[pos(1) pos(2) pos(3)*.95 pos(4)]);
pos=get(ax(2),'pos');
hc2=colorbar('location','eastoutside','position',[pos(1)+pos(3)+.01 pos(2)+pos(4)/2 colorbarsize pos(4)/2]);
set(hc2,'linewidth',2,'FontSize',12,'FontName','Helvetica');
ylabel(hc2,scalelabel,'FontSize',15,'FontName','Helvetica');

% adjust the top axis accordingly

pos=get(ax(1),'pos');
set(ax(1),'pos',[pos(1) pos(2) pos(3)*.95 pos(4)]);

% make sure the time axes are synced

linkaxes(ax,'x');
if length(TFIMAGE.t)>1 & TFIMAGE.t(end)>TFIMAGE.t(1)
	xlim([TFIMAGE.t(1) TFIMAGE.t(end)]);
end
