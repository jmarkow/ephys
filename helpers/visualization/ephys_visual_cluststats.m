function fig_num=ephys_visual_cluststats(SPIKEWINDOWS,SPIKETIMES,SPIKEDATA,varargin)
%cluster statistics, include Fisher projection and other quality metrics
%
%
%

%
% perhaps include plot of stability across trials (maybe threshold?)
%

nparams=length(varargin);

if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs');
end

fig_num=[];
maxlag=.2;
xres=.00075;
spacingvert=.15;
spacinghor=.15;

colors=[...	
	1 .6445 0;... % orange	
	1 0 1;... % magenta
	0 1 1;... % cyan	
	1 0 0;... % red				
	0 0 1;... % blue		
	.5 0 .5; ... % purple		
	.6445 .1641 .1641; ... % brown
	1 1 0;... % yellow
	.1953 .8008 .1953;... % lime-green
	.1328 .5430 .1328;... % forest green
	0 0 .5;... % navy
	0 .5 .5;... % teal
	.5430 .2695 .0742;... % saddle-brown
	1 .2695 0;... % orange red
	]; 

crosscolor=[0 0 0];
spike_fs=200e3;

% just in case add the hot colormap at the end

stats=[];

for i=1:2:nparams
	switch lower(varargin{i})
		case 'fig_num'
			fig_num=varargin{i+1};
		case 'maxlag'
			maxlag=varargin{i+1};
		case 'xres'
			xres=varargin{i+1};
		case 'spike_fs'
			spike_fs=varargin{i+1};
		case 'stats'
			stats=varargin{i+1};
	end
end

colorappend=colormap('hot');
colors=[colors;colorappend];

% fisher LDA for each cluster
% now compare all combinations

K=length(SPIKETIMES);

multi_clust=0;
if K>1
	clustercombos=nchoosek(1:K,2);
	ncombos=size(clustercombos,1);
	clustercombos=sortrows(clustercombos,[1 2]);
	multi_clust=1;
end

if isempty(fig_num)
	fig_num=figure('Visible','on','position',[0 0 200+300*K 200+300*K],'renderer','painters');
end

% grid for figure as follows:
% top left 2x2 axes (mean wave)
% top right KxK axes (fisher proj)
% bottom left Kx2 axes (autocorr)
% bottom right KxK axes (crosscorr)

nrows=max(2,K-1)+K-1;
ncols=2+K-1;

left_spacinghor=.15.*(1/(1+log(K)));
left_spacingvert=.05.*(1/(1+log(K)));
left_margin=.15.*(1/(1+log(K)));

legends={};
%if ~isempty(stats)
%	for i=1:length(stats.isod)
%		legends{i}=[ 'L ' sprintf('%.2f',stats.lratio(i)) ...
%			' ID ' sprintf('%.0f',stats.isod(i)) ];
%	end
%end

spikeheight=max(1,K-2);

subaxis(nrows,ncols,1,1,2,spikeheight,'margin',left_margin,'spacingvert',left_spacingvert,'spacinghor',left_spacinghor);
ephys_visual_clustresults(SPIKEWINDOWS,'spike_fs',spike_fs,'fig_num',fig_num,'legend_labels',legends);
xlabel('');
ylabel('Amplitude ($\mu$Volts)','interpreter','latex');
set(gca,'ticklength',[0 0],'linewidth',1);

% move to the top right for the Fisher plot
% top right, Fisher projection

row=0;
col=2;

right_margin=.05;
right_spacingvert=(1/K).*spacingvert;
right_spacinghor=(1/K).*spacinghor;

if multi_clust
	for i=1:ncombos

		[density1,density2,xi]=fisher_projection(SPIKEDATA{clustercombos(i,1)},...
			SPIKEDATA{clustercombos(i,2)});

		xcoord= [ xi fliplr(xi) ];
		ycoord1=[zeros(size(xi)) fliplr(density1)];
		ycoord2=[zeros(size(xi)) fliplr(density2)];

		subaxis(nrows,ncols,clustercombos(i,2)+col-1,clustercombos(i,1)+row,1,1,...
			'margin',right_margin,'spacingvert',right_spacingvert,'spacinghor',right_spacinghor);

		% on a Mac transparency is completely borked :(

		patch(xcoord,ycoord1,1,'facecolor',colors(clustercombos(i,1),:),...
			'edgecolor','none');
		hold on
		patch(xcoord,ycoord2,1,'facecolor',colors(clustercombos(i,2),:),...
			'edgecolor','none');

		prettify_axis(gca,'FontSize',12,'FontName','Helvetica','linewidth',1);
		set(gca,'layer','top','linewidth',1,'ticklength',[0 0]);
		box off;
		axis tight;

		if i==1
			ylabel('P');
			xlabel('Fisher Score');
		end

		prettify_axislabels(gca,'FontSize',12,'FontName','Helvetica','linewidth',1);
		axis off;


	end
end

% bottom left, autocorr

row=spikeheight;
col=0;

for i=1:K

	subaxis(nrows,ncols,1+col,i+row,2,1,...
		'margin',left_margin,'spacingvert',left_spacingvert,'spacinghor',left_spacinghor);

	get_spike_correlogram(SPIKETIMES{i},SPIKETIMES{i},...
		'fig_num',fig_num,'type','auto','maxlag',maxlag,'xres',xres,...
		'color',colors(i,:));

	prettify_axis(gca,'FontSize',12,'FontName','Helvetica','linewidth',1);
	set(gca,'xtick',[-maxlag*1e3:50:maxlag*1e3],'ticklength',[0 0],'layer','top');

	if i==K
		xlabel('Lag (ms)');
		ylabel('Autocorr (Hz)');
	else
		set(gca,'xtick',[]);
	end

	box off
	axis tight
	prettify_axislabels(gca,'FontSize',12,'FontName','Helvetica');
	
end

row=K-1;
col=2;
spacinghor=.025;

if multi_clust
	for i=1:ncombos


		subaxis(nrows,ncols,clustercombos(i,2)+col-1,clustercombos(i,1)+row,1,1,...
			'margin',right_margin,'spacingvert',right_spacingvert,'spacinghor',right_spacinghor);

		get_spike_correlogram(SPIKETIMES{clustercombos(i,1)},SPIKETIMES{clustercombos(i,2)},...
			'fig_num',fig_num,'type','cross','maxlag',maxlag,'xres',xres,...
			'color',crosscolor);

		prettify_axis(gca,'FontSize',12,'FontName','Helvetica','linewidth',1);
		set(gca,'xtick',[-maxlag*1e3 maxlag*1e3],'ticklength',[0 0],'layer','top');

		if i==1
			xlabel('Lag (ms)');
			ylabel('Crosscorr (Hz)');
		else
			set(gca,'xtick',[]);
		end

		prettify_axislabels(gca,'FontSize',12,'FontName','Helvetica');

		% shift zero ytick up

		box off
		axis tight


	end
end


