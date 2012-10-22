function fig_num=auto_correlogram(SPIKETIMES,varargin)
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

% just in case add the hot colormap at the end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'fig_num'
			fig_num=varargin{i+1};
		case 'maxlag'
			maxlag=varargin{i+1};
		case 'xres'
			xres=varargin{i+1};
	end
end

colorappend=colormap('hot');
colors=[colors;colorappend];

% fisher LDA for each cluster
% now compare all combinations

K=length(SPIKETIMES);

clustercombos=nchoosek(1:K,2);
ncombos=size(clustercombos,1);

if isempty(fig_num)
	fig_num=figure('Visible','on','position',[0 0 500 200*K]);
end

for i=1:K

	subaxis(K,1,1,i,'margin',.15,'spacingvert',.05);
	get_spike_correlogram(SPIKETIMES{i},SPIKETIMES{i},...
		'fig_num',fig_num,'type','auto','maxlag',maxlag,'xres',xres,...
		'color',colors(i,:));

	prettify_axis(gca,'FontSize',17,'FontName','Helvetica');
	set(gca,'xtick',[-maxlag*1e3 maxlag*1e3],'ticklength',[0 0],'layer','top','linewidth',2);

	if i==K
		xlabel('Lag (ms)');
		ylabel('Autocorr (Hz)');
	else
		set(gca,'xtick',[]);
	end

	% shift zero ytick up

	box off
	axis tight

end

