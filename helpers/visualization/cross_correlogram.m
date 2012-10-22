function fig_num=cross_correlogram(SPIKETIMES,varargin)
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
color=[0 0 0];

for i=1:2:nparams
	switch lower(varargin{i})
		case 'fig_num'
			fig_num=varargin{i+1};
		case 'maxlag'
			maxlag=varargin{i+1};
		case 'xres'
			xres=varargin{i+1};
		case 'color'
			color=varargin{i+1};
	end
end

% fisher LDA for each cluster
% now compare all combinations

K=length(SPIKETIMES);

clustercombos=nchoosek(1:K,2);
ncombos=size(clustercombos,1);

if isempty(fig_num)
	fig_num=figure('Visible','on','position',[0 0 200*K 200*K]);
end

for i=1:ncombos

	subaxis(K,K,clustercombos(i,2),clustercombos(i,1),'margin',.1,'spacingvert',.05,'spacinghor',.05);

	get_spike_correlogram(SPIKETIMES{clustercombos(i,1)},SPIKETIMES{clustercombos(i,2)},...
		'fig_num',fig_num,'type','cross','maxlag',maxlag,'xres',xres,...
		'color',color);

	prettify_axis(gca,'FontSize',17,'FontName','Helvetica');
	set(gca,'xtick',[-maxlag*1e3 maxlag*1e3],'ticklength',[0 0],'layer','top','linewidth',2);

	if i==1
		xlabel('Lag (ms)');
		ylabel('Crosscorr (Hz)');
	else
		set(gca,'xtick',[]);
	end

	box off
	axis tight

end



