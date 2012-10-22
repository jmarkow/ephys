function fig_num=ephys_visual_clustresults(SPIKEWINDOWS,varargin)
%spike statistics figure, include 2D histogram and ISI distribution
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
spike_fs=50e3;

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

legend_labels=[];

for i=1:2:nparams
	switch lower(varargin{i})
		case 'spike_fs'
			spike_fs=varargin{i+1};
		case 'fig_num'
			fig_num=varargin{i+1};
		case 'colors'
			colors=varargin{i+1};
		case 'legend_labels'
			legend_labels=varargin{i+1};
	end
end

[nclust]=length(SPIKEWINDOWS);
[samples,trials]=size(SPIKEWINDOWS{1});

timevec=([1:samples]./spike_fs)*1e3;

if isempty(fig_num)
	fig_num=figure('Visible','on','renderer','painters');
end

% organize cluster by p2p

for i=1:nclust
	mean_waveform=mean(SPIKEWINDOWS{i},2);
	minpeak=min(mean_waveform);
	maxpeak=max(mean_waveform);
	p2p(i)=maxpeak-minpeak;
end

% organize waveforms by p2p

[val plotorder]=sort(p2p,'descend');

% instead simply plot mean and std, then color code everything
% and show--(1) waveforms, (2) fisher projections, (3) auto-correlations, and (4) cross-correlations

patchx=[ timevec fliplr(timevec) ];
hsvcolors=rgb2hsv(colors);

for i=plotorder

	meanwave=mean(SPIKEWINDOWS{i},2)';
	varwave=std(SPIKEWINDOWS{i},0,2)';
	patchy=[ meanwave-varwave fliplr(meanwave+varwave) ];

	patchcolor=hsv2rgb(hsvcolors(i,:));
	edgecolor=hsv2rgb(hsvcolors(i,:)-[0 0 .3]);
	linecolor=hsv2rgb(hsvcolors(i,:)-[0 0 .4]);

	h(i)=patch(patchx,patchy,1,'facecolor',...
		patchcolor,'edgecolor',edgecolor,'facealpha',1);

	hold on;
	plot(timevec,meanwave,'-','color',linecolor,'linewidth',2);

end

ylabel('microVolts');
xlabel('Time (ms)');
box off
axis tight;

if ~isempty(legend_labels)
	L=legend(h,legend_labels);
	legend boxoff;
end

prettify_axis(gca,'FontSize',17,'FontName','Helvetica');
prettify_axislabels(gca,'FontSize',15,'FontName','Helvetica');

if ~isempty(legend_labels)
	set(L,'location','SouthEast','FontSize',10);
end

set(gca,'layer','top');
