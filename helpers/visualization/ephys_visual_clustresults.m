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
	0 .5 0;... % green	
	1 0 1;... % magenta	
	0 1 1;... % cyan		
	1 0 0;... % red		
	1 1 0;... % yellow
	0 0 1;... % blue		
	.6445 .1641 .1641; ... % brown
	1 .6445 0;... % orange		
	0 0 0;... % black
	.1953 .8008 .1953;... % lime-green
	0 0 .5;... % navy
	0 .5 .5;... % teal	
	.5 0 .5; ... % purple		
	.9542 .6406 .3750;... % sandy brown
	1 .2695 0;... % orange red	
	.5 .5 0;... % olive	
	.5 .5 .5;... % silver	
	.1328 .5430 .1328;... % forest green
	.4 .4 .4;... % gray	
	.2930 0 .5078;... % indigo	
	.1797 .5430 .3398;... % sea green
	.5430 0 .5430;... % dark magenta
	.2 .2 .2;... % dark gray		
	.6 .6 .6;... % lighter
	.7188 .5234 .0430;... % dark goldenrod
	.5938 .9833 .5938;... % pale green
	]; 



for i=1:2:nparams
	switch lower(varargin{i})
		case 'spike_fs'
			spike_fs=varargin{i+1};
		case 'fig_num'
			fig_num=varargin{i+1};
		case 'colors'
			colors=varargin{i+1};
	end
end

[nclust]=length(SPIKEWINDOWS);

[samples,trials]=size(SPIKEWINDOWS{1});

timevec=([1:samples]./spike_fs)*1e3;

if isempty(fig_num)
	fig_num=figure('Visible','on');
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

for i=1:nclust
	plot(timevec,SPIKEWINDOWS{plotorder(i)},'-','color',colors(i,:));
	hold on;
end

ylabel('Voltage (in $\mu$V)','FontName','Helvetica','FontSize',13,'Interpreter','Latex');
xlabel('Time (ms)','FontName','Helvetica','FontSize',13);
box off
axis tight;

prettify_axis(gca,'FontSize',12,'FontName','Helvetica');
prettify_axislabels(gca,'FontSize',15,'FontName','Helvetica');
