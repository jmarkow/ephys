function ephys_visual_simpleraster(CLUSTER,CLUSTERNUM,varargin)
%for cluster quality visualization, plots mean/median waveforms against variance
%
%	ephys_visual_waveplot(WINDOWS,varargin)
%
%	WINDOWS
%	samples x trials matrix of waveforms (loaded from sua_channels x.mat, stored in clusterwindows)
%
%	the following may be passed as parameter/value pairs:
%	
%		fs
%		sampling frequency of spikes (normally twice Intan sampling rate, spikes are interpolated by default)
%
%		snr
%		snr to display as a title
%
%
%
%

nparams=length(varargin);
if mod(nparams,2)>0
	error('ephysPipeline:argChk','Parameters must be specified as parameter/value pairs!');
end

%%%
spike_height=.5;
max_time=[];


for i=1:2:nparams
	switch lower(varargin{i})
		case 'spike_height'
			spike_height=varargin{i+1};	
		case 'max_time'
			max_time=varargin{i+1};
		otherwise

	end
end

spikefig=figure();

% grab the spikes from a particular cluster and triasl

trials=min(CLUSTER.trials{CLUSTERNUM}):max(CLUSTER.trials{CLUSTERNUM});

maxt=-inf;
for i=trials
	
	% grab the spike times for each trial

	newspikes{i}=CLUSTER.times{CLUSTERNUM}(CLUSTER.trials{CLUSTERNUM}==i)./CLUSTER.parameters.fs;

	maxspike=max(newspikes{i});
	if maxspike>maxt
		maxt=maxspike;
	end
end

newtrialvec=[];
newspikevec=[];

% smooth ifr with Gaussian window per Hahnloser

for i=1:length(newspikes)
	newspikevec=[newspikevec [ newspikes{i};newspikes{i}] ];
	newtrialvec=[newtrialvec [ ones(1,length(newspikes{i})).*i+spike_height; ones(1,length(newspikes{i})).*i-spike_height] ];
end

if isempty(max_time)
	max_time=maxt;
end

fig=figure();

plot(newspikevec,newtrialvec,'-','color','k');

xlabel({'Time (in s)'},'FontSize',18,'FontName','Helvetica');
ylabel('Trial','interpreter','latex','FontSize',18,'FontName','Helvetica');


xlim([0 max_time]);
set(gca,'tickdir','out','FontSize',15,'FontName','Helvetica','linewidth',1.25,'TickLength',[.025 .025],'ydir','reverse');
box off


