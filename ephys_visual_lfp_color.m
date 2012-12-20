function ephys_visual_lfpcolor(EPHYS_DATA,CHANNELS,varargin)
%
%
%%


% 5 panels, one for each electrode, color-coded, last panel simply plot the fields on top of each other


%%%% 

% gather field data, zero-phase filter, collect the firing rates

lfp_channels=[6 9 13 15];
spike_channel=[];
spike_cluster=[1];
fs=25e3;
proc_fs=1e3;
freq_range=[10 50];
filt_order=3;
medfilt_scale=1.5;
filedir=pwd;
sort_type='pipeline';
reref_channels=[];
fig_num=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nparams=length(varargin);

if mod(nparams,2)>0
	error('ephysPipeline:argChk','Parameters must be specified as parameter/value pairs!');
end



for i=1:2:nparams
	switch lower(varargin{i})
		case 'fs'
			fs=varargin{i+1};
		case 'spike_channel'
			spike_channel=varargin{i+1};
		case 'lfp_channels'
			lfp_channels=varargin{i+1};
		case 'spike_cluster'
			spike_cluster=varargin{i+1};
		case 'freq_range'
			freq_range=varargin{i+1};
		case 'reref_channels'
			reref_channels=varargin{i+1};
		case 'fig_num'
			fig_num=varargin{i+1};
	end
end

if length(spike_channel)>1
	error('ephysPipeline:lfpcolor:toomanyspikechannels',...
		'More than one spike channel not supported');
end

downfact=fs/proc_fs;

% anti-alias

[b,a]=butter(3,[200/(25e3/2)],'low');

[nsamples,ntrials,nchannels]=size(EPHYS_DATA);

to_del=[];
for i=1:length(lfp_channels)
	if ~any(lfp_channels(i)==CHANNELS)
		to_del=[to_del i];
	end
end

lfp_channels([to_del])=[];

if ~isempty(reref_channels)

	ref_data=zeros(nsamples,ntrials);

	for i=1:reref_channels	
		ref_data=ref_data+EPHYS_DATA(:,:,find(CHANNELS==reref_channels(i)))./length(reref_channels);
	end

	for i=1:nchannels
		EPHYS_DATA(:,:,i)=EPHYS_DATA(:,:,i)-ref_data;
	end

end

proc_data=zeros(nsamples,ntrials,length(lfp_channels));

for i=1:length(lfp_channels)
	proc_data(:,:,i)=filtfilt(b,a,double(EPHYS_DATA(:,:,find(CHANNELS==lfp_channels(i)))));
end

% downsample

clear EPHYS_DATA;
proc_data=downsample(proc_data,downfact);

% filter, median filter, demean, detrend...

proc_data=ephys_condition_signal(proc_data,'l','freq_range',freq_range,'medfilt_scale',medfilt_scale,'medfilt',1,...
	'fs',proc_fs,'filt_order',filt_order);

[nsamples,ntrials,nchannels]=size(proc_data);
spike_data=zeros(nsamples,ntrials);
binedges=[1:nsamples]./proc_fs;

if isempty(spike_channel)
	spike_data=[]
else
	load(fullfile(pwd,'sua',sort_type,[ 'sua_channels ' num2str(spike_channel) '.mat']),'clust_spike_vec','subtrials');

	for j=subtrials
		spike_data(:,j)=ephys_ifr(clust_spike_vec{1}{spike_cluster}{j}.*proc_fs,nsamples,proc_fs);
	end

	proc_data=proc_data(:,subtrials,:);

end

% take the mean in each case, plot against firing rate of each neuron

[nsamples,ntrials,nchannels]=size(proc_data)
mean_lfp=zeros(nsamples,nchannels);

for i=1:length(lfp_channels)
	mean_lfp(:,i)=mean(proc_data(:,:,i),2);
end

if isempty(fig_num)
	fig_num=figure('Visible','on');
end

nplots=length(lfp_channels);

% get the firing rate of each neuron
%subplots(nplots,1,1);

bones=bone;
bones(1:4,:)=[];

if isempty(spike_data)
	col=zeros(2,nsamples);
else
	col=mean(spike_data,2)';
	col=[col;col];
end

min_ylim=inf;
max_ylim=-inf;

for i=1:length(lfp_channels)

	ax(i)=subaxis(nplots,1,i,'spacingvert',0,'paddingbottom',0);
	
	x=[1:nsamples]./proc_fs;
	x=[x;x];

	y=mean_lfp(:,i)';
	y=[y;y];

	z=zeros(size(x));

	surface(x,y,z,col,'facecol','no','edgecol','interp','linew',2);
	colormap(1-bones)


	axis tight;
	box off;

	ylimits=ylim();
	set(gca,'TickDir','out','ytick',[]);

	% where to put the ylabel...right side?

	ylabel({[ 'CH' num2str(lfp_channels(i))];...
		sprintf('%.0f',ylimits(1));...
		sprintf('%.0f',ylimits(2))});

	if ylimits(1)<min_ylim
		min_ylim=ylimits(1);
	end

	if ylimits(2)>max_ylim
		max_ylim=ylimits(2);
	end

end

%for i=1:length(ax)
%	set(fig_num,'CurrentAxes',ax(i));
%	ylim([min_ylim max_ylim])
%end

%  plot the lfps together

%ax(length(lfp_channels)+1)=subaxis(nplots,1,length(lfp_channels)+1);
%colors=lines;
%
%%
%
%for i=1:length(lfp_channels)
%
%	plot([1:nsamples]./proc_fs,mean_lfp(:,i),'k-','color',colors(i,:),'linewidth',1.2);
%	hold on;
%end

linkaxes(ax,'x');

