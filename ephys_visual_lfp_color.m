function [mean_lfp lfp_channels proc_ref_data]=ephys_visual_lfpcolor(EPHYS_DATA,CHANNELS,varargin)
%generates a panel with average LFPs
%
%	ephys_visual_lfpcolor(EPHYS_DATA,HISTOGRAM,CHANNELS,varargin)
%
%	EPHYS_DATA
%	sound-aligned voltage traces from extracted_data.mat (should be the variable ephys_data)
%
%	CHANNELS
%	channel labels (i.e. the channel that corresponds to a given element in the cell array
%	ephys_data) from extracted_data.mat%
%
%	the following may be specified as parameter/value pairs:
%
%		fs
%		sample rate (default: 25e3)
%
%		proc_fs
%		processing sampling rate (default: 1e3, downsamples data to this frequency)
%
%

% 5 panels, one for each electrode, color-coded, last panel simply plot the fields on top of each other


%%%%%%%%%%%%%%%%%%% 

% gather field data, zero-phase filter, collect the firing rates

lfp_channels=CHANNELS;
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
use_ave=0; % take an average of the lfp_channels and reference
           % against the average out?

%%%%%%%%%%%%%%%%%%% Parameter collection

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
		case 'use_ave'
			use_ave=varargin{i+1};
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if length(spike_channel)>1
	error('ephysPipeline:lfpcolor:toomanyspikechannels',...
		'More than one spike channel not supported');
end


downfact=fs/proc_fs;

% anti-alias

[b,a]=butter(3,[100/(25e3/2)],'low');

[nsamples,ntrials,nchannels]=size(EPHYS_DATA);

%%%%%% delete any non-existent channels 

to_del=[];
for i=1:length(lfp_channels)
	if ~any(lfp_channels(i)==CHANNELS)
		to_del=[to_del i];
	end
end

lfp_channels([to_del])=[];

to_del=[];
for i=1:length(reref_channels)
	if ~any(reref_channels(i)==CHANNELS)
		to_del=[to_del i];
	end
end

reref_channels([to_del])=[];

% compute the reference if it's defined
proc_ref_data=[];

if ~isempty(reref_channels)

	disp(['Rereferencing with channels ' num2str(reref_channels) ]);
	ref_data=zeros(nsamples,ntrials);

	for i=1:length(reref_channels)
		ref_data=ref_data+EPHYS_DATA(:,:,find(CHANNELS==reref_channels(i)))./length(reref_channels);
	end

	% for visualization plot the mean of the reference

end

proc_data=zeros(nsamples,ntrials,length(lfp_channels));

if use_ave

	for i=1:length(lfp_channels)
		proc_data(:,:,i)=EPHYS_DATA(:,:,find(CHANNELS==lfp_channels(i)));
	end

	proc_data=mean(proc_data,3);

	if ~isempty(reref_channels)
		proc_data=proc_data-ref_data;
	end

	proc_data=filtfilt(b,a,double(proc_data));
	lfp_channels=0;

else

	if ~isempty(reref_channels)
		EPHYS_DATA(:,:,i)=EPHYS_DATA(:,:,i)-ref_data;
	end

	for i=1:length(lfp_channels)
		proc_data(:,:,i)=filtfilt(b,a,double(EPHYS_DATA(:,:,find(CHANNELS==lfp_channels(i)))));
	end
end

nplots=length(lfp_channels);

if ~isempty(reref_channels)
	
	ref_data=downsample(filtfilt(b,a,double(ref_data)),downfact);

	proc_ref_data=ephys_condition_signal(ref_data,'l','freq_range',freq_range,'medfilt_scale',medfilt_scale,'medfilt',0,...
		'fs',proc_fs,'filt_order',filt_order);
	
	ref_mean=mean(proc_ref_data,2);
	nplots=nplots+1;

end

% downsample

clear EPHYS_DATA;
proc_data=downsample(proc_data,downfact);

% filter, median filter, demean, detrend...

proc_data=ephys_condition_signal(proc_data,'l','freq_range',freq_range,'medfilt_scale',medfilt_scale,'medfilt',0,...
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

for i=1:nchannels
	mean_lfp(:,i)=mean(proc_data(:,:,i),2);
end

if isempty(fig_num)
	fig_num=figure('Visible','on');
end


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


%%%%%%%%%%%%%%%%%%% Plotting code

min_ylim=inf;
max_ylim=-inf;

for i=1:nchannels

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

% add a plot if we're rereferencing

if ~isempty(reref_channels)
	
	ax(length(lfp_channels+1))=subaxis(nplots,1,length(lfp_channels)+1);
	plot([1:nsamples]./proc_fs,ref_mean,'linewidth',1.5);

	ylimits=ylim();
	set(gca,'TickDir','out','ytick',[]);

	% where to put the ylabel...right side?

	ylabel({'REF';...
		sprintf('%.0f',ylimits(1));...
		sprintf('%.0f',ylimits(2))});
	
	axis tight;

end

linkaxes(ax,'x');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
