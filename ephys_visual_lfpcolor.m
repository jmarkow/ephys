function ephys_visual_lfpcolor(EPHYS_DATA,CHANNELS,varargin)
%
%
%%


% 5 panels, one for each electrode, color-coded, last panel simply plot the fields on top of each other


%%%% 

% gather field data, zero-phase filter, collect the firing rates

lfp_channels=[6 9 13 15];
spike_channels=lfp_channels;
spike_clusters=[1 1 1 1];
fs=25e3;
proc_fs=1e3;
freq_range=[30 50];
filt_order=3;
medfilt_scale=1.5;
filedir=pwd;
sort_type='pipeline';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nparams=length(varargin);

if mod(nparams,2)>0
	error('ephysPipeline:argChk','Parameters must be specified as parameter/value pairs!');
end



for i=1:2:nparams
	switch lower(varargin{i})
		case 'fs'
			fs=varargin{i+1};
		case 'spike_channels'
			spike_channels=varargin{i+1};
		case 'lfp_channels'
			lfp_channels=varargin{i+1};
		case 'spike_clusters'
			spike_clusters=varargin{i+1};
	end
end




downfact=fs/proc_fs;

% anti-alias

[b,a]=butter(3,[200/(25e3/2)],'low');

[nsamples,ntrials,nchannels]=size(EPHYS_DATA);
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
spike_data=zeros(nsamples,ntrials,nchannels);
binedges=[1:nsamples]./proc_fs;

for i=1:length(spike_channels)

	load(fullfile(pwd,'sua',sort_type,[ 'sua_channels ' num2str(spike_channels(i)) '.mat']),'clust_spike_vec','subtrials');

	for j=subtrials
		spike_data(:,j,i)=ephys_ifr(clust_spike_vec{1}{spike_clusters(i)}{j}.*proc_fs,nsamples,proc_fs);
	end

end


% take the mean in each case, plot against firing rate of each neuron

% resynthesize based on contour image
proc_data=proc_data(:,subtrials,:);
[nsamples,ntrials,nchannels]=size(proc_data)
mean_lfp=zeros(nsamples,nchannels);

for i=1:length(lfp_channels)
	mean_lfp(:,i)=mean(zscore(proc_data(:,:,i)),2);
end

fig=figure('Visible','on');
nplots=length(lfp_channels)+1;

% get the firing rate of each neuron
%subplots(nplots,1,1);

for i=1:length(lfp_channels)

	ax(i)=subaxis(nplots,1,i,'spacingvert',0,'paddingbottom',0)
	
	x=[1:nsamples]./proc_fs;
	x=[x;x];

	y=mean_lfp(:,i)';
	y=[y;y];

	z=zeros(size(x));

	col=mean(spike_data(:,:,i),2)';
	col=[col;col];

	size(x)
	size(y)
	size(col)
	surface(x,y,z,col,'facecol','no','edgecol','interp','linew',2);

	%grays=colormap('gray');
	%grays(1:10,:)=[];

	ylabel([ 'LFP ' num2str(lfp_channels(i)) ' SU ' num2str(spike_channels(i))]);

	bones=colormap('bone');
	bones(1:5,:)=[];
	colormap(1-bones)
	%colormap(1-grays);

	axis tight;
	box off;

	set(gca,'TickDir','out','ytick',[]);
end

% plot the forcing functions together

ax(length(lfp_channels)+1)=subaxis(nplots,1,length(lfp_channels)+1);
colors=lines;

for i=1:length(lfp_channels)
	plot([1:nsamples]./proc_fs,mean_lfp(:,i),'k-','color',colors(i,:),'linewidth',1.2);
	hold on;
end

linkaxes(ax,'x');

