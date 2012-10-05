function ephys_quicklook(MIC_DATA,EPHYS_DATA,CHANNELS,varargin)
%generates a figure to glance at data extracted by the pipeline
%
%	ephys_quicklook(MIC_DATA,EPHYS_DATA,CHANNELS,varargin)
%
%	MIC_DATA
%	microphone time-series
%
%	EPHYS_DATA
%	ephys time-series aligned to microphone data
%
%	CHANNELS
%	channel IDs
%
%	the following can be passed as parameter/value pairs:
%
%		colors
%		colormap for spectrogram (string, default: jet)
%
%		channels
%		channels to display (default: 1)
%
%		min_f
%		minimum spectrogram FS for display (default: 1e3)
%
%		max_f
%		maximum spectrogram FS for display (default: 8e3)
%	
%		noise
%		noise reduction technique ('none','NN' for nearest neighbor,'CAR' for common average reference, default: 'none')
%
%		type
%		how to display data ('s' for single-unit,'m' for multi-unit,'l' for LFP, default: 's')
%
%		freq_range
%		frequency range for filtering (leave blank for auto-selection by display type, default: [])
%
%		n
%		spectrogrm Hanning window size (default: 800)
%
%		overlap
%		spectrogram overlap (default: 750)
%
%		nfft
%		spectrogram nfft (default: 1024)
%
%		ylim_match
%		if set to 1, matches all y axes to min and max, if set to a two element vector sets all axes
%		to min=ylim_match(1) and max=ylim_match(2), otherwise defaults to standard axes (default: 0)
%
%
% see also ephys_condition_data.m,ephys_denoise_data.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% TODO finish documenting


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETER COLLECTION  %%%%%%%%%%%%%%

if nargin<3
	error('ephysPipeline:quicklook:notenoughparams','Need 3 arguments to continue, see documentation');
end

[samples,nchannels]=size(EPHYS_DATA);

channels=[1];
colors='jet';
fs=25e3;
interp_fs=[];
t=[1:length(MIC_DATA)]./fs;
min_f=0;
max_f=8e3;
exclude=[];
figtitle=[];
figsave=[];
noise='none';
car_exclude=[];
freq_range=[800];
n=1024;
overlap=1e3;
nfft=1024;
type='s';
ylim_match=0;
filt_type='high';

% high pass the mic trace

[b,a]=butter(4,1e3/(fs/2),'high');
MIC_DATA=filtfilt(b,a,MIC_DATA);

nparams=length(varargin);

if mod(nparams,2)>0
	error('ephysPipeline:argChk','Parameters must be specified as parameter/value pairs!');
end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'channels'
			channels=varargin{i+1};
		case 'colors'
			colors=varargin{i+1};
		case 'fs'
			fs=varargin{i+1};
		case 'interp_fs'
			interp_fs=varargin{i+1};
		case 'n'
			n=varargin{i+1};
		case 'overlap'
			overlap=varargin{i+1};
		case 'highpass'
			highpass=varargin{i+1};
		case 'exclude'
			exclude=varargin{i+1};
		case 'figtitle'
			figtitle=varargin{i+1};
		case 'figsave'
			figsave=varargin{i+1};
		case 'noise'
			noise=varargin{i+1};
		case 'freq_range'
			freq_range=varargin{i+1};
		case 'type'
			type=varargin{i+1};
		case 'ylim_match'
			ylim_match=varargin{i+1};
		otherwise
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SIGNAL CONDITIONING  %%%%%%%%%%%%%%%

for i=1:length(channels)
	if isempty(find(channels(i)==CHANNELS))
		error('ephysPipeline:quicklook:suchanneldne','SU channel %g does not exist',channels(i));
	end
end

if isempty(nfft)
	nfft=2^nextpow2(n);
end

% denoise according to user's preference

plot_data=ephys_denoise_signal(EPHYS_DATA,CHANNELS,channels,'method',noise,'car_exclude',car_exclude);

% condition the signal according to data type

if ~isempty(freq_range)
	plot_data=ephys_condition_signal(plot_data,type,'freq_range',freq_range,'filt_type','high');
else
	plot_data=ephys_condition_signal(plot_data,type);
end

% delete any non-existent channels from the list

todel=[];
for i=1:length(channels)
	if ~any(channels(i)==CHANNELS)
		todel=[todel i];
	end	
end
channels(todel)=[];

% quick multi-taper spectrogram

[spect,f,t2]=pretty_sonogram(MIC_DATA,fs,'n',n,'overlap',overlap,'nfft',nfft);

startidx=max([find(f<=min_f)]);
stopidx=min([find(f>=max_f)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTTING CODE  %%%%%%%%%%%%%%%%%%%%%

% find min and max across channels

totalmin=min(plot_data(:));
totalmax=max(plot_data(:));

% 3 plots for sound + n channels

nplots=3+length(channels);

quickfig=figure('Visible','off','Units','Pixels','Position',[0 0 800 1.2e3]);

[spect,f,t2]=pretty_sonogram(MIC_DATA,fs,'n',n,'overlap',overlap,'nfft',nfft);

startidx=max([find(f<=min_f)]);
stopidx=min([find(f>=max_f)]);

% 3 plots for sound + n channels
nplot=3+length(channels);
if strcmp(lower(noise),'car')
	nplots=nplots+1;
end

quickfig=figure('Visible','off','units','pixels','position',[ 0 0 700 length(channels)*400]);

ax(1)=subaxis(nplots,1,1,'margin',.1,'spacingvert',0);
imagesc(t2,f(startidx:stopidx),spect(startidx:stopidx,:));
colormap(colors);
set(gca,'ydir','normal','tickdir','out','xtick',[],'ytick',[min_f max_f],'xcolor',get(quickfig,'color'));
ylim([min_f max_f]);
ylabel('Hz');
title([figtitle]);
prettify_axislabels(gca,'FontSize',15,'FontName','Helvetica');
prettify_axis(gca,'FontSize',15,'FontName','Helvetica','Linewidth',2,'TickLength',[.025 .025]);
box off;

ax(2)=subaxis(nplots,1,2,'spacingvert',0,'margin',0.1,'paddingbottom',.025);
plot(t,MIC_DATA,'-','color',[0 1 1]);
ylabel('Osc.');
prettify_axis(gca,'FontSize',15,'FontName','Helvetica','Linewidth',2,'TickLength',[.025 .025]);
prettify_axislabels(gca,'FontSize',15,'FontName','Helvetica');

set(gca,'tickdir','out','xtick',[],'ytick',[],'xcolor',get(quickfig,'color'));

if ~isempty(interp_fs)
	
	interpfactor=interp_fs/fs;

	disp(['Interpolating ephys vector to ' num2str(interp_fs) ' samples/s ']);

	% interpolation points

	oldt=t;
	t=linspace(t(1),t(end),length(t)*interpfactor);

	tmp=plot_data;
	clear plot_data;


	for i=1:length(channels)
		plot_data(:,i)=spline(oldt,tmp(:,i),t);
	end

end

for i=1:length(channels)

	ax(2+i)=subaxis(nplots,1,2+i,'spacingvert',0.025,'margin',0.1,'paddingbottom',0);
	plot(t,plot_data(:,i),'-','color',[0 .5430 .5430]); % took out CAR 6/20/12

	if i<length(channels)
		set(gca,'xtick',[],'xcolor',get(quickfig,'color'));
	end

	% add a label to the right

	ylabel('Voltage ($\mu$V)','interpreter','latex');
	prettify_axislabels(gca,'FontSize',15,'FontName','Helvetica');
	prettify_axis(gca,'FontSize',15,'FontName','Helvetica','Linewidth',2,'TickLength',[.025 .025]);

	box off;
	axis tight;
	
	%ylabel('V (in microvolts)');

	if length(ylim_match)>1
		ylim([ylim_match(1) ylim_match(2)]);
	elseif ylim_match
		ylim([totalmin totalmax]);
	end

	rightax(i)=axes('position',get(gca,'Position'),'color','none',...
		'xtick',[],'ytick',[],'yaxislocation','right','box','off','linewidth',2);
	prettify_axislabels(gca,'FontSize',15,'FontName','Helvetica');
	ylabel(rightax(i),[ 'Channel ' num2str(channels(i))])

end


if strcmp(lower(noise),'car')
    car=mean(plot_data,2);
	ax(end+1)=subaxis(nplots,1,3+i,'spacingvert',0.025,'margin',0.1);
	plot(t,car);
	rightax(i)=axes('position',get(gca,'Position'),'color','none',...
		'xtick',[],'ytick',[],'yaxislocation','right','box','off');

	ylabel('CAR');
end

xlabel(ax(end),'Time (in secs)');
linkaxes([ax rightax],'x');
xlim([t2(1) t2(end)]);
prettify_axislabels(gca,'FontSize',15,'FontName','Helvetica');

set(quickfig,'visible','on','PaperPositionMode','auto');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
