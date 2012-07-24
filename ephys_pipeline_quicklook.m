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
%
% see also ephys_condition_data.m,ephys_denoise_data.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% TODO finish documenting


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETER COLLECTION  %%%%%%%%%%%%%%


[samples,nchannels]=size(EPHYS_DATA);

channels=[1];
colors='jet';
sr=25e3;
t=[1:length(MIC_DATA)]./sr;
min_f=1e3;
max_f=8e3;
exclude=[];
figtitle=[];
figsave=[];
noise='none';
car_exclude=[];
freq_range=[];
n=800;
overlap=750;
nfft=1024;
type='s';

% high pass the mic trace

[b,a]=butter(4,1e3/(sr/2),'high');
MIC_DATA=filtfilt(b,a,MIC_DATA);

nparams=length(varargin);

if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs');
end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'channels'
			channels=varargin{i+1};
		case 'colors'
			colors=varargin{i+1};
		case 'sr'
			sr=varargin{i+1};
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
		otherwise
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SIGNAL CONDITIONING  %%%%%%%%%%%%%%%


if isempty(nfft)
	nfft=2^nextpow2(n);
end

% denoise according to user's preference

plot_data=ephys_denoise_signal(EPHYS_DATA,CHANNELS,channels,'method',noise,'car_exclude',car_exclude);

% condition the signal according to data type

if ~isempty(freq_range)
	plot_data=ephys_condition_signal(plot_data,type,'freq_range',freq_range);
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

[spect,f,t2]=pretty_sonogram(MIC_DATA,sr,'n',n,'overlap',overlap,'nfft',nfft);

startidx=max([find(f<=min_f)]);
stopidx=min([find(f>=max_f)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTTING CODE  %%%%%%%%%%%%%%%%%%%%%



% 3 plots for sound + n channels

nplots=3+length(channels);

quickfig=figure('Visible','off','Units','Pixels','Position',[0 0 800 1.2e3]);

[spect,f,t2]=pretty_sonogram(MIC_DATA,sr,'n',n,'overlap',overlap,'nfft',nfft);

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
box off;

ax(2)=subaxis(nplots,1,2,'spacingvert',0,'margin',0.1,'paddingbottom',.025);
plot(t,MIC_DATA,'-k');
ylabel('Osc.');
set(gca,'tickdir','out','xtick',[],'ytick',[],'xcolor',get(quickfig,'color'));

for i=1:length(channels)

	ax(2+i)=subaxis(nplots,1,2+i,'spacingvert',0.025,'margin',0.1,'paddingbottom',0);
	plot(t,plot_data(:,i)); % took out CAR 6/20/12

	if i<length(channels)
		set(gca,'xtick',[],'xcolor',get(quickfig,'color'));
	end

	% add a label to the right

	set(gca,'tickdir','out');
	box off;
	axis tight;
	
	%ylabel('V (in microvolts)');

	rightax(i)=axes('position',get(gca,'Position'),'color','none',...
		'xtick',[],'ytick',[],'yaxislocation','right','box','off');
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

set(quickfig,'visible','on');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
