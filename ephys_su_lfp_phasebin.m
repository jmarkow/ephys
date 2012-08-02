function [BINNED_PHASE,PHASE_MEAN]=ephys_su_lfp_phasebin(LFPCHANNEL,HISTOGRAM,varargin)
%computes coherograms between LFPs and single units
%
%	[BINNED_PHRASE PHASE_MEAN]=ephys_su_lfp_phasebin(LFPCHANNEL,HISTOGRAM,varargin)
%	
%	LFPCHANNEL
%	LFPCHANNEL to use
%
%	HISTOGRAM
%	time frequency histogram structure (result of ephys_visual_histogram or loading histogram.mat)
%
%	the following may be passed as parameter/value pairs:
%
%		savedir
%		directory to save results in (default: 'coherence')
%
%		filedir
%		directory that contains the single-unit data and LFP subdirectory (default: 'pwd')
%
%		colors
%		colormap for phase plot (default: '')
%
%		min_f
%		minimum frequency to display (default: 10)
%		
%		max_f
%		maximum frequency to display (default: 100)
%
%		freq_range
%		frequency range for LFP filtering (default: 300)
%
%		lfp_fs
%		lfp sampling rate (default: 25e3)
%
%		trial_range
%		only compute the coherence for these trials (leave blank for all, default: [])
%
%		medfilt_scale
%		median filter timescale in ms (default: 1.5)
%
%		phasebins
%		edges for phase bins
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETER COLLECTION %%%%%%%%%%%%%%%%%

if nargin<2
	error('ephysPipeline:phasebin:notenoughparams','Need 2 arguments to continue, see documentation');
end

nparams=length(varargin);

savedir=pwd;
filedir=pwd;
fig_title='noname';
phase_colors=[1 0 0;1 .84 0;0 1 0;0 0 1]; % let's do a custom colormap with colors matched to nphasebins
hist_colors='jet';
window_size=25; % in samples how large a window size do we want for the non-overlapping average?

ntapers=[];
freq_range=[10 20 40 50]; % let's just refilter
lfp_fs=25e3; % default Intan sampling rate
trial_range=[];
medfilt_scale=1.5; % median filter scale (in ms)
phase_bins=[-inf pi/4 2*pi/4 3*pi/4 inf]; % define phase bin edges (radians)
filt_name='kaiser';
min_f=1;
max_f=10e3;
colorbarsize=.02; % normalized units, height of the colorbar

if mod(nparams,2)>0
	error('ephysPipeline:argChk','Parameters must be specified as parameter/value pairs!');
end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'filedir'
			filedir=varargin{i+1};
		case 'savedir'
			savedir=varargin{i+1};
		case 'fig_title'
			fig_title=varargin{i+1};
		case 'lfp_fs'
			lfp_fs=varargin{i+1};
		case 'phase_bins'
			phase_bins=varargin{i+1};
		case 'filt_name'
			filt_name=varargin{i+1};
		case 'freq_range'
			freq_range=varargin{i+1};
		case 'trial_range'
			trial_range=varargin{i+1};
		case 'min_f'
			min_f=varargin{i+1};
		case 'max_f'
			max_f=varargin{i+1};
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DATA COLLECTION %%%%%%%%%%%%%%%%%%%%

% where to grab the files from?

load(fullfile(filedir,'aggregated_data.mat'),'CHANNELS','EPHYS_DATA'); % get the channel map and LFPs

% frequency range, may need to add Kaiser to condition_signal for tighter filter cutoffs

lfp_data=ephys_denoise_signal(EPHYS_DATA,CHANNELS,LFPCHANNEL);
clear EPHYS_DATA;
lfp_data=ephys_condition_signal(lfp_data,'l','freq_range',freq_range,'medfilt_scale',medfilt_scale,'filt_name',filt_name);
lfp_data=squeeze(lfp_data);

% need to account for subset of trials if used in single unit data 

if ~isempty(trial_range)
	disp(['Truncating trials to ' num2str([trial_range(1) trial_range(end)])]);
	lfp_data=lfp_data(:,trial_range);
	spike_data=spike_data(trial_range);
end

[samples,trials]=size(lfp_data);
lfp_data=lfp_data';
time=[1:samples]./lfp_fs;

% should get a reasonable estimation with Hilbert as long as we've filtered appropriately 
% (zero phase distortions is pretty important here)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PHASE ESTIMATION %%%%%%%%%%%%%%%%%%%


% take the mean phase in a sliding window and subsequently bin
% each time point is simply the middle of the window

% number of hops

indxs=[1:window_size:samples];

% the time of each bin is defined as the CENTER

colidx=1+(0:(length(indxs)-1))*window_size;
time=((colidx-1)+(window_size/2))/lfp_fs;

PHASE_MEAN=zeros(trials,length(indxs)-1);
BINNED_PHASE=zeros(trials,length(indxs)-1);

for i=1:trials

	currphase=angle(hilbert(lfp_data(i,:)));
	currphase=mod(unwrap(currphase),2*pi); % set from 0 to 2pi

	for j=1:length(indxs)-1
		phase_win=currphase(indxs(j):indxs(j+1));
		PHASE_MEAN(i,j)=mean(phase_win);
	end

	[density BINNED_PHASE(i,:)]=histc(PHASE_MEAN(i,:),phase_bins);

end

% image with the histogram
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%


startidx=max([find(HISTOGRAM.f<=min_f)]);

if isempty(startidx)
	startidx=1;
end

stopidx=min([find(HISTOGRAM.f>=max_f)]);

if isempty(stopidx)
	stopidx=length(HISTOGRAM.f);
end


binimage=figure();
ax(1)=subaxis(4,1,1,1,1,2,'margin',.12,'spacingvert',0,'paddingbottom',0);
imagesc(HISTOGRAM.t,HISTOGRAM.f(startidx:stopidx),HISTOGRAM.imask(startidx:stopidx,:));
ylabel('Hz');
set(gca,'XTick',[]);
prettify_axis(ax(1),'FontName','Helvetica','FontSize',15);
prettify_axislabels(ax(1));
box off;
axis xy;
colormap(hist_colors);
freezeColors;

ax(2)=subaxis(4,1,1,3,1,2,'spacingvert',0.1,'margin',0.12,'paddingbottom',0);
imagesc(time,1:trials,BINNED_PHASE);
xlabel('Time (in s)');
ylabel('Trial');
prettify_axis(ax(2),'FontName','Helvetica','FontSize',15);
prettify_axislabels(ax(2));
box off;
colormap(phase_colors);

pos=get(ax(2),'pos');
set(ax(2),'pos',[pos(1) pos(2) pos(3)*.95 pos(4)]);
pos=get(ax(2),'pos');
hc2=colorbar('location','eastoutside','position',[pos(1)+pos(3)+.01 pos(2)+pos(4)/2 colorbarsize pos(4)/2]);
set(hc2,'linewidth',2);
set(hc2,'YTick',[1 2.5 4],'YTickLabel',{'0','p','2p'},'FontName','Symbol','FontSize',15); % these will need to be changed with the phase edges
%ylabel(hc2,'phase','FontSize',15,'FontName','Helvetica');

pos=get(ax(1),'pos');
set(ax(1),'pos',[pos(1) pos(2) pos(3)*.95 pos(4)]);

linkaxes(ax,'x');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


