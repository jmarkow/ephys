function [binned_phase,phase_mean]=ephys_su_lfp_phasebin(LFPCHANNEL,SUCHANNEL,SUCLUSTER,HISTOGRAM,varargin)
%computes coherograms between LFPs and single units
%
%	[abscoh,t,f]=ephys_su_lfp_coherence_tf(LFPCHANNEL,SUCHANNEL,SUCLUSTER,HISTOGRAM,varargin)
%	
%	LFPCHANNEL
%	LFPCHANNEL to use
%
%	SUCHANNEL
%	channel with single units (already clustered and saved in 'sua')
%
%	SUCLUSTER
%	number of cluster with single unit
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
%		colormap for coherence graph (default: 'jet')
%		
%		n
%		spectrogram window length (default: 6.25e3 samples)
%		
%		nfft
%		spectrogram nfft (default: 10e3 samples)
%
%		overlap
%		spectrogram overlap (default: 6e3 samples)
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
%		nphasebins
%		number of phase bins for amplitude-weighted frequency/phase 2D histogram (default: 50)
%
% see also ephys_visual_histogra.m,ephys_su_lfp_coherence_spect.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETER COLLECTION %%%%%%%%%%%%%%%%%


nparams=length(varargin);

savedir=pwd;
filedir=pwd;
fig_title='noname';
phasecolors='hsv'; % let's do a custom colormap with colors matched to nphasebins
histcolors='jet'; % colors for histogram
window_size=25; % in samples how large a window size do we want for the non-overlapping average?

ntapers=[];
freq_range=[300]; % let's just refilter
lfp_fs=25e3; % default Intan sampling rate
trial_range=[];
medfilt_scale=1.5; % median filter scale (in ms)
phase_bins=[];

if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs');
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
		case 'freq_range'
			freq_range=varargin{i+1};
		case 'trial_range'
			trial_range=varargin{i+1};
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DATA COLLECTION %%%%%%%%%%%%%%%%%%%%

% where to grab the files from?

savename=[ fig_title '_lfpch_' num2str(LFPCHANNEL) '_such' num2str(SUCHANNEL) '_clust' num2str(SUCLUSTER)];
sua_mat=fullfile(filedir,'sua',['sua_channels ' num2str(SUCHANNEL) '.mat']);

load(fullfile(filedir,'aggregated_data.mat'),'CHANNELS','EPHYS_DATA'); % get the channel map and LFPs
load(sua_mat,'smooth_spikes','clust_spike_vec','subtrials'); % smooth spikes

% frequency range, may need to add Kaiser to condition_signal for tighter filter cutoffs

lfp_data=ephys_denoise_signal(EPHYS_DATA,CHANNELS,LFPCHANNEL);
clear EPHYS_DATA;
lfp_data=ephys_condition_signal(lfp_data,'l','freq_range',freq_range,'medfilt_scale',medfilt_scale);
lfp_data=squeeze(lfp_data);

% need to account for subset of trials if used in single unit data 

lfp_data=lfp_data(:,subtrials);
spike_data=clust_spike_vec{1}{SUCLUSTER}; % spike times (in ms)

if ~isempty(trial_range)
	disp(['Truncating trials to ' num2str([trial_range(1) trial_range(end)])]);
	lfp_data=lfp_data(:,trial_range);
	spike_data=spike_data(trial_range);
end

[samples,trials]=size(lfp_data);
lfp_data=lfp_data';

% should get a reasonable estimation with Hilbert as long as we've filtered appropriately 
% (zero phase distortions is pretty important here)

% set up the phase bins to count with histc

%phase_bins=linspace(-pi,pi,nphasebins)
phase_bins=[-inf pi/4 2*pi/4 3*pi/4 inf];

% take the mean phase in a sliding window and subsequently bin

for i=1:trials
	indxs=[1:window_size:samples];
	currphase=angle(hilbert(lfp_data(i,:)));
	%currphase=unwrap(currphase);

	for j=1:length(indxs)-1
		phase_win=currphase(indxs(j):indxs(j+1));
		phase_mean(i,j)=mean(phase_win);
	end

	[density binned_phase(i,:)]=histc(phase_mean(i,:),phase_bins);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


