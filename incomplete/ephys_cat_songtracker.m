function [sono_store]=ephys_ratetracker_sono(STORE_AUDIO,STORE_PLAYBACK,varargin)
%generates song-aligned single-unit rasters
%
%	ephys_visual_sua(EPHYS.data,HISTOGRAM,EPHYS.labels,varargin)
%
%	EPHYS
%	structure with the following fields
%
%	EPHYS.data
%	sound-aligned voltage traces from extracted_data.mat (should be the variable ephys_data)
%
%	EPHYS.labels
%	channel labels (i.e. the channel that corresponds to a given element in the cell array
%	ephys_data) from extracted_data.mat
%
%	HISTOGRAM
%	contour histogram returned by ephys_visual_histogram.m (or loaded from histogram.mat)
%
%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETER COLLECTION %%%%%%%%%%%%%%%%%

nparams=length(varargin);

if mod(nparams,2)>0
	error('ephysPipeline:argChk','Parameters must be specified as parameter/value pairs!');
end

subtrials=[];

audio_range=[2e3 4.5e3];
audio_thresh=.1;
audio_scale=.08;

savedir='';
savefile='cat_song.mat';

bound=.5;

% need sufficient resolution 

nwin=256;
noverlap=200;
nfft=256;

clipping=-3;

trialwin=100;
trialoverlap=20;

% remove eps generation, too slow here...

for i=1:2:nparams
	switch lower(varargin{i})
		case 'savedir'
			savedir=varargin{i+1};
		case 'bound'
			bound=varargin{i+1};
		case 'savedir'
			savedir=varargin{i+1};
		case 'savefile'
			savefile=varargin{i+1};
		case 'nwin'
			nwin=varargin{i+1};
		case 'nfft'
			nfft=varargin{i+1};
		case 'noverlap'
			noverlap=varargin{i+1};
		case 'clipping'
			clipping=varargin{i+1};
		case 'trialwin'
			trialstep=varargin{i+1};
		case 'trialoverlap'
			trialoverlap=varargin{i+1};
	end
end

[b,a]=ellip(5,.2,60,[800]/(STORE_AUDIO.fs/2),'high');

STORE_AUDIO.data=filtfilt(b,a,STORE_AUDIO.data);
STORE_PLAYBACK.data=filtfilt(b,a,STORE_PLAYBACK.data);

[nsamples,ntrials]=size(STORE_AUDIO.data);

% take the song data in the notch, use to determine syllable boundaries

[t,f,minpt,maxpt]=getspecgram_dim(nsamples,nwin,noverlap,nfft,STORE_AUDIO.fs,audio_range(1),audio_range(2));

mint=1/STORE_AUDIO.fs;
maxt=nsamples/STORE_AUDIO.fs;

smooth_win=audio_scale.*(1/(t(2)-t(1)));
sono_store=zeros(length(f),length(t),length(ntrials));

for i=1:ntrials

	%ha=adaptfilt.lms(500,.6,.4);
	%[y,e]=filter(ha,STORE_PLAYBACK.data(:,i),STORE_AUDIO.data(:,i));

	%[s]=spectrogram(STORE_AUDIO.data(:,i)-y,nwin,noverlap,nfft,STORE_AUDIO.fs);
	%
	
	[sono_store(:,:,i)]=spectrogram(STORE_AUDIO.data(:,i),nwin,noverlap,nfft,STORE_AUDIO.fs);

	% take energy in the relevant band

	%song_energy=max(clipping,log(mean(abs(s(minpt:maxpt,:)))));

	% smooth, moving average
	%smooth_energy=song_energy;
	%%smooth_energy=mfcc_deltacoef(smooth_energy,3);
	%%
	%smooth_energy=smooth(song_energy,smooth_win,'loess');
	%
	%%smooth_t=[1:length(smooth_energy)];

	%idx=1:length(smooth_energy)-1;

	%% get rising and falling edges
	%
	%rising_edges=(smooth_energy(idx)<audio_thresh)&(smooth_energy(idx+1)>audio_thresh);
	%falling_edges=(smooth_energy(idx)>audio_thresh)&(smooth_energy(idx+1)<audio_thresh);

	%im_s=abs(s);
	%im_s=im_s./max(im_s(:));
	%im_s=max(clipping,log(im_s));

	%figure(1);
	%ax(1)=subplot(2,1,1);imagesc(t,f,im_s);
	%axis xy;
	%hold on;
	%ax(2)=subplot(2,1,2);plot(t,smooth_energy);

	%t(rising_edges)
	%t(falling_edges)

	%linkaxes(ax,'x');

	%pause();
end

trialstep=trialwin-trialoverlap;
counter=1;
steps=1:trialstep:ntrials-trialwin;

for i=1:length(steps)

	ledge=steps(i);
	redge=ledge+trialwin-1;

	ave_sono=mean(sono_store(:,:,ledge:redge),3);
	%ave_abs=max(clipping,log(abs(ave_sono)));
	ave_abs=log(abs(ave_sono));

	song_energy=mean(ave_abs(minpt:maxpt,:));
	smooth_energy=smooth(song_energy,smooth_win,'loess');
	%smooth_energy=zscore(smooth_energy);
	%smooth_energy=smooth_energy-smooth(smooth_energy,.3*(1/(t(2)-t(1))));
	im_s=max(clipping,log(abs(ave_sono)));

	figure(1);
	ax(1)=subplot(2,1,1);imagesc(t,f,im_s);
	axis xy;
	hold on;
	ax(2)=subplot(2,1,2);plot(t,smooth_energy);

	linkaxes(ax,'x');

	pause();

end

% windowed average of the complex sonogram, use as a crude denoiser and recover duration
%
%



