function [sono_store,env_xcorr,dur_pdist]=ephys_ratetracker_sono(STORE_AUDIO,STORE_PLAYBACK,varargin)
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

nwin=128;
noverlap=0;
nfft=128;

clipping=-3;

trialwin=100;
trialoverlap=80;

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
%[b2,a2]=ellip(5,.2,60,[audio_range]/(STORE_AUDIO.fs/2),'bandpass');

STORE_AUDIO.data=filtfilt(b,a,STORE_AUDIO.data);
%STORE_PLAYBACK.data=filtfilt(b,a,STORE_PLAYBACK.data);

[nsamples,ntrials]=size(STORE_AUDIO.data);

% take the song data in the notch, use to determine syllable boundaries

[t,f,minpt,maxpt]=getspecgram_dim(nsamples,nwin,noverlap,nfft,STORE_AUDIO.fs,audio_range(1),audio_range(2));

mint=1/STORE_AUDIO.fs;
maxt=nsamples/STORE_AUDIO.fs;
sono_fs=1./(t(2)-t(1));

smooth_win=round(audio_scale.*sono_fs);
sono_store=zeros(length(f),length(t),length(ntrials));


[nblanks formatstring]=fb_progressbar(100);
fprintf(1,['Progress:  ' blanks(nblanks)]);

for i=1:ntrials

	fprintf(1,formatstring,round((i/ntrials)*100));	

	%p0=10*eye(200);
	%ha=adaptfilt.rls(200,.99,p0);
	%[y,e]=filter(ha,STORE_PLAYBACK.data(:,i),STORE_AUDIO.data(:,i));

	%[s]=spectrogram(STORE_AUDIO.data(:,i)-y,nwin,noverlap,nfft,STORE_AUDIO.fs);
	%
	
	%figure();imagesc(log(abs(s)));
	%axis xy;


	[sono_store(:,:,i)]=spectrogram(STORE_AUDIO.data(:,i),nwin,noverlap,nfft,STORE_AUDIO.fs);

	%figure();imagesc(log(abs(sono_store(:,:,i))));
	%axis xy;

	%pause();

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

fprintf(1,'\n');

trialstep=trialwin-trialoverlap;
counter=1;
steps=trialwin:trialstep:ntrials-trialwin;

% use the first 50 trials as a template, check change in zero crossings within song bounds, perhaps also correlate
% the whole envelope


% get the env template
%

ledge=1;
redge=trialwin;

[smooth_energy,pos_zcross]=get_smooth_energy(sono_store(:,:,ledge:redge),smooth_win,minpt,maxpt,t,bound,maxt);

[~,sono_bounds(1)]=min(abs(t-bound));
[~,sono_bounds(2)]=min(abs(t-(maxt-bound)));

env_template=smooth_energy(sono_bounds(1):sono_bounds(2));
dur_template=pos_zcross;

env_xcorr=zeros(1,length(steps));
dur_pdist=zeros(length(dur_template),length(steps));

for i=1:length(steps)

	ledge=steps(i)+1;
	redge=ledge+trialwin;

	[smooth_energy,pos_zcross]=get_smooth_energy(sono_store(:,:,ledge:redge),smooth_win,minpt,maxpt,t,bound,maxt);

	% get all zero crossings, both directions

	im_s=max(clipping,log(abs(mean(sono_store(:,:,ledge:redge),3))));

	env_check=smooth_energy(sono_bounds(1):sono_bounds(2));

	tmp=corrcoef(env_template,env_check);
	env_xcorr(i)=tmp(2,1);
	
	for j=1:length(dur_template)
		[~,tmp]=min(abs(pos_zcross-dur_template(j)))/sono_fs;
		dur_pdist(j,i)=pos_zcross(tmp)/sono_fs;
	end
	
	%figure(1);
	%ax(1)=subplot(2,1,1);imagesc(t,f,im_s);
	%axis xy;
	%hold on;
	%ax(2)=subplot(2,1,2);plot(t,smooth_energy);
	%hold on
	%plot(t(pos_zcross),smooth_energy(pos_zcross),'r*');

	%linkaxes(ax,'x');

	%pause();


end

end

function [SMOOTH_ENERGY,PTS]=get_smooth_energy(SONO,WIN,MINPT,MAXPT,T,BOUND,MAXT)
%
%
%
ave_sono=mean(SONO,3);
ave_abs=log(abs(ave_sono));

song_energy=mean(ave_abs(MINPT:MAXPT,:));

SMOOTH_ENERGY=smooth(song_energy,WIN,'loess');	
SMOOTH_ENERGY=zscore(SMOOTH_ENERGY);
SMOOTH_ENERGY=SMOOTH_ENERGY-smooth(SMOOTH_ENERGY,.3*(1/(T(2)-T(1))));

idx=1:length(SMOOTH_ENERGY)-1;

pos_zcross=find(SMOOTH_ENERGY(idx)<0&SMOOTH_ENERGY(idx+1)>0);
%neg_zcross=find(SMOOTH_ENERGY(idx)>0&SMOOTH_ENERGY(idx+1)<0);

% for each pos_zcross find the turning point 

for i=1:length(pos_zcross)

	% start walking back 
	%
	for j=pos_zcross(i):-1:2

		deriv=SMOOTH_ENERGY(j-1)-SMOOTH_ENERGY(j);

		if deriv>0
			break;

		end

	end

	turning_point=j;
	
	% what's the trough value?
	
	trough_val(i)=SMOOTH_ENERGY(turning_point);

end

pos_zcross_t=T(pos_zcross);
to_del=find(pos_zcross_t<BOUND|pos_zcross_t>(MAXT-BOUND)|trough_val>-.2);

pos_zcross(to_del)=[];
trough_val(to_del)=[];


trough_val

PTS=pos_zcross;

end
% windowed average of the complex sonogram, use as a crude denoiser and recover duration
%
%



