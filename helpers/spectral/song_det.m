function [song_idx power f t song_detvec]=song_det(audio,fs,minfs,maxfs,window,noverlap,songduration,ratio_thresh,song_thresh,pow_thresh)
%based on Andalmann's code

if nargin<10 | isempty(pow_thresh)
    pow_thresh=0;
end

[s,f,t]=spectrogram(audio,window,noverlap,[],fs);

% take the power and find our fs band

power=abs(s);
min_idx=max(find(f<=minfs));
max_idx=min(find(f>=maxfs));

% take the song/nonsong power ratio

song=mean(power(min_idx:max_idx,:),1);
nonsong=mean(power([1:min_idx-1 max_idx+1:end],:),1)+eps;

song_ratio=song./nonsong;
%song_detvec=smooth(double(song_ratio>ratio_thresh),round((fs*songduration)/(window-noverlap)));

% convolve with a moving average filter

filt_size=round((fs*songduration)/(window-noverlap));
mov_filt=ones(1,filt_size)*1/filt_size;
song_detvec=conv(double(song_ratio>ratio_thresh),mov_filt,'same');

% where is the threshold exceeded for both the raw power and the ratio?

pow_idx=song>pow_thresh;
ratio_idx=song_detvec>song_thresh;

%%%%

song_idx=pow_idx&ratio_idx;

