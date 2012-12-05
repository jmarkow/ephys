function [dxsub sonosub]=chirp2(slopec,sigma,signal,SAMPLING,OVERLAP,NFFT);

%seq_length = NFFT; 
%time_halfbandwidth = 3.5;
%num_seq = 2*(2.5)-1;

%Obtain DPSSs

%[dps_seq,lambda] = dpss(seq_length,time_halfbandwidth,num_seq);

N=NFFT;

t=-N/2:N/2-1; % Corrected to shift a window to be on the center (sep 25.2009)

%Gaussian and first derivative as windows.

sigma=(sigma/1000)*SAMPLING;
w=exp(-(t/sigma).^2);
dw=(w).*((t)/(sigma^2))*-2;

% chirplet if slopec>0

alpha=1/sigma^2-1j*slopec;
z=power(2/sigma^2,.25);

wv=spectrogram(signal,z*exp(-pi*alpha*(t.*t)),OVERLAP);
dwv=spectrogram(signal,z*(-2*pi*t).*exp(-pi*alpha*(t.*t)),OVERLAP);

dxsub=-dwv./wv;
sonosub=wv;

    

