function [STFT DX]=chirp_stft(SIGNAL,varargin)
%
%
%
%

nparams=length(varargin);

if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs');
end


fs=25e3;
n=500;
overlap=450;
nfft=512;
slopec=0; % 0 for Gabor
wsigma=.15;

for i=1:2:nparams
	switch lower(varargin{i})
		case 'fs'
			fs=varargin{i+1};
		case 'overlap'
			overlap=varargin{i+1};
		case 'nfft'
			nfft=varargin{i+1};
		case 'slopec'
			slopec=varargin{i+1};
		case 'n'
			n=varargin{i+1};
		case 'wsigma'
			wsigma=varargin{i+1};
	end
end

t=-n/2:n/2-1; 
wsigma=wsigma*fs;

% gauss and Dgauss

%w=exp(-(t/sigma).^2);
%dw=(w).*((t)/(sigma^2))*-2;

% chirplet if slopec>0

alpha=1/wsigma^2-1j*slopec;
z=power(2/wsigma^2,.25);

wv=spectrogram(SIGNAL,z*exp(-pi*alpha*(t.*t)),overlap,nfft);
dwv=spectrogram(SIGNAL,z*(-2*pi*t).*exp(-pi*alpha*(t.*t)),overlap,nfft);

DX=-dwv./wv;
STFT=wv;

    

