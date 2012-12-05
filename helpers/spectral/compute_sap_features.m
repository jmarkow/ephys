function [T FEATURES SPECDERIV]=compute_sap_features(SIGNALS,varargin)
%contour_approx computes the contour approximation via Chris' method.
%
%	[REMASK IMASK F T CONTOURS]=contour_histogram(SIGNAL,varargin)
%
%	SIGNAL
%	vector that contains the signal of interest
%
%	the following additional parameters may be specified as parameter/value pairs
%
%		N
%		length of the window for spectrogram (default: 2048)
%
%		overlap
%		window overlap (overlap<N) (default: 2030)
%
%		tscale
%		timescale of filtering in ms (default: 1.5ms)
%
%		fs
%		sampling rate (default: 48e3)
%

N=409;
overlap=365;
fs=48e3;
nfft=1024;
spect_thresh=.2;
mask_only=0;

nparams=length(varargin);

if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs');
end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'n'
			N=varargin{i+1};
		case 'overlap'
			overlap=varargin{i+1};
		case 'fs'
			fs=varargin{i+1};
		case 'spect_thresh'
			spect_thresh=varargin{i+1};
		case  'mask_only'
			mask_only=varargin{i+1};
		otherwise
	end
end

[samples,trials]=size(SIGNALS);
[t,f]=getspecgram_dim(samples,N,overlap,nfft,fs);
T=t;

features={'amplitude','am','fm','entropy','gravityc','pitchgood','pitch','pitchweight','pitchchose'};

for i=1:length(features)
	FEATURES.(features{i})=zeros(length(t),trials);
end


spectmp=sap_features(SIGNALS(:,1),fs,'N',N,'overlap',overlap,'nfft',nfft);

[m,n]=size(spectmp);
SPECDERIV=zeros(m,n,trials);

for i=1:trials
	[SPECDERIV(:,:,i) FEATURES.am(:,i) FEATURES.fm(:,i) FEATURES.entropy(:,i) FEATURES.amplitude(:,i) ...
	       FEATURES.gravityc(:,i) FEATURES.pitchgood(:,i) FEATURES.pitch(:,i) FEATURES.pitchchose(:,i) FEATURES.pitchweight(:,i)] ...
		=sap_features(SIGNALS(:,i),fs,'N',N,'overlap',overlap,'nfft',nfft);

end

