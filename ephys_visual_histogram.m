function HISTOGRAM=ephys_visual_histogram(MIC_DATA,varargin)
%computes a contour histogram for a group of sounds
%
%	HISTOGRAM=ephys_visual_histogram(MIC_DATA,varargin)
%
%	MIC_DATA
%	samples x trials matrix of aligned sounds
%	
%	the following may be specified as parameter/value pairs:%
%		
%		fs
%		sampling frequency (default: 25e3)
%
%		tscale
%		time scale for Gaussian window for the Gabor transform (in ms, default: 1.5)
%	
%		n
%		window length
%
%		nfft
%		number of points in fft
%
%		overlap
%		window overlap
%
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETER COLLECTION %%%%%%%%%%%%%%%%%

if nargin<1
	error('ephysPipeline:tfhistogram:notenoughparams','Need 1 argument to continue, see documentation');
end

nparams=length(varargin);

if mod(nparams,2)>0
	error('ephysPipeline:argChk','Parameters must be specified as parameter/value pairs!');
end

fs=25e3;
tscale=1.5;
N=1024;
nfft=1024;
overlap=1e3;
mic_filtering=500; % highpass for mic trace
mask_only=0;
spect_thresh=.2;

for i=1:2:nparams
	switch lower(varargin{i})
		case 'fs'
			fs=varargin{i+1};
		case 'tscale'
			tscale=varargin{i+1};
		case 'n'
			N=varargin{i+1};
		case 'nfft'
			nfft=varargin{i+1};
		case 'overlap'
			overlap=varargin{i+1};
		case 'mic_filtering'
			filtering=varargin{i+1};
		case 'mask_only'
		    mask_only=varargin{i+1};
		case 'spect_thresh'
		    spect_thresh=varargin{i+1};

	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% compute the contour histogram
% normalize the mic trace

% 

[nsamples,ntrials]=size(MIC_DATA);

if ~isempty(mic_filtering)
	[b,a]=butter(2,[mic_filtering/(fs/2)],'high');
	MIC_DATA=filtfilt(b,a,MIC_DATA);
end

MIC_DATA=MIC_DATA./repmat(max(abs(MIC_DATA),[],1),[nsamples 1]);

[HISTOGRAM.rmask HISTOGRAM.imask HISTOGRAM.f HISTOGRAM.t]=contour_histogram(MIC_DATA,'fs',fs,...
	'tscale',tscale,'nfft',nfft,'n',N,'overlap',overlap,'mask_only',mask_only,'spect_thresh',spect_thresh);

% mean oscillogram

HISTOGRAM.mean_osc=mean(MIC_DATA,2);
