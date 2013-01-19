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
%		savedir
%		if defined, saves the results in savedir in histogram.mat (leave blank to skip, default: pwd)
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
savedir=pwd;
N=1024;
nfft=1024;
overlap=1e3;


for i=1:2:nparams
	switch lower(varargin{i})
		case 'fs'
			fs=varargin{i+1};
		case 'tscale'
			tscale=varargin{i+1};
		case 'savedir'
			savedir=varargin{i+1};
		case 'n'
			N=varargin{i+1};
		case 'nfft'
			nfft=varargin{i+1};
		case 'overlap'
			overlap=varargin{i+1};

	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% compute the contour histogram

[HISTOGRAM.rmask HISTOGRAM.imask HISTOGRAM.f HISTOGRAM.t]=contour_histogram(MIC_DATA,'fs',fs,...
	'tscale',tscale,'nfft',nfft,'n',N,'overlap',overlap);

% mean oscillogram

HISTOGRAM.mean_osc=mean(MIC_DATA,2);

% save results if specified by the user

if ~isempty(savedir)
	save(fullfile(savedir,'histogram.mat'),'HISTOGRAM');
end

