function [RMASK IMASK F T CONTOURS]=contour_histogram(SIGNALS,varargin)
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
%		SR
%		sampling rate (default: 48e3)
%

N=1025;
overlap=1000;
tscale=1.5; % in msecs
fs=48e3;
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
		case 'tscale'
			tscale=varargin{i+1};
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

[rmask_pre imask_pre spect]=contour_approx(SIGNALS(:,1),'fs',fs,'N',N,'overlap',overlap,'tscale',tscale);

if mask_only
	RMASK=rmask_pre./trials;
	IMASK=imask_pre./trials;
else
	RMASK=((rmask_pre.*abs(spect))>spect_thresh)./trials;
	IMASK=((imask_pre.*abs(spect))>spect_thresh)./trials;
end

[rows,columns]=size(rmask_pre);

re_contours=zeros(rows,columns,trials,'uint8');
im_contours=zeros(rows,columns,trials,'uint8');

re_contours(:,:,1)=uint8(rmask_pre);
im_contours(:,:,1)=uint8(imask_pre);

parfor i=1:trials
	
	[rmask_pre imask_pre spect]=contour_approx(SIGNALS(:,i),'fs',fs,'N',N,'overlap',overlap,'tscale',tscale);

	re_contours(:,:,i)=uint8(rmask_pre);
	im_contours(:,:,i)=uint8(imask_pre);

	if mask_only
		RMASK=RMASK+rmask_pre./trials;
		IMASK=IMASK+imask_pre./trials;
	else
		RMASK=RMASK+(((rmask_pre.*abs(spect))>spect_thresh))./trials;
		IMASK=IMASK+(((imask_pre.*abs(spect))>spect_thresh))./trials;
	end

end

[rows,cols]=size(RMASK);

CONTOURS.re=re_contours;
CONTOURS.im=im_contours;

% shamelessly cribbed from MATLAB's computation for spectrogram

% should scale 1:fbins * nyquist

F=((1:rows)./rows).*(fs/2);

% starting at 1 one hop is n-overlap samples

col_idx=1+(0:(cols-1))*(N-overlap);

% then in time each step is samples + window/2 offset /SR

T=((col_idx-1)+((N/2)'))/fs;



