function [RMASK IMASK q]=contour_approx(SIGNAL,varargin)
%contour_approx computes the contour approximation via Chris' method.
%
%	[REMASK IMASK q]=contour_approx(SIGNAL,varargin)
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

N=2048;
overlap=2030;
tscale=1.5;
fs=25e3;

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
		otherwise
	end
end

t=-N/2+1:N/2;

sigma=(tscale/1e3)*fs;

w = exp(-(t/sigma).^2);
dw = -2*w.*(t/(sigma^2));

q = specgram(SIGNAL,N,[],w,overlap) + eps;
q2 = specgram(SIGNAL,N,[],dw,overlap) + eps;
dx = (q2./q)/(2*pi);

redx = real(dx)./abs(real(dx));
imdx = imag(dx)./abs(imag(dx));

RMASK = redx - circshift(redx,[1 0]) ~= 0 | redx - circshift(redx,[0 1]) ~= 0;
IMASK = imdx - circshift(imdx,[1 0]) ~= 0 | imdx - circshift(imdx,[0 1]) ~= 0;





