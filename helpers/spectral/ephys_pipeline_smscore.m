function [features]=ephys_pipeline_smscore(s,fs,varargin)
%computes spectral features of a given signal
%
%


nparams=length(varargin);

if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs');
end

n=1024; % spectrogram window size
overlap=1000; % spectrogram overlap
sigma=1.5; % Gaussian timescale (in ms)
downsampling=5; % downsampling factor (skip columns)
filter_scale=10; % disk filter scale (samples)

for i=1:2:nparams
	switch lower(varargin{i})
		case 'n'
			n=varargin{i+1};
		case 'overlap'
			overlap=varargin{i+1};
		case 'sigma'
			sigma=varargin{i+1};
		case 'filter_scale'
			filter_scale=varargin{i+1};
		case 'downsampling'
			downsampling=varargin{i+1};
	end
end

% TODO remove dynamic allocation of feature matrix

disp('Computing score');
t=-n/2+1:n/2;

sigma=(sigma/1000)*fs;

%Gaussian and first derivative as windows.
% let's remove redundant angles and gradients, maybe just cos

w=exp(-(t/sigma).^2);
dw=(w).*((t)/(sigma^2))*-2;
q=specgram(s,n,[],w,overlap)+eps; %gaussian windowed spectrogram
q2=specgram(s,n,[],dw,overlap)+eps; %deriv gaussian windowed spectrogram

% add FM and pitch?

dx=(q2./q)/(2*pi); %displacement according to the remapping algorithm

sonogram=flipdim(q,1);
dx=flipdim(dx,1);

% take subset of frequencies to focus on

dx=dx(300:480,:);
sonogram=sonogram(300:480,:);

% larger disks really slows down compute time

H = fspecial('disk',filter_scale);

%compute local angles

s1=abs(cos(angle(dx)));

%filter the angle images.

blurred = imfilter(s1,H,'circular');
blurredcm = imfilter(log(abs(sonogram)),H,'circular');

[fx,fy]=gradient(abs(cos(angle(dx))));

sfx=imfilter(abs(fx),H,'circular');
sfy=imfilter(abs(fy),H,'circular');

v{1}=blurred;
v{2}=sfx;
v{3}=sfy;
v{4}=blurredcm;

[a,b]=size(v{1});

for k=1:length(v)
	jj=1;
	
	%subsample the image by grouping columns

	features{k}=[];

	for i=1:downsampling:b-downsampling
		features{k}(:,jj)=sum(v{k}(:,i:i+downsampling),2);
		jj=jj+1;
	end
end

