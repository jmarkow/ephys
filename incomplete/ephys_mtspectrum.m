function [SPECT,F,TRIALNUM,TAPERNUM,TAPERS]=ephys_mtspectrum(SIGNAL,varargin)
%Slepian-based time-frequency coherence
%
%
%

nparams=length(varargin);

nfft=1024;
w=2;

ntapers=[];
fs=25e3;

for i=1:2:nparams
	switch lower(varargin{i})
		case 'nfft'
			nfft=varargin{i+1};
		case 'overlap'
			overlap=varargin{i+1};
		case 'ntapers'
			ntapers=varargin{i+1};
		case 'fs'
			fs=varargin{i+1};
		case 'w'
			w=varargin{i+1};
	end
end

[ntrials,nsamples]=size(SIGNAL);
n=nsamples;

if isempty(ntapers)
    ntapers=2*(w)-1;
end

[tapers,lambda]=dpss(n,w,ntapers);

if isempty(nfft)
	nfft=max([n 2^nextpow2(n)]);
else
	nfft=2^nextpow2(nfft);
end

resolution=w*1/(n/fs);
disp(['Resolution:  ' num2str(resolution)  ' Hz']);
disp(['NFFT:  ' num2str(nfft)]);

F=linspace(0,1,nfft/2+1).*fs/2;
columns=length(F);
rows=ntrials*ntapers;

SPECT=zeros(rows,columns);
TRIALNUM=zeros(rows,1);
TAPERNUM=zeros(size(TRIALNUM));

TAPERS=tapers(:,1:ntapers);

counter=1;
for i=1:ntrials

	disp(['Trial ' num2str(i)]);

	currtrial=detrend(SIGNAL(i,:)-mean(SIGNAL(i,:)));
		
	% multi-taper estimate
	% accumulate across tapers to cut down compt time, then we can use parfor

	for j=1:ntapers

		spect=fft(currtrial.*tapers(:,j)',nfft)./n; % two sided spectrum
		power=spect.*conj(spect); % get the energy spectrum
		SPECT(counter,:)=2.*power(1:nfft/2+1); % take one sided (multiply by 2 to conserve energy)
		counter=counter+1;

	end
end

% phase correction

