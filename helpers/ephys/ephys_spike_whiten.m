function NEWSPIKES=ephys_spike_whiten(SPIKEWINDOWS,NOISE_DATA,varargin)
%
%
%
%
%



nparams=length(varargin);

if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs!');
end

maxnoisetraces=1e6;
regularize=.01;

for i=1:2:nparams
	switch lower(varargin{i})
		case 'maxnoisetraces'	
			maxnoisetraces=varargin{end+1};
		case 'regularize'
			regularize=varargin{end+1};
	end
end

nsamples=length(NOISE_DATA);
noisetrials=floor(length(NOISE_DATA)/nsamples);

if noisetrials>maxnoisetraces
	noisetrials=maxnoisetraces;
end

disp(['Noise trials ' num2str(noisetrials)]);

noisematrix=zeros(noisetrials,nsamples);

counter=0;
for i=1:noisetrials
	noisematrix(i,:)=NOISE_DATA(counter+1:counter+nsamples);
	counter=counter+nsamples;
end

noisematrix=noisematrix+regularize.*randn(size(noisematrix));
noisecov=cov(noisematrix);

R=chol(noisecov);
invR=inv(R);

NEWSPIKES=[SPIKEWINDOWS'*invR]';
