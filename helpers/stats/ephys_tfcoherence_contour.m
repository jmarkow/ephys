function [COH,T,F,SPECT1,SPECT2,STATS,FORCINGFUNCTION]=ephys_tfcoherence_contour(SIGNAL1,SIGNAL2,varargin)
%Contour-based time-frequency coherence
%
%
%

STATS=[];

nparams=length(varargin);

nfft=10000;
n=6250;
overlap=6000;
w=2;

ntapers=[];
alpha=[.0001 .001 .01 .05];

angles=-pi/4:pi/16:pi/4; % for contour image
ts=50:20:150; %timescales in milliseconds;
contour_thresh=98;

fs=25e3;

for i=1:2:nparams
	switch lower(varargin{i})
		case 'nfft'
			nfft=varargin{i+1};
		case 'n'
			n=varargin{i+1};
		case 'overlap'
			overlap=varargin{i+1};
		case 'ntapers'
			ntapers=varargin{i+1};
		case 'fs'
			fs=varargin{i+1};
		case 'alpha'
			alpha=varargin{i+1};
		case 'w'
			w=varargin{i+1};
		case 'angles'
			angles=varargin{i+1};
		case 'ts'
			ts=varargin{i+1};
		case 'contour_thresh'
			contour_thresh=varargin{i+1};
	end
end

[trials1,samples1]=size(SIGNAL1);
[trials2,samples2]=size(SIGNAL2);
ntimescales=length(ts);
nangles=length(angles);

if trials1~=trials2
	error('Need same number of trials for both signals!');
end

if samples1~=samples2
	error('Need same number of samples for both signals!');
end

nsamples=samples1;
ntrials=trials1;

% set up spectrogram parameters

if n>nsamples
	difference=n-overlap;
	n=round(nsamples/5);
	overlap=n-difference;
	disp('Reset window size and overlap, longer than nsamples');
	disp(['Window size:  ' num2str(n) ' samples']);
	disp(['Overlap:  ' num2str(overlap) ' samples']);
end

if isempty(ntapers)
    ntapers=2*(w)-1;
end

if isempty(nfft)
	nfft=max([n 2^nextpow2(n)]);
else
	nfft=2^nextpow2(nfft);
end

if nfft~=n
	disp('Setting nfft to n');
	nfft=n;
end

if overlap~=n-1
	disp('Setting overlap to 1 sample');
	overlap=n-1;
end

[T,F]=getspecgram_dim(nsamples,n,overlap,nfft,fs);

rows=length(F);
columns=length(T);

% compute any statistics necessary for interpreting the results

% STATS will need to be computed non-parametrically
% pre-allocate our reduction matrices

cross_spect_mean=zeros(rows,columns);
spect1_mean=zeros(rows,columns);
spect2_mean=zeros(rows,columns);

parfor i=1:ntrials

	disp(['Trial ' num2str(i)]);

	% should be zero mean in the first place

	curr1=SIGNAL1(i,:)-mean(SIGNAL1(i,:));
	curr2=SIGNAL2(i,:)-mean(SIGNAL2(i,:));

	% spectra

	% multi-taper estimate
	% accumulate across tapers to cut down compt time, then we can use parfor

	cross_spect_tmp=zeros(rows,columns);
	spect1_tmp=zeros(rows,columns);
	spect2_tmp=zeros(rows,columns);
	consensus_tmp=zeros(rows,columns);

	% use the first signal for consensus weighting

	for j=1:ntimescales

		[dxsub1 spect1]=chirp2(0,ts(j),curr1,fs,overlap,nfft);
		[dxsub2 spect2]=chirp2(0,ts(j),curr2,fs,overlap,nfft);

		% cross spectrum
		
		cross=spect1.*conj(spect2);

		% take power of two signals

		power1=spect1.*conj(spect1);
		power2=spect2.*conj(spect2);

		% average across tapers

		for jj=1:nangles

			a_consensus=chirp_consensus(curr1,fs,contour_thresh,...
				angles(jj),nfft,dxsub1,0,jj,1);	
			consensus_tmp=consensus_tmp+(cross.*a_consensus)./length(angles);

		end

		cross_spect_tmp=cross_spect_tmp+consensus_tmp./(ntapers);
		spect1_tmp=spect1_tmp+power1./(ntapers);
		spect2_tmp=spect2_tmp+power2./(ntapers);
	end

	% average across trials

	cross_spect_mean=cross_spect_mean+cross_spect_tmp./ntrials;
	spect1_mean=spect1_mean+spect1_tmp./ntrials;
	spect2_mean=spect2_mean+spect2_tmp./ntrials;

end

SPECT1=spect1_mean;
SPECT2=spect2_mean;

COH=cross_spect_mean./sqrt(spect1_mean.*spect2_mean);

disp('Computing forcing function');

forcing=zeros(rows,columns);
parfor i=1:ntrials

	curr1=SIGNAL1(i,:)-mean(SIGNAL1(i,:));

	for j=1:ntimescales
		[dxsub1 spect]=chirp2(0,ts(j),curr1,fs,overlap,nfft);
		forcing=forcing+(spect.*abs(COH))./(ntrials*ntimescales);
	end
end

absv=abs(forcing);
preang=angle(forcing);

% phase shift every other line by pi

for jj=1:2:size(preang,1)
	preang(jj,:) = preang(jj,:) + pi;
end

preang = mod(preang,2*pi);
sono = absv.*(cos(preang)+1i*sin(preang));

FORCINGFUNCTION=sum(real(sono));

