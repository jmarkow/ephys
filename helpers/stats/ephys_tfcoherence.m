function [COH,T,F,SPECT1,SPECT2,STATS]=ephys_tfcoherence(SIGNAL1,SIGNAL2,varargin)
%Slepian-based time-frequency coherence
%
%
%

nparams=length(varargin);

nfft=10000;
n=6250;
overlap=6000;
w=2;

ntapers=[];
alpha=[.0001 .001 .01 .05];

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
	end
end

[trials1,samples1]=size(SIGNAL1);
[trials2,samples2]=size(SIGNAL2);

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

STATS.dof=2*ntapers*ntrials; % correction added 10/2/12 (forgot 2!)
[tapers,lambda]=dpss(n,w,ntapers);

if isempty(nfft)
	nfft=max([n 2^nextpow2(n)]);
else
	nfft=2^nextpow2(nfft);
end

resolution=w*1/(n/fs);
disp(['Resolution:  ' num2str(resolution)  ' Hz']);
disp(['NFFT:  ' num2str(nfft)]);

if n>nsamples
	difference=n-overlap;
	n=round(nsamples/5);
	overlap=n-difference;
	disp('Reset window size and overlap, longer than nsamples');
	disp(['Window size:  ' num2str(n) ' samples']);
	disp(['Overlap:  ' num2str(overlap) ' samples']);
end

[T,F]=getspecgram_dim(nsamples,n,overlap,nfft,fs);

rows=length(F);
columns=length(T);

% compute any statistics necessary for interpreting the results

STATS.ncomp=length(T)*length(F);

STATS.alpha=alpha;

% compute the null distribution cutoffs for particular alpha values

STATS.null=sqrt(1-(alpha).^(1/(STATS.dof/2-1)));

% bonferonni corrected values

STATS.null_boncorrected=sqrt(1-(alpha./STATS.ncomp).^2).^(1/(STATS.dof/2-1));

% compute the variance according to the asymptotic measure

STATS.var=1.96/(sqrt(STATS.dof/2));
STATS.tapers=tapers;
STATS.lambda=lambda;
STATS.ntapers=ntapers;

% pre-allocate our reduction matrices

cross_spect_mean=zeros(rows,columns);
spect1_mean=zeros(rows,columns);
spect2_mean=zeros(rows,columns);

parfor i=1:ntrials

	disp(['Trial ' num2str(i)]);

	% should be zero mean in the first place

	curr1=SIGNAL1(i,:)-mean(SIGNAL1(i,:));
	curr2=SIGNAL2(i,:)-mean(SIGNAL2(i,:));


	% multi-taper estimate
	% accumulate across tapers to cut down compt time, then we can use parfor

	cross_spect_tmp=zeros(rows,columns);
	spect1_tmp=zeros(rows,columns);
	spect2_tmp=zeros(rows,columns);

	for j=1:ntapers

		spect1=spectrogram(curr1,tapers(:,j)',overlap,nfft);
		spect2=spectrogram(curr2,tapers(:,j)',overlap,nfft);
		
		% cross spectrum
		
		cross=spect1.*conj(spect2);

		% take power of two signals

		power1=spect1.*conj(spect1);
		power2=spect2.*conj(spect2);

		% average across tapers

		cross_spect_tmp=cross_spect_tmp+cross./(ntapers);
		spect1_tmp=spect1_tmp+power1./(ntapers);
		spect2_tmp=spect2_tmp+power2./(ntapers);
	end

	% average across trials

	cross_spect_mean=cross_spect_mean+cross_spect_tmp./ntrials;
	spect1_mean=spect1_mean+spect1_tmp./ntrials;
	spect2_mean=spect2_mean+spect2_tmp./ntrials;

end

% phase correction

SPECT1=spect1_mean;
SPECT2=spect2_mean;

% z-transform as well

COH=cross_spect_mean./sqrt(spect1_mean.*spect2_mean);

% display is user selects to display

% return p-values

abscoh=abs(COH);
STATS.pcoh=(STATS.dof-2).*abscoh.*(1-abscoh.^2).^(STATS.dof/2-2);

