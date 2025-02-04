function [coh,spect1,spect2,freqs,store_cross,store_spect1,store_spect2,stats]=ephys_su_lfp_coherence_spect(SIGNAL1,SIGNAL2,varargin)
%computes coherency spectra between fields and single units
%
%
%
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETER COLLECTION %%%%%%%%%%%%%%%%%

nparams=length(varargin);

nfft=[];
w=2;

% about 200 trials seems to be our limit here...

ntapers=[];
beta=1.5; % parameter for z-transforming coherence
freq_range=[300]; % let's just refilter
fs=1e3; % default Intan sampling rate
trial_range=[];
alpha=[.001 .01 .05]; % alpha for null hypothesis line 
n=500;
nfft=1024;
overlap=450;

if mod(nparams,2)>0
	error('ephysPipeline:argChk','Parameters must be specified as parameter/value pairs!');
end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'debug'
			debug=varargin{i+1};
		case 'nfft'
			nfft=varargin{i+1};
		case 'alpha'
			alpha=varargin{i+1};
		case 'fs'
			fs=varargin{i+1};
	        case 'w'
            		w=varargin{i+1};
		case 'ntapers'
			ntapers=varargin{i+1};
		case 'n'
			n=varargin{i+1};
		case 'nfft'
			nfft=varargin{i+1};
		case 'overlap'
			overlap=varargin{i+1};
	end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DATA COLLECTION %%%%%%%%%%%%%%%%%%%%

[trials,samples]=size(SIGNAL1);
[trials,samples]=size(SIGNAL2);

if isempty(nfft)
	nfft=max([samples 2^nextpow2(samples)]);
end

if isempty(ntapers)
    ntapers=2*(w)-1;
end

resolution=w*1/(samples/fs);
disp(['Resolution:  ' num2str(resolution)  ' Hz']);
disp(['NFFT:  ' num2str(nfft)]);

% take the average spectra to compute our cross spectrum

% normalize spectrum by samples and sampling rate

freqs=fs/2*linspace(0,1,nfft/2+1);

% smooth with multi-taper

n=samples;
[tapers,lambda]=dpss(n,w,ntapers);

%ntapers
%w
% pre allocate cell arrays to store taper and trial estimates

cross_spect_mean=zeros(1,nfft);
spect1_mean=zeros(1,nfft);
spect2_mean=zeros(1,nfft);

store_cross=zeros(trials,nfft);
store_spect1=zeros(trials,nfft);
store_spect2=zeros(trials,nfft);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% AVERAGE SPECTRA %%%%%%%%%%%%%%%%%%%%

counter=1;
for i=1:trials

	% demean the data and compute the cross spectrum with the fields

	sig1=detrend(SIGNAL1(i,:));
	sig2=detrend(SIGNAL2(i,:));

	% multi-taper estimate
	% index by trial and taper, then we can estimate jackknife error bars

	cross_spect_tmp=zeros(1,nfft);
	spect1_tmp=zeros(1,nfft);
	spect2_tmp=zeros(1,nfft);

	for j=1:ntapers

		spect1=fft(sig1.*tapers(:,j)',nfft);
		spect2=fft(sig2.*tapers(:,j)',nfft);

		%[spect1]=spectrogram(sig1,n,overlap,nfft);
		%[spect2]=spectrogram(sig2,n,overlap,nfft);

		cross=spect1.*conj(spect2);
		
		power1=spect1.*conj(spect1);
		power2=spect2.*conj(spect2);
		
		cross_spect_tmp=cross_spect_tmp+cross./(ntapers);
		spect1_tmp=spect1_tmp+power1./(ntapers);
		spect2_tmp=spect2_tmp+power2./(ntapers);

		store_cross(counter,:)=cross;
		store_spect1(counter,:)=power1;
		store_spect2(counter,:)=power2;

		counter=counter+1;

	end

	cross_spect_mean=cross_spect_mean+cross_spect_tmp./trials;
	spect1_mean=spect1_mean+spect1_tmp./trials;	
	spect2_mean=spect2_mean+spect2_tmp./trials;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% STATS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

spect1=spect1_mean./fs;
spect2=spect2_mean./fs;
spect1=spect1(1:nfft/2+1);
spect2=spect2(1:nfft/2+1);

dof=ntapers*trials;
coh=cross_spect_mean./sqrt(spect1_mean.*spect2_mean);
coh=coh(1:(nfft/2+1));
abscoh=abs(coh);

% asymptotic errorbars

abscoh=abscoh(1:nfft/2+1);

m_coherence_err=1.96./(sqrt(dof));

upperconf=abscoh+m_coherence_err;
lowerconf=abscoh-m_coherence_err;

% standard z-transform (beta parameter from Pesaran et al. 2008)

q=sqrt(-(dof-2)*log(1-abscoh.^2));

% z-transform, coherence now in units of standard deviation from null

stats.zcoh=beta*(q-beta);
stats.alpha=alpha;
stats.var=(1.96)./(sqrt(dof));
stats.nullline=sqrt(1-(alpha).^(1/(((2*dof)/2)-1)));

%%%%%
