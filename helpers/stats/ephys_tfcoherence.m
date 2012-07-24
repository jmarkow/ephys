function [COH,T,F,SPECT1,SPECT2]=ephys_tfcoherence(SIGNAL1,SIGNAL2,varargin)
%
%
%
%

nparams=length(varargin);

nfft=10000;
n=6250;
overlap=6000;
min_f=15;
max_f=100;
w=2;

ntapers=[];
freq_range=[300]; % let's just refilter

min_f=15;
max_f=100;

sr=25e3;

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
		case 'sr'
			sr=varargin{i+1};
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

if isempty(nfft)
	nfft=max([n 2^nextpow2(n)]);
else
	nfft=2^nextpow2(nfft);
end

resolution=w*1/(n/sr);
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

if isempty(ntapers)
    ntapers=2*(w)-1;
end

[tapers,lambda]=dpss(n,w,ntapers);

if mod(nfft,2)==0
	rows=(nfft/2+1);
else
	rows=(nfft+1)/2;
end

columns=fix((nsamples-overlap)/(n-overlap));

% axes

F=((1:rows)./rows).*(sr/2);
col_idx=1+(0:(columns-1))*(n-overlap);
T=((col_idx-1)+((n/2)'))/sr;

cross_spect_mean=zeros(rows,columns);
spect1_mean=zeros(rows,columns);
spect2_mean=zeros(rows,columns);

parfor i=1:ntrials

	disp(['Trial ' num2str(i)]);

	curr1=SIGNAL1(i,:);
	curr2=SIGNAL2(i,:);

	% spectra

	% multi-taper estimate

	% accumulate across tapers to cut down compt time, then we can use parfor

	cross_spect_tmp=zeros(rows,columns);
	spect1_tmp=zeros(rows,columns);
	spect2_tmp=zeros(rows,columns);

	for j=1:ntapers

		spect1=spectrogram(curr1,tapers(:,j)',overlap,nfft);
		spect2=spectrogram(curr2,tapers(:,j)',overlap,nfft);

		% take absolute value of the cross power after trial and taper
		% averaging

		% don't need the normalizations

		cross=(spect1.*conj(spect2));

		% take power of two signals

		power1=(abs(spect1).^2);
		power2=(abs(spect2).^2);

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

SPECT1=spect1_mean;
SPECT2=spect2_mean;
COH=cross_spect_mean./sqrt(spect1_mean.*spect2_mean);

% display is user selects to display

