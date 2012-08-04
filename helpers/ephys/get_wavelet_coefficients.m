function SELCOEFFS=get_wavelet_coefficients(WAVEFORMS,NCOEFFS,varargin)
%gets wavelet coefficients and sorts according to ks stat, negentropy, or coefficient of bimodality
%
%
%
%

% sort either according to negentropy or KS test (deviation from normality)

nparams=length(varargin);

if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs!');
end

method='ks'; % ks with mpca seem to be the best choice here
mpca=1;

for i=1:2:nparams
	switch lower(varargin{i})
		case 'method'
			method=varargin{i+1};
		case 'mpca'
			mpca=varargin{i+1};
	end
end

[samples,trials]=size(WAVEFORMS);

% tried db8, sym7, and coif3, haar is still the most reliable!

parfor i=1:trials
	[wavecoef(:,i),l]=wavedec(WAVEFORMS(:,i),4,'haar'); % 4-level wavelet decomposition
end

[coeffs,trials]=size(wavecoef);

for i=1:coeffs

	% remove points > +/- 3 std per quiroga et al 2004, probably irrelevant for
	% coeff of bimodality

	testpoints=wavecoef(i,:);

	%cutoff=3*std(testpoints);

	% apparently overlapping spikes create outliers that corrupt the KS stat
	% this doesn't appear to work well in practice, let's stick with 
	% either negentropy or coefficient of bimodality

	%testpoints(testpoints<-cutoff)=[];
	%testpoints(testpoints>cutoff)=[];

	% need to compute deviation from normality take max of empirical cdf and normcdf 
	% with same mean and variance

	% samplecdf

	if length(testpoints(testpoints~=NaN))<1
		normdev(i)=0;
		negentropy(i)=0;
		coeffbimodal(i)=0;
		continue;
	end

	[fx,x]=ecdf(testpoints);

	% evaluate normal cdf at same mean and variance as data

	% use robust mean and variance estimators per Takekawa et al. (2012)

	samplemedian=median(testpoints);
	samplevar=median(abs(testpoints-samplemedian))./.6745;

	gx=normcdf(x,samplemedian,samplevar);

	% deviation, ks statistic

	normdev(i)=[max(abs(fx-gx))];

	normentropy=-sum(gx.*log2(gx+eps));
	obsentropy=-sum(fx.*log2(fx+eps));

	negentropy(i)=[normentropy-obsentropy]; % an alternative measure

	% consider the coefficient of bimodality here... should incorporate into compiled auto clustering

	coeffbimodal(i)=(1+skewness(testpoints)^2)/(kurtosis(testpoints)+3);

end

% with bimodality we could add a simple cutoff...if no one passes then let everyone in

normdev(normdev==1)=0;
normdev(isnan(normdev))=0;
coeffbimodal(isnan(coeffbimodal))=0;

switch lower(method(1))

	case 'k'
		[val,loc]=sort(normdev,'descend');
	case 'n'
		[val,loc]=sort(negentropy,'descend');
	case 'b'
		[val,loc]=sort(coeffbimodal,'descend');
	otherwise
		error('ephysPipeline:getwaveletcoeffs:badmethod','Did not understand coefficient sorting method');
end

if mpca
	% take components, weight by KS (or CB) and take PCA
	% scale KS test by L2 norm, use as a weight for coefficient

	for i=1:coeffs
		wavecoeff(i,:)=(normdev(i)/sqrt(sum(wavecoef(i,:).^2))).*wavecoef(i,:);
	end

	% get PCA scores

	[pcacoeff pcascore]=princomp(wavecoef');

	%multimodal PCA

	SELCOEFFS=pcascore(:,1:NCOEFFS);
else

	if length(loc)>=NCOEFFS
		sortcoeffs=loc(1:NCOEFFS);
	else
		sortcoeffs=loc;
	end

	for i=1:length(sortcoeffs)
		SELCOEFFS(:,i)=wavecoef(sortcoeffs(i),:);
	end
end



