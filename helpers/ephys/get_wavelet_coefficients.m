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
mpca=1; % mpca improves results dramatically

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
to_del=[];

for i=1:coeffs

	% remove points > +/- 3 std per quiroga et al 2004, probably irrelevant for
	% coeff of bimodality

	% first normalize the data (zscore)

	testpoints=wavecoef(i,:);

	samplemedian=median(testpoints);
	samplevar=median(abs(testpoints-samplemedian))./.6745;
	testpoints=(testpoints-samplemedian)./samplevar;

	% store the normalized data

	wavecoef(i,:)=testpoints;

	% need to compute deviation from normality take max of empirical cdf and normcdf 
	% with same mean and variance

	if any(isnan(testpoints))
		normdev(i)=0;
		negentropy(i)=0;
		coeffbimodal(i)=0;
		to_del=[to_del i];
		continue;
	end

	% get the empirical cdf

	[fx,x]=ecdf(testpoints);

	% evaluate normal cdf at same mean and variance as data
	% use robust mean and variance estimators per Takekawa et al. (2012)

	samplemedian=median(testpoints);

	% robust variance

	samplevar=median(abs(testpoints-samplemedian))./.6745;

	% evaluate normcdf

	gx=normcdf(x,samplemedian,samplevar);

	% deviation, ks statistic

	normdev(i)=[max(abs(fx-gx))];

	% negentropy, has the entropy been reduced relative to the normal distribution?

	normentropy=-sum(gx.*log2(gx+eps));
	obsentropy=-sum(fx.*log2(fx+eps));
	negentropy(i)=[normentropy-obsentropy]; % an alternative measure

	% coefficient of bimodality...
	% consider the coefficient of bimodality here... should incorporate into compiled auto clustering

	coeffbimodal(i)=(1+skewness(testpoints)^2)/(kurtosis(testpoints)+3);

end

normdev(normdev==1)=0;
normdev(isnan(normdev))=0;
coeffbimodal(isnan(coeffbimodal))=0;

normdev(to_del)=[];
coeffbimodal(to_del)=[];
negentropy(to_del)=[];
wavecoef(to_del,:)=[];

[coeffs,trials]=size(wavecoef);

switch lower(method(1))

	case 'k'
		[val,loc]=sort(normdev,'descend');
		weights=normdev;
	case 'n'
		[val,loc]=sort(negentropy,'descend');
		weights=negentropy;
	case 'b'
		[val,loc]=sort(coeffbimodal,'descend');
		weights=coeffbimodal;
	otherwise
		error('ephysPipeline:getwaveletcoeffs:badmethod','Did not understand coefficient sorting method');
end

if mpca

	% take components, weight by KS (or CB) and take PCA
	% scale KS test by L2 norm, use as a weight for coefficient

	for i=1:coeffs
		wavecoeff(i,:)=(weights(i)/sqrt(sum(wavecoef(i,:).^2))).*wavecoef(i,:);
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



