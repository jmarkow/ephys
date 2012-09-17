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

method='neg';
mpca=0; % enable at your own risk

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

% center and scale the coefficients

[coeffs,trials]=size(wavecoef);
to_del=[];

nbins=sqrt(trials);

for i=1:coeffs

	testpoints=wavecoef(i,:);

	samplemedian=median(testpoints);
	samplevar=median(abs(testpoints-samplemedian))./.6745;

	% normalizing does not seem to make much sense to me...maybe we should leave out for now

	testpoints=testpoints-samplemedian;
	wavecoef(i,:)=testpoints;

	samplemedian=median(testpoints);

	% need to compute deviation from normality take max of empirical cdf and normcdf 
	% with same mean and variance

	if any(isnan(testpoints)) | nbins==0
		ks_stat(i)=0;
		negentropy(i)=0;
		coeffbimodal(i)=0;
		to_del=[to_del i];
		continue;
	end

	% get the empirical pdf, kernel density estimate

	x=linspace(min(testpoints),max(testpoints),nbins);

	count=ksdensity(testpoints,x);

	fx_pdf=count./sum(count);
	fx_cdf=cumsum(fx_pdf);

	% evaluate normpdf

	gx_pdf=normpdf(x,samplemedian,samplevar);
	gx_pdf=gx_pdf./sum(gx_pdf);
	gx_cdf=cumsum(gx_pdf);

	% deviation, ks statistic

	ks_stat(i)=max(abs(fx_cdf-gx_cdf));

	normentropy=-sum(gx_pdf.*log2(gx_pdf+eps));
	obsentropy=-sum(fx_pdf.*log2(fx_pdf+eps));

	negentropy(i)=[normentropy-obsentropy];

	% is obsentropy==0 then the distribution has collapsed onto zero

	if negentropy(i)==normentropy
		negentropy(i)=NaN;
	end

	% an alternative measure

	% coefficient of bimodality...
	% consider the coefficient of bimodality here... should incorporate into compiled auto clustering

	coeffbimodal(i)=(1+skewness(testpoints)^2)/(kurtosis(testpoints)+((3*(trials-1)^2)/((trials-2)*(trials-3))));

	if abs(skewness(testpoints)>.4)
		coeffbimodal(i)=0;
	end


end

% remove any NaNs

ks_stat(ks_stat==1)=0;
ks_stat(isnan(ks_stat))=0;
coeffbimodal(isnan(coeffbimodal))=0;
negentropy(isnan(negentropy))=-inf;

% remove stray points

ks_stat(to_del)=[];
coeffbimodal(to_del)=[];
negentropy(to_del)=[];
wavecoef(to_del,:)=[];

% what's left?

[coeffs,trials]=size(wavecoef);

switch lower(method(1))

	case 'k'
		[val,loc]=sort(ks_stat,'descend');
		weights=ks_stat;
	case 'n'
		[val,loc]=sort(negentropy,'descend');
		weights=negentropy;
	case 'b'
		[val,loc]=sort(coeffbimodal,'descend');
		weights=coeffbimodal;
	otherwise
		error('ephysPipeline:getwaveletcoeffs:badmethod','Did not understand coefficient sorting method');
end


% mpca works surprisingly well, just scale the data by a multi-modality measure, let it rip...

if mpca

	% take components, weight by KS (or CB) and take PCA
	% scale KS test by L2 norm, use as a weight for coefficient

	% remove any negative values (only need to worry about negentropy here)

	weights(weights<0)=0;

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



