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

method='bimodal';

for i=1:2:nparams
	switch lower(varargin{i})
		case 'method'
			method=varargin{i+1};
	end
end

[samples,trials]=size(WAVEFORMS);

% may want to consider a new family here, perhaps coiflet instead of Haar

parfor i=1:trials
	[wavecoef(:,i),l]=wavedec(WAVEFORMS(:,i),5,'haar'); % 4-level wavelet decomposition
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
		normdev(i)=-inf;
		negentropy(i)=-inf;
		coeffbimodal(i)=-inf;
		continue;
	end

	[fx,x]=ecdf(testpoints);

	% evaluate normal cdf at same mean and variance as data

	samplemean=mean(testpoints);
	samplestd=std(testpoints);

	gx=normcdf(x,samplemean,samplestd);

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
coeffbimodal(isnan(coeffbimodal))=0;

if strcmp(lower(method),'ks')
	[val,loc]=sort(normdev,'descend');
elseif strcmp(lower(method),'neg')
	[val,loc]=sort(negentropy,'descend');
else

	[val,loc]=sort(coeffbimodal,'descend');

	% strip anything below .3

	cutoff=find(val<.3);

	% if everything is below .333 then include it all

	if length(cutoff)==length(val)
		cutoff=[];
	end

	val(cutoff)=[];
	loc(cutoff)=[];

end

if length(loc)>=NCOEFFS
	sortcoeffs=loc(1:NCOEFFS);
else
	sortcoeffs=loc;
end

for i=1:length(sortcoeffs)
	SELCOEFFS(:,i)=wavecoef(sortcoeffs(i),:);
end


