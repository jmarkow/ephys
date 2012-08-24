function [COEFFS SCORES]=pp_negentropy(DATA,varargin)
% simple projection pursuit (core of fixed point ICA algorithm)
% NOTE:  this function has not been tested and should be considered INCOMPLETE
%
%
%

% data should be samples x channels


% negentropy based projection pursuit
% usually only the first 2-3 dimensions are meaningful

% first, sphere the data

projections=2; % number of dimensions to return, should retest for 
         % negentropy of coefficient of bimodality

maxiter=400;
epsilon=1e-4;

nparams=length(varargin);

if mod(nparams,2)>0
	error('ephysPipeline:argChk','Parameters must be specified as parameter/value pairs!');
end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'maxiter'
			maxiter=varargin{i+1};
		case 'projections'
			projections=varargin{i+1};
		case 'epsilon'
			epsilon=varargin{i+1};
	end
end

DATA=DATA';

% now observations x dimensions

[m,n]=size(DATA);
data_mean=mean(DATA);

% center the data

demeaned_data=DATA-ones(m,1)*data_mean;
sphered_data=demeaned_data./repmat(std(demeaned_data),m,1);

sphered_data=sphered_data';

% return to dimensions x observations

for i=1:projections+1

	w=randn(1,m);
	w=w/norm(w);

	for j=1:maxiter

		y=sphered_data*w';

		% add options for other functions

		g=repmat(y.^3,1,m);
		gprime=repmat(3.*y.^2,1,m);
		
		% calculate the gradient

		leftside=mean(sphered_data.*g);

		% derivative on w

		rightside=mean(gprime).*w;
		
		wplus=leftside-rightside;
		wnew=wplus/norm(wplus);

		tmp=[];

		% orthogonalize
		
		for k=1:i-1

			% inner product of wnew and other vectors
			% .* COEFFS(:,k)

			tmp(k,:)=(COEFFS(k,:)*wnew').*COEFFS(k,:);
		end

		if i>1
			wnew=wnew-sum(tmp,1);
			wnew=wnew/norm(wnew);	
		end

		if 1-(sum(wnew.*w))<epsilon
			w=wnew;
			break;
		end

		w=wnew;


	end

	COEFFS(i,:)=w;

end

% the scores are simply the inner product of each sample with the projection matrix
% arrange scores similar to princomp, observations x dims

SCORES=zeros(m,projections);

% take the projections with maximum non-normality

for i=1:projections+1

	testpoints=COEFFS(i,:);
	[fx,x]=ecdf(testpoints);

	% evaluate normal cdf at same mean and variance as data
	% use robust mean and variance estimators per Takekawa et al. (2012)

	samplemedian=median(testpoints);

	% robust variance

	samplevar=median(abs(testpoints-samplemedian))./.6745;

	% evaluate normcdf

	gx=normcdf(x,samplemedian,samplevar);

	% deviation, ks statistic

	% negentropy, has the entropy been reduced relative to the normal distribution?

	normentropy=-sum(gx.*log2(gx+eps));
	negentropy(i)=[normentropy-obsentropy]; % an alternative measure

end

% rank dimensions by negentropy

% TODO:  finish rank sorting, perhaps put in exclusion criteria

[val,loc]=sort(negentropy,'descending');
rank_proj=loc(1:projections)


for i=1:n
	for j=1:projections
		SCORES(i,j)=COEFFS(j,:)*sphered_data(i,:)';
	end
end

