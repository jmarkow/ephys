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

projections=5; % number of dimensions to return, should retest for 
         % negentropy of coefficient of bimodality

maxiter=1000;
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

%sphered_data=sphered_data';

% return to dimensions x observations

for i=1:projections+1

	w=randn(m,1);
	w=w/norm(w);

	for j=1:maxiter

		y=w'*sphered_data;

		% add options for other functions

		leftside=mean(sech(y).^2).*w;

		% calculate the gradient

		% derivative on w

		size(y)
		pause();

		rightside=mean(sphered_data.*repmat(tanh(y),m,1),2);
	
		wplus=leftside-rightside;
		
		wnew=wplus/norm(wplus);

		tmp=[];

		% orthogonalize
		
		for k=1:i-1

			% inner product of wnew and other vectors
			% .* COEFFS(:,k)

			tmp(:,k)=(wnew'*COEFFS(:,k)).*COEFFS(:,k);
		end

		if i>1

			wnew=wnew-sum(tmp,2);
			wnew=wnew/norm(wnew);	
		end

	
		i	
		[w wnew]
		pause();

		1-(sum(wnew.*w))

		if 1-(sum(wnew.*w))<epsilon
			w=wnew;
			break;
		end

		w=wnew;


	end

	COEFFS(:,i)=w;

end

% the scores are simply the inner product of each sample with the projection matrix
% arrange scores similar to princomp, observations x dims

SCORES=zeros(m,projections);

% take the projections with maximum non-normality

