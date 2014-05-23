function [newmodel]=free_gmem(DATA,INIT,NCLUST,varargin)
%EM for Gaussian mixture with support for outlier distribution
%
%
%
% all credit to Daniel Wagenaar's code
% http://www.danielwagenaar.net/res/papers/00-Wage2.pdf
%
% algorithm derived from Maneesh Sahani's thesis

% data should be passed as a matrix observations x variables

if nargin<3
	NCLUST=2;
end

[datapoints,D]=size(DATA);
nparams=length(varargin);

if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs');
end

%%%

maxiter=100;
regularize=1e-6;
epsilon=1e-10;
lambda=.05; % changed from .01 to .05 9/19/13
garbage=1;
merge=1;
splitepsi=2; % noise scale for splits
maxcand=20; % maximum number of smem candidates
smemiter=20; % maximum smem iterations
debug=0;

for i=1:2:nparams
	switch lower(varargin{i})
		case 'regularize'
			regularize=varargin{i+1};
		case 'maxiter'
			maxiter=varargin{i+1};
		case 'lambda'
			lambda=varargin{i+1};
		case 'epsilon'
			epsilon=varargin{i+1};
		case 'maxiter'
			maxiter=varargin{i+1};
		case 'garbage'
			garbage=varargin{i+1};
		case 'beta'
			beta=varargin{i+1};
		case 'merge'
			merge=varargin{i+1};
		case 'debug'
			debug=varargin{i+1};
		otherwise
	end
end

if nargin<2 | isempty(INIT)
	INIT=randinit(DATA,NCLUST,regularize);
end

mu=INIT.mu;
mixing=INIT.mixing;
sigma=INIT.sigma;

% make sure dimensions are correct

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% CHECK DIMENSIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if size(mu,2)~=D
	error('Mean dimensionality %g does not match data point dimensions %g',size(mu,2),D);
end

if length(mixing)~=NCLUST
	error('Wrong number of mixing proportions %g, NCLUST=%g',length(mixing),NCLUST);
end

if size(sigma,1)~=size(sigma,2)
	error('Covariance matrix must have the same dimensions...');
end

if size(sigma,3)~=NCLUST
	error('Must have covariance matrices for all clusters...');
end

if size(sigma,1)~=D
	error('Covariance matrix has wrong dimensionality %g, D=%g',size(sigma,1),D);
end

% number of parameters (full covariance)

nparams=(NCLUST*D*((D+1)/2));
nparams=nparams+NCLUST-1+NCLUST*D;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% FULL EM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% return the parameters


if garbage
	[unip,mixing]=initgarbage(DATA,NCLUST,INIT);
	INIT.mixing=mixing;
else
	unip=[];
end

if debug
	fig=figure();
	gaussvis(INIT,DATA,'fig_num',fig);
	pause(.5);
end

[newmodel]=fullem(DATA,INIT,unip,maxiter,epsilon,lambda);

if debug
	gaussvis(newmodel,DATA,'fig_num',fig);
	drawnow;
	pause(.5);
end

newmodel.garbage=zeros(1,size(newmodel.R,2));

if garbage
	newmodel.garbage(end)=1;
end

% get the merge candidates
% perform SMEM

newmodel=compute_icl(DATA,newmodel);
mergemodel1=newmodel;
	
% keep merging until BIC no longer improves

state=1; % 1 for merge, 2 for split, 3 for quit
fail=0;

% initialize the model
for ii=1:smemiter
	
	state=1;
	fail=0;
	%mergemodel1=newmodel;
	NCLUST=size(mergemodel1.mu,1);

	while fail<2

		% get the merge candidates

		if state==1

			% can't merge with only 1 cluster

			if NCLUST==1
				state=2;
				fail=fail+1;
				continue;
			end

			[~,pairs]=mergemerit(mergemodel1);
			npairs=size(pairs,1);

			% try pairs until we get a better result

			for i=1:min(maxcand,npairs)

				currpair=pairs(i,:);

				fprintf(1,'Merging %g and %g\n',currpair(1),currpair(2));
				
				currpair=pairs(i,:);
				
				tmpmodel=mergemodel1;
				tmpmodel=mergeclust(mergemodel1,currpair(1),currpair(2));

				tmpmodel=partialem(DATA,tmpmodel,unip,maxiter,...
					epsilon,lambda,currpair(1));
				tmpmodel=fullem(DATA,tmpmodel,unip,maxiter,epsilon,lambda);	

				if garbage
					tmpmodel.garbage=mergemodel1.garbage;	
					tmpmodel.garbage(currpair(2))=[];
				end

				tmpmodel=compute_icl(DATA,tmpmodel);

				fprintf(1,'Old icl %g, new icl %g\n',mergemodel1.ICL,tmpmodel.ICL);

				if tmpmodel.ICL<mergemodel1.ICL
					break;
				end
			end

			if tmpmodel.ICL<mergemodel1.ICL
				mergemodel1=tmpmodel;
				NCLUST=NCLUST-1;
				state=1;
				fail=0;
			else
				state=2;
				fail=fail+1;
			end

		else

			[~,splits]=splitmerit(DATA,mergemodel1,unip);
			nsplits=length(splits);

			for i=1:min(maxcand,nsplits)

				fprintf(1,'Splitting %g\n',splits(i));

				tmpmodel=splitclust(mergemodel1,splits(i),NCLUST+1,splitepsi);
				tmpmodel=partialem(DATA,tmpmodel,unip,maxiter,epsilon,lambda,[splits(i) NCLUST+1]);
				tmpmodel=fullem(DATA,tmpmodel,unip,maxiter,epsilon,lambda);
				tmpmodel=compute_icl(DATA,tmpmodel);

				if garbage
					tmpmodel.garbage=zeros(1,NCLUST+2);
					tmpmodel.garbage(end)=1;
				end

				fprintf(1,'Old icl %g, new icl %g\n',mergemodel1.ICL,tmpmodel.ICL);

				if tmpmodel.ICL<mergemodel1.ICL
					break;
				end

			end

			if tmpmodel.ICL<mergemodel1.ICL
				mergemodel1=tmpmodel;
				NCLUST=NCLUST+1;
				state=1; % set state to merge
				fail=0;
			else
				state=1; % go back to merge
				fail=fail+1;
			end

		end

		% run partial em on the new cluster

	end

end

newmodel=mergemodel1;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% RANDOM INITILIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function NEWMODEL=randinit(DATA,NCLUST,regularize)

[datapoints,D]=size(DATA);

% get the variance of all dimensions, set up
% a diagonal covariance

% randinit

fprintf(1,'Initializing %g cluster(s) with random points\n',NCLUST);

initpoints=randsample(datapoints,NCLUST);

NEWMODEL.mu=DATA(initpoints,:);

% initialize all covariance matrices

datavar=var(DATA);
initsigma=diag(datavar);

for i=1:NCLUST
	NEWMODEL.sigma(:,:,i)=initsigma+eye(D).*regularize;
	NEWMODEL.mixing(i)=1/NCLUST;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% GARBAGE UNIFORM DENSITY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [P,mixing]=initgarbage(DATA,NCLUST,MODEL)

mixing=MODEL.mixing;
[datapoints,D]=size(DATA);

% set the uniform density

datarange=range(DATA);

% attempt to compute the volume of the space

p=prod(datarange);

if p==inf

	% issue warning here

	warning('gmem:volumetoolarge',...
		'Cannot compute volume, setting to : %e',p);
	p=1e30;
end

P=(1/p).*ones(datapoints,1);

% set the initial mixing proportion

mixing(end+1)=1/(NCLUST);

% renormalize the mixing proportions

mixing=mixing./sum(mixing);

end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% FULL EM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [NEWMODEL]=fullem(DATA,MODEL,unip,maxiter,epsilon,lambda)

%
%
%
%
%

[datapoints,D]=size(DATA);

mu=MODEL.mu;
sigma=MODEL.sigma;
mixing=MODEL.mixing;
NCLUST=size(mu,1);

garbage=0;
if ~isempty(unip)
	garbage=1;
end

prev_likelihood=1e-9;

for i=1:maxiter

	% compute likelihoods
	% responsibilities

	if garbage
		R=zeros(datapoints,NCLUST+1);
	else
		R=zeros(datapoints,NCLUST);
	end

	px=zeros(datapoints,1);
	den=zeros(datapoints,1);

	% e step, get the responsibilities

	for j=1:NCLUST

		pxtheta=mvnpdf(DATA,mu(j,:),sigma(:,:,j)); % point x 1 vector
		mixprob=mixing(j)*pxtheta;
		px=px+mixprob;
		den=den+mixprob;
		R(:,j)=mixprob;

	end

	if garbage
		mixprob=mixing(NCLUST+1)*unip;
		R(:,NCLUST+1)=mixprob;
		den=den+mixprob;
		px=px+mixprob;
	end

	% likelihood

	for j=1:NCLUST
		R(:,j)=R(:,j)./(den+1e-300);
	end

	if garbage
		R(:,NCLUST+1)=R(:,NCLUST+1)./(den+1e-300);
	end

	mixing=mean(R);
	likelihood=sum(log(px+1e-300));

	deltalikelihood=(likelihood-prev_likelihood);

	% break if we've reached our stopping criterion

	if deltalikelihood>=0 && deltalikelihood<epsilon*abs(prev_likelihood)
		break;
	end

	prev_likelihood=likelihood;
	
	% update mu, sigma and mixing probabilities

	for j=1:NCLUST

		% need the total r for normalization

		totalR=sum(R(:,j)); 

		% recompute mu, the inner product between all datapoints and their
		% responsibilities within the cluster, normalized by the totalR

		mu(j,:)=(DATA'*R(:,j))./(totalR+1e-300); 

		% get the deviation from the mean for each datapoint

		dx=(DATA-repmat(mu(j,:),[datapoints 1]))';

		% transpose so we have D x datapoints

		Rdx=repmat(R(:,j)',[D 1]).*dx;

		% now R for the cluster is repeated so we have D x datapoints

		% take the inner product between the mean deviation D x datapoints and datapoints x D 
		% for responsibilities

		% add the regularization constant and normalize

		sigma(:,:,j)=(Rdx*dx'+lambda*eye(D))/(totalR+lambda);
	end

	% store the likelihood

end

NEWMODEL.R=R;
NEWMODEL.sigma=sigma;
NEWMODEL.mu=mu;
NEWMODEL.mixing=mixing;
NEWMODEL.likelihood=likelihood;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% PARTIAL EM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [NEWMODEL]=partialem(DATA,MODEL,unip,maxiter,epsilon,lambda,idx)

%
%
%
%
%

garbage=0;
if ~isempty(unip)
	garbage=1;
end

% total responsibilities for normalization

mu=MODEL.mu;
sigma=MODEL.sigma;
mixing=MODEL.mixing;

R=estep(DATA,MODEL,unip);

NCLUST=size(mu,1);

if garbage
	normalizationR=sum(R(:,[idx NCLUST+1]),2);
else
	normalizationR=sum(R(:,[idx]),2);
end

[datapoints,D]=size(DATA);
prev_likelihood=1e-9;

for i=1:maxiter

	% compute likelihoods
	% responsibilities

	if garbage
		R=zeros(datapoints,NCLUST+1);
	else
		R=zeros(datapoints,NCLUST);
	end

	px=zeros(datapoints,1);
	den=zeros(datapoints,1);

	% e step, get the responsibilities

	for j=idx

		pxtheta=mvnpdf(DATA,mu(j,:),sigma(:,:,j)+1e-6*eye(D)); % point x 1 vector
		mixprob=mixing(j)*pxtheta;
		px=px+mixprob;
		den=den+mixprob;
		R(:,j)=mixprob;

	end

	if garbage
		mixprob=mixing(NCLUST+1)*unip;
		R(:,NCLUST+1)=mixprob;
		den=den+mixprob;
		px=px+mixprob;
	end

	% likelihood

	den=normalizationR./(den+1e-300);

	for j=idx
		R(:,j)=R(:,j).*den;
	end

	if garbage
		R(:,NCLUST+1)=R(:,NCLUST+1).*den;
	end

	likelihood=sum(log(px+1e-300));
	deltalikelihood=(likelihood-prev_likelihood);

	% break if we've reached our stopping criterion

	if deltalikelihood>=0 && deltalikelihood<epsilon*abs(prev_likelihood)
		break;
	end

	prev_likelihood=likelihood;
	% update mu, sigma and mixing probabilities

	for j=idx

		% need the total r for normalization

		totalR=sum(R(:,j)); 

		% recompute mu, the inner product between all datapoints and their
		% responsibilities within the cluster, normalized by the totalR

		mu(j,:)=(DATA'*R(:,j))./totalR; 

		% get the deviation from the mean for each datapoint

		dx=(DATA-repmat(mu(j,:),[datapoints 1]))';

		% transpose so we have D x datapoints

		Rdx=repmat(R(:,j)',[D 1]).*dx;

		% now R for the cluster is repeated so we have D x datapoints

		% take the inner product between the mean deviation D x datapoints and datapoints x D 
		% for responsibilities

		% add the regularization constant and normalize

		sigma(:,:,j)=(Rdx*dx'+lambda*eye(D))/(totalR+lambda);
		mixing(j)=mean(R(:,j));
	end

	if garbage
		mixing(NCLUST+1)=mean(R(:,NCLUST+1));
	end

end

NEWMODEL.R=R;
NEWMODEL.sigma=sigma;
NEWMODEL.mu=mu;
NEWMODEL.mixing=mixing;
NEWMODEL.likelihood=likelihood;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% MERGE MERIT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [merit,pairs]=mergemerit(MODEL)

%
%
%
%

clusterids=find(~MODEL.garbage);
pairs=nchoosek(clusterids,2);

for i=1:size(pairs,1)
	merit(i)=MODEL.R(:,pairs(i,1))'*MODEL.R(:,pairs(i,2));
	merit(i)=merit(i)./(norm(MODEL.R(:,pairs(i,1)))*norm(MODEL.R(:,pairs(i,2))));
end

[val,idx]=sort(merit,'descend');
merit=merit(idx);
pairs=pairs(idx,:);


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% SPLIT MERIT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [split,splitidx]=splitmerit(DATA,MODEL,unip)
%
%
%

clusterids=find(~MODEL.garbage);

R=estep(DATA,MODEL,unip);

% take the responsibilities, compare with the model PDFs

mu=MODEL.mu;
sigma=MODEL.sigma;

for i=clusterids

	pxtheta=mvnpdf(DATA,mu(i,:),sigma(:,:,i)); % point x 1 vector
	
	% "empirical" density
	
	f=R(:,i)./sum(R(:,i));

	pxtheta=pxtheta+1e-5;
	idx=find(f>1e-5);

	% get KL divergence

	split(i)=sum(f(idx).*log(f(idx)./pxtheta(idx)));

end

[split splitidx]=sort(split,'descend');

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% MERGING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [NEWMODEL]=mergeclust(MODEL,c1,c2)
%
%
%
%


% merge into c1


mu=MODEL.mu;
mixing=MODEL.mixing;
sigma=MODEL.sigma;

mu(c1,:)=(mixing(c1)*mu(c1,:)+mixing(c2)*mu(c2,:))./(mixing(c1)+mixing(c2));
mixing(c1)=mixing(c1)+mixing(c2);
sigma(:,:,c1)=(sigma(:,:,c1)+sigma(:,:,c2))./2;

% set the mixing proportion of the merged cluster to 0

mu(c2,:)=[];
mixing(c2)=[];
sigma(:,:,c2)=[];

NEWMODEL.mu=mu;
NEWMODEL.sigma=sigma;
NEWMODEL.mixing=mixing;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% SPLITTING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [NEWMODEL]=splitclust(MODEL,sc,sc1,splitepsi)


% take the largest eigenvalue of the covariance

mu=MODEL.mu;
mixing=MODEL.mixing;
sigma=MODEL.sigma;

[nclust,d]=size(mu);
len=length(mixing);

if len>nclust
	garbage_mix=MODEL.mixing(len);
end

[U DD V]=svd(MODEL.sigma(:,:,sc));

% put the new cluster in the beginning

mixing(sc)=mixing(sc)/2;
mixing(sc1)=mixing(sc);

if len>nclust
	mixing(sc1+1)=garbage_mix;
end

sd=sqrt(diag(d));

mu(sc,:)=mu(sc,:)+(splitepsi*U*(sd*randn(d,1)))';
mu(sc1,:)=mu(sc,:)+(splitepsi*U*(sd*randn(d,1)))';

sigma(:,:,sc)=DD(1,1)*eye(d);
sigma(:,:,sc1)=sigma(:,:,sc);

NEWMODEL.mu=mu;
NEWMODEL.sigma=sigma;
NEWMODEL.mixing=mixing;


end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% ESTEP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function R=estep(DATA,MODEL,unip)
%
%
%
%
%
%
%
% performs the e step

[datapoints,D]=size(DATA);

garbage=0;
if ~isempty(unip)
	garbage=1;
end

mu=MODEL.mu;
sigma=MODEL.sigma;
mixing=MODEL.mixing;

NCLUST=size(mu,1);

px=zeros(datapoints,1);
den=zeros(datapoints,1);

if ~garbage
	R=zeros(datapoints,NCLUST);
else
	R=zeros(datapoints,NCLUST+1);
end

for i=1:NCLUST

	pxtheta=mvnpdf(DATA,mu(i,:),sigma(:,:,i)); % point x 1 vector
	mixprob=mixing(i)*pxtheta;
	px=px+mixprob;
	den=den+mixprob;
	R(:,i)=mixprob;

end

if garbage

	mixprob=mixing(NCLUST+1)*unip;
	R(:,NCLUST+1)=mixprob;
	den=den+mixprob;
	px=px+mixprob;

end

for i=1:NCLUST
	R(:,i)=R(:,i)./(den+1e-300);
end

if garbage
	R(:,NCLUST+1)=R(:,NCLUST+1)./(den+1e-300);
end


end

function MODEL=compute_icl(DATA,MODEL)

[datapoints,d]=size(DATA);

nclust=size(MODEL.mu,1);

nparams=(nclust*d*((d+1)/2));
nparams=nparams+nclust-1+nclust*d;

MODEL.BIC=-2*MODEL.likelihood+log(datapoints)*nparams;

% get the total entropy

entropy=zeros(1,nclust);
for i=1:nclust
	pxtheta=mvnpdf(DATA,MODEL.mu(i,:),MODEL.sigma(:,:,i));
	tmp=MODEL.R(:,i).*pxtheta;
	tmp(tmp<0)=[];
	entropy(i)=sum(tmp.*log(tmp+1e-300));
end

MODEL.ICL=MODEL.BIC-2*sum(entropy);

end
