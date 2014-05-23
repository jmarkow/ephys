function [newmodel]=gmem(DATA,INIT,NCLUST,varargin)
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
splitepsi=1; % noise scale for splits
maxcand=5; % maximum number of smem candidates
smemiter=100; % maximum smem iterations
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


if merge & NCLUST>2
	
	% keep merging until BIC no longer improves

	breakflag=0;

	for ii=1:smemiter

		% get the merge candidates

		[merges,pairs]=mergemerit(newmodel);

		% get the split candidates

		[~,splits]=splitmerit(DATA,newmodel,unip);

		% go through each pair find a non-matching split
		% now for each pair we take the split candidates ~==pair


		triplet=[];
		for i=1:size(pairs,1)

			currpair=pairs(i,:);

			for j=1:length(splits)
				if ~any(splits(j)==currpair)
					triplet=[triplet;currpair splits(j)];
				end
			end
		end

		for i=1:maxcand

			if i>size(triplet,1)
				break;
			end

			currtrip=triplet(i,:);
			fprintf(1,'Merging %g and %g, splitting %g\n',currtrip(1),currtrip(2),currtrip(3));


			mergemodel1=mergeclust(newmodel,...
				currtrip(1),currtrip(2));
			mergemodel1=splitclust(mergemodel1,...
				currtrip(3),currtrip(2),splitepsi);

			% run partial em on our new merged cluster

			mergemodel1=partialem(DATA,mergemodel1,...
				unip,maxiter,epsilon,lambda,[currtrip]);
			mergemodel1=fullem(DATA,mergemodel1,...
				unip,maxiter,epsilon,lambda);

			if debug
				gaussvis(mergemodel1,DATA,'fig_num',fig);
				drawnow;
				pause(.5);
			end

			mergemodel1.garbage=newmodel.garbage;


			if mergemodel1.likelihood<newmodel.likelihood
				
				fprintf(1,'No improvement in likelihood trying another candidate\n');
				continue;

			else

				fprintf(1,'Candidate improved, breaking out of inner SMEM loop\n');
				newmodel=mergemodel1;
				
				% break out of the candidate loop if we've improved

				break;
			end

		end	

		% if we've ended and there's no improvement, break out of the main loop

		if mergemodel1.likelihood<newmodel.likelihood

			fprintf(1,'No improvements found in any candidates, breaking out of outer SMEM loop\n');
			break;
		end

		% run partial em on the new cluster

	end

end

if debug
	gaussvis(newmodel,DATA,'fig_num',fig);
	title('Final');
	drawnow;
	pause(.5);
end


newmodel.BIC=-2*newmodel.likelihood+log(datapoints)*nparams;

% get the total entropy

for i=1:NCLUST
	pxtheta=mvnpdf(DATA,newmodel.mu(i,:),newmodel.sigma(:,:,i));
	tmp=newmodel.R(:,i).*pxtheta;
	tmp(tmp<0)=[];
	entropy(i)=sum(tmp.*log(tmp+1e-300));
end

newmodel.ICL=newmodel.BIC-2*sum(entropy);

% mml

nparams=1+D+.5*D*(D+1);

rightterm=-.5*nparams*sum(log((newmodel.mixing.*datapoints)/12))+...
	(NCLUST/2)*log(datapoints/12)+(NCLUST*(nparams+1))/2;
leftterm=newmodel.likelihood;
newmodel.MML=leftterm+rightterm;

%newmodel.MML2=.5*sum(log(newmodel.mixing))+((NCLUST*nparams+nparams)/2)*log(datapoints)-newmodel.likelihood;

fprintf(1,'NComponents %g, Likelihood %5.4e, BIC %5.4e, MML %5.4e, ICL %5.4e\n',NCLUST,...
	newmodel.likelihood,newmodel.BIC,newmodel.MML,newmodel.ICL);

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

[datapoints,D]=size(DATA);
NCLUST=size(mu,1);

if garbage
	normalizationR=sum(R(:,[idx NCLUST+1]),2);
else
	normalizationR=sum(R(:,[idx]),2);
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

	for j=idx

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

mixing(c2)=0;

NEWMODEL.mu=mu;
NEWMODEL.sigma=sigma;
NEWMODEL.mixing=mixing;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% SPLITTING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [MODEL]=splitclust(MODEL,sc,sc1,splitepsi)


% take the largest eigenvalue of the covariance

[NCLUST,D]=size(MODEL.mu);

[U DD V]=svd(MODEL.sigma(:,:,sc));

% put the new cluster in the beginning

MODEL.mixing(sc)=MODEL.mixing(sc)/2;
MODEL.mixing(sc1)=MODEL.mixing(sc);

sd=sqrt(diag(D));

MODEL.mu(sc,:);
MODEL.mu(sc,:)=MODEL.mu(sc,:)+(splitepsi*U*(sd*randn(D,1)))';
MODEL.mu(sc1,:)=MODEL.mu(sc,:)+(splitepsi*U*(sd*randn(D,1)))';

MODEL.sigma(:,:,sc)=DD(1,1)*eye(D);
MODEL.sigma(:,:,sc1)=MODEL.sigma(:,:,sc);

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

