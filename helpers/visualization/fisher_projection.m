function [DENSITY1,DENSITY2,XI]=fisher_projection(DATAMAT1,DATAMAT2,varargin)
%cluster statistics, include Fisher projection and other quality metrics
%
%
%

%
% perhaps include plot of stability across trials (maybe threshold?)
%

nparams=length(varargin);

if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs');
end

% just in case add the hot colormap at the end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'fig_num'
			fig_num=varargin{i+1};
	end
end

% compute the discriminant

mu1=mean(DATAMAT1)';
mu2=mean(DATAMAT2)';

cov1=cov(DATAMAT1,1);
cov2=cov(DATAMAT2,1);

cov1=cov1+eye(size(cov1)).*1e-20;
cov2=cov2+eye(size(cov2)).*1e-20;

% pool the covariance

w_proj=(cov1+cov2)\(mu2-mu1);

% project both sets of points

points1=DATAMAT1*w_proj;
points2=DATAMAT2*w_proj;

% between class separation over within class separation

separation=((w_proj'*(mu1-mu2))^2)/(w_proj'*(cov1+cov2)*w_proj);

% define mean and covariance for samples
% plot the histogram

XI=linspace(min([points1;points2]),max([points1;points2]),...
	sqrt(length(points1)+length(points2)));

[DENSITY1]=ksdensity(points1,XI);
[DENSITY2]=ksdensity(points2,XI);
