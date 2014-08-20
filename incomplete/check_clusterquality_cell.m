function [STATS]=check_clusterquality_cell(CLUSTER)
%
%
%
%

SPIKEDATA=CLUSTER.spikedata;
MODEL=CLUSTER.model;

nmodels=size(MODEL.mu,1);
nclust=nmodels;

STATS.lratio=zeros(1,nclust);
STATS.isod=zeros(1,nclust);
nfeatures=size(SPIKEDATA{1},2);
npoints=1e5; % point from each cluster

for i=1:nclust

	% get the feature data

	%clusterpoints=SPIKEDATA(clusterlocs,:);
	
	clusterpoints=SPIKEDATA{i};
	otherpoints=cat(1,SPIKEDATA{setdiff(1:nclust,i)});

	% l ratio is the sum inv chi2cdf of mahal distance of all other points over
	% n spikes

	nclustpoints=size(clusterpoints,1);

	if size(otherpoints,1)<=size(otherpoints,2) || size(clusterpoints,1)<=size(clusterpoints,2)
		warning('ephysPipeline:spikesort:toofewrowsmahal',...
			'Too few rows for Mahal distance calculation');
		STATS.lratio(i)=NaN;
		STATS.isod(i)=NaN;
		continue;
	end

	% regularize the covariance matrix to make sure it's invertible

	covm=MODEL.sigma(:,:,i)+eye(size(MODEL.sigma(:,:,i)))*1e-5;
	invcovm=inv(covm);	
	refmean=MODEL.mu(i,:);

	mahaldist=zeros(1,size(otherpoints,1));

	% take the square mahal distance

	%for j=1:size(otherpoints,1)
		%mahaldist(j)=(otherpoints(j,:)-refmean)*invcovm*(otherpoints(j,:)-refmean)';
	%end
	
	mahaldist=mahal(otherpoints,clusterpoints);

	% l-ratio

	if length(otherpoints>0)
		STATS.lratio(i)=sum(1-chi2cdf(mahaldist,nfeatures))/nclustpoints;
	else
		STATS.lratio(i)=NaN;
	end

	% iso distance

	if length(mahaldist)>=nclustpoints
		sortmahal=sort(mahaldist,'ascend');
		STATS.isod(i)=sortmahal(nclustpoints);
	else
		STATS.isod(i)=NaN;
	end

end

% or generate points from each cluster
%

if any(MODEL.garbage)
	garbage=1;
	nclust=nclust+1;
	alldata=cat(1,SPIKEDATA{:});
	datarange=range(alldata);
	datamin=min(alldata);
	datamax=max(alldata);
	p=prod(datarange);
	P=(1/p);
end

total_points=nclust*npoints;
labels=zeros(total_points,1);
D=size(SPIKEDATA{1},2);
points=zeros(total_points,D);

for i=1:nclust
	startpoint=((i-1)*npoints);
	labels(startpoint+1:startpoint+npoints)=i;

	% generate using uniform density spanning entire space

	if i==nclust & garbage
		for j=1:D
			points(startpoint+1:startpoint+npoints,j)=(rand(npoints,1).*datamax(j))-datamin(j);
		end
	else
		points(startpoint+1:startpoint+npoints,:)=mvnrnd(MODEL.mu(i,:),MODEL.sigma(:,:,i),npoints);
	end
end

for i=1:nclust

	if i==nclust & garbage
		prob(:,i)=P.*ones(size(points,1),1);
	else
		prob(:,i)=mvnpdf(points,MODEL.mu(i,:),MODEL.sigma(:,:,i));
	end
end

prob=prob.*repmat(MODEL.mixing,[total_points 1]);
[~,classification]=max(prob,[],2);

for i=1:size(MODEL.mu,1)

	% how many false positives and false negatives
	%
	fp=sum((classification==i)&(labels~=i));

	% false negatives
	%
	fn=sum((classification~=i)&(labels==i));

	STATS.misclassification(i)=(fp/npoints)+(fn/npoints);
end
