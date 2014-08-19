function STATS=check_clusterquality_cell(CLUSTER)
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
