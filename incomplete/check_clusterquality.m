function [WINDOWS TIMES TRIALS CLUSTDATA ISI STATS]=check_clusterquality(SPIKEWINS,SPIKETIMES,SPIKEDATA,LABELS,TRIALNUM,MODEL)
%
%
%
%

nmodels=size(MODEL.mu,1);

clusters=unique(LABELS(LABELS>0)); % how many clusters?
nclust=length(clusters);

STATS.lratio=zeros(1,nclust);
STATS.isod=zeros(1,nclust);
nfeatures=size(SPIKEDATA,2);

for i=1:nclust

	clusterlocs=find(LABELS==i);
	otherlocs=find((LABELS~=i)&(LABELS>0));

	% get the feature data

	clusterpoints=SPIKEDATA(clusterlocs,:);
	otherpoints=SPIKEDATA(otherlocs,:);

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

	for j=1:size(otherpoints,1)
		mahaldist(j)=(otherpoints(j,:)-refmean)*invcovm*(otherpoints(j,:)-refmean)';
	end
	
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

% resort labels by number of spikes, first vector is cluster id, and the second the 
% point where those labels end (e.g. 1 1 1 2 2 2 2 would return 1 2 and 3 7)

% return labels, and windows and ISI sorted by cluster IDX

[uniq_trial trial_boundary trial_group]=unique(TRIALNUM);

if length(uniq_trial)==1
	trial_boundary=[0;length(SPIKETIMES)];
else
	trial_boundary=[0;trial_boundary];
end

OUTLIERS=SPIKEWINS(:,LABELS==0);

size(TRIALNUM)
size(LABELS)
size(SPIKETIMES)

for i=1:nclust

	WINDOWS{i}=SPIKEWINS(:,LABELS==i);
	TIMES{i}=SPIKETIMES(LABELS==i);
	CLUSTDATA{i}=SPIKEDATA(LABELS==i,:);

	spikeisitmp=[];
	trialtmp=[];

	for j=1:length(uniq_trial)
		
		% all spike times in this trial	

		trial_idx=find(TRIALNUM==uniq_trial(j));

		%currtrial=SPIKETIMES(trial_boundary(j)+1:trial_boundary(j+1));

		currtrial=SPIKETIMES(trial_idx);

		% now all spike ids from this trial

		currlabels=LABELS(trial_idx);

		% spike times for this cluster

		currtrial=currtrial(currlabels==clusters(i));
		trialtmp=[trialtmp;uniq_trial(j)*ones(length(currtrial),1)];
		currisi=(diff(currtrial));
		spikeisitmp=[spikeisitmp;currisi(:)];
	end

	TRIALS{i}=trialtmp';
	ISI{i}=spikeisitmp'; 
end

if ~isempty(OUTLIERS)
	CLUSTDATA{end+1}=SPIKEDATA(LABELS==0,:);
end

STATS.contamination=check_contamination(SPIKEDATA,MODEL);
% if we have any outliers, assign them to the final cluster



