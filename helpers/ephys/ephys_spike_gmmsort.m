function [LABELS SPIKE_DATA MODEL]=ephys_pipeline_gmmsort(SPIKE_DATA,varargin)
%
%
%
%
%
%

% calculate merge and accept thresholds based on estimate of the noise
% match the noise sampling rate to the spike sampling rate and 

nparams=length(varargin);

interpolate_fs=200e3;
fs=25e3;
proc_fs=25e3; % downsample the spikes by a factor of 4 and upsample the noise by
mergecut=3;
maxnoisetraces=1e6;
cluststart=10;
pcs=2;
pcareplicates=5;
clustreplicates=1;
garbage=1;
workers=1;
modelselection='icl';
smem=1;
outliercut=.2; % exclude outliers from robpca

if mod(nparams,2)>0
	error('ephysPipeline:argChk','Parameters must be specified as parameter/value pairs!');
end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'cluststart'
			cluststart=varargin{i+1};
		case 'interpolate_fs'
			interpolate_fs=varargin{i+1};
		case 'fs'
			fs=varargin{i+1};
		case 'proc_fs'
			proc_fs=varargin{i+1};
		case 'pcareplicates'
			pcareplicates=varargin{i+1};
		case 'clustreplicates'
			clustreplicates=varargin{i+1};
		case 'pcs'
			pcs=varargin{i+1};
		case 'garbage'
			garbage=varargin{i+1};
		case 'workers'
			workers=varargin{i+1};
		case 'modelselection'
			modelselection=varargin{i+1};
		case 'smem'
			smem=varargin{i+1};
		case 'outliercut'
			outliercut=varargin{i+1};
	end
end

% consider setting higher depending on results

[nsamples,ntrials]=size(SPIKE_DATA);

% resample the noise data to match the spike FS?
% get the covariance matrix of the noise

disp(['Sort fs ' num2str(proc_fs)]);
disp(['Starting clusters ' num2str(cluststart)]);
disp(['PCS:  ' num2str(pcs)]);
disp(['Garbage collection: ' num2str(garbage)]);
disp(['SMEM:  ' num2str(smem)]);
disp(['Workers (deployed only):  ' num2str(workers)]);
disp(['Model selection ' modelselection]);

% downsample spikes back to original FS

downfact=interpolate_fs/fs;

if mod(downfact,1)~=0
	error('ephyspipeline:templatesortexact:baddownfact',...
		'Need to downsample by an integer factor');
end

SPIKE_DATA=downsample(SPIKE_DATA,downfact);

[nsamples,ntrials]=size(SPIKE_DATA);

% do we need to resample the noise data to match the sorting fs?

if isdeployed
	matlabpool('open',workers);
end

likelihood=zeros(1,pcareplicates);
parfor i=1:pcareplicates
	tmpnewmodel{i}=gmem(SPIKE_DATA',[],1,'garbage',1,'merge',0,'debug',0);
	likelihood(i)=tmpnewmodel{i}.likelihood;
end

% choose the model with the highest likelihood

[~,loc]=max(likelihood);
newmodel=tmpnewmodel{loc(1)};
[v,d]=eigs(newmodel.sigma(:,:,1));
newscore=-SPIKE_DATA'*v;
rankcut=pcs;

% these comprise the outliers before the projection...

outlierpoints=newmodel.R(:,2)>=outliercut;

disp(['Robust PCA outliers:  ' num2str(sum(outlierpoints))]);

SPIKE_DATA=newscore(:,1:rankcut);

% only cluster the non-outliers

newscore=newscore(~outlierpoints,:);

% might consider removing outliers before we project to the new space (Sahani,1999)

% commented out, but need to consider later...
% outlierpoints=find(newmodel.R(:,2)>.98);

% outliers are any points with garbage components R>.9

oldstate1=warning('off','stats:gmdistribution:FailedToConverge');
oldstate2=warning('off','stats:kmeans:FailedToConvergeRep');
oldstate3=warning('off','stats:kmeans:FailedToConverge');

parfor i=1:1:length(cluststart)

	tmpclustobj={};
	startmu=[];
	startcov=[];
	mixing=[];

	startobj=struct('mu',startmu,'sigma',startcov,'mixing',mixing);

	loglikelihood=zeros(1,clustreplicates);
	idx=kmeans(newscore(:,1:rankcut),cluststart(i),'replicates',5);

	%% set up initial model

	mu=[];
	for j=1:cluststart(i)
		startmu(j,:)=mean(newscore(idx==j,1:rankcut))';
		startcov(:,:,j)=diag(var(newscore(:,1:rankcut)));
	end

	startobj.mu=startmu;
	startobj.sigma=startcov;

	for j=1:cluststart(i)
		startobj.mixing(j)=sum(idx==j)/length(idx);
	end

	
	for j=1:clustreplicates
		tmpclustobj{j}=gmem(newscore(:,1:rankcut),startobj,cluststart(i),...
			'garbage',garbage,'merge',smem,'debug',0);
		loglikelihood(j)=tmpclustobj{j}.likelihood;
	end

	% only keep the clustobj with the best likelihood

	[~,loc]=max(loglikelihood);
	clustobj{i}=tmpclustobj{loc(1)};
	BIC(i)=clustobj{i}.BIC;
	MML(i)=clustobj{i}.MML;
	ICL(i)=clustobj{i}.ICL;

end

if isdeployed
	matlabpool('close');
end

warning(oldstate1);
warning(oldstate2);
warning(oldstate3);

switch lower(modelselection(1))
	case 'b'
		[~,loc]=min(BIC);
	case 'm'
		[~,loc]=max(MML);
	case 'i'
		[~,loc]=min(ICL);
	otherwise
end

clustermodel=clustobj{loc(1)};
MODEL=clustermodel;

% get the labels from the responsibilities (probability of each cluster given each datapoint)

idx=[];
for i=1:size(clustermodel.R,1)
	posteriors=clustermodel.R;
	[~,idx(i)]=max(posteriors(i,:));
end

if garbage
	garbageidx=find(clustermodel.garbage);
	idx(idx==garbageidx)=NaN;
end

LABELS=zeros(ntrials,1);

% what did we label through clustering

LABELS(~outlierpoints)=idx;

% pre-pca outliers

LABELS(outlierpoints)=NaN;

grps=unique(LABELS(LABELS>0));
nclust=length(grps);

% ensure the labeling is contiguous

idx=LABELS;
for i=1:nclust
	idx(LABELS==grps(i))=i;
end

LABELS=idx;
clear idx;

disp(['N Clust:  ' num2str(length(unique(LABELS(LABELS>0))))]);
