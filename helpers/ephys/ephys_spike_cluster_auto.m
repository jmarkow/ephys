function [WINDOWS TIMES TRIALS ISI STATS OUTLIERS SPIKEDATA MODEL]=...
		ephys_spike_cluster_auto(SPIKES,NOISEDATA,varargin)
%automated spike clustering using a GMM with EM
%
%
%



% spikewindows', rows x samples, each row is a windowed spike waveform

if nargin<2
	error('ephysPipeline:suavis:notenoughparams','Need 2 arguments to continue, see documentation');
end

SPIKEDATA=[];
CLUSTERPOINTS=[];
LABELS=[];
TRIALS=[];
ISI=[];
WINDOWS=[];
OUTLIERS=[];
MODEL=[];

fs=25e3;
interpolated_fs=200e3;
interpolate_fs=200e3;
proc_fs=25e3;
maxnoisetraces=1e6;
cluststart=10;
pcs=4;
workers=1;
garbage=1;
smem=1;
modelselection='icl';
use_spiketime=0; % use spiketime as a clustering feature (usually helps if SNR is low)
nparams=length(varargin);
align='min';
regularize=.01;

%outlier_cutoff=.05; % posterior probability cutoff for outliers (.6-.8 work well) [0-1, high=more aggresive]

nfeatures=10; % number of features to use, ranked by dimreduction technique

if mod(nparams,2)>0
	error('ephysPipeline:argChk','Parameters must be specified as parameter/value pairs!');
end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'fs'
			fs=varargin{i+1};
		case 'interpolated_fs'
			interpolated_fs=varargin{i+1};
		case 'proc_fs'
			proc_fs=varargin{i+1};
		case 'use_spiketime'
			use_spiketime=varaargin{i+1};
		case 'spikecut'
			spikecut=varargin{i+1};
		case 'merge'
			merge=varargin{i+1};
		case 'maxnoisetraces'
			maxnoisetaces=varargin{i+1};
		case 'ranklimit'
			ranklimit=varargin{i+1};
		case 'cluststart'
			cluststart=varargin{i+1};
		case 'pcs'
			pcs=varargin{i+1};
		case 'garbage'
			garbage=varargin{i+1};
		case 'smem'
			smem=varargin{i+1};
		case 'workers'
			workers=varargin{i+1};
		case 'modelselection'
			modelselection=varargin{i+1};
		case 'align'
			align=varargin{i+1};
	end
end

% TODO: chi2 test for dealing with outliers (simply take residuals and use chi2 CDF)
% need to deal with cell input (multiple trials), convert input to big matrix
% and spit out trial number

spikewindows=[];
spiketimes=[];
spikeifr=[];
trialnum=[];
spike_data=[];

% use ifr as a clustering feature

% string the channels together for clustering

disp('Whitening spikes');

NEWSPIKES=SPIKES;

% get the covariance matrices for whitening

[nsamples,ntrials,nchannels]=size(SPIKES(1).windows);
for i=1:length(NOISEDATA)

	noisetrials=floor(length(NOISEDATA{i})/nsamples);

	if noisetrials>maxnoisetraces
		noisetrials=maxnoisetraces;
	end

	disp(['Noise trials ' num2str(noisetrials)]);

	noisematrix=zeros(noisetrials,nsamples);

	counter=0;
	for j=1:noisetrials
		noisematrix(j,:)=NOISEDATA{i}(counter+1:counter+nsamples);
		counter=counter+nsamples;
	end

	noisematrix=noisematrix+regularize.*randn(size(noisematrix));
	noisecov=cov(noisematrix);

	R=chol(noisecov);
	invR{i}=inv(R);
end

for j=1:length(NEWSPIKES)

	[samples,trials,nchannels]=size(NEWSPIKES(j).windows);

	NEWSPIKES(j).oldwindows=NEWSPIKES(j).windows;
	for k=1:nchannels
		NEWSPIKES(j).windows(:,:,k)=[NEWSPIKES(j).windows(:,:,k)'*invR{k}]';
	end


	%figure(1);
	%for k=1:size(NEWSPIKES(j).windows,3)

	%	subplot(size(NEWSPIKES(j).windows,3),1,k);plot(NEWSPIKES(j).windows(:,:,k));
	%	hold on;

	%end

	%pause();

	% upsample and realign

	UPNEWSPIKES(j)=ephys_spike_upsample_align(NEWSPIKES(j),'interpolate_fs',interpolate_fs,'align',align);

	%figure(2);
	%for k=1:size(UPNEWSPIKES(j).windows,3)

	%	subplot(size(UPNEWSPIKES(j).windows,3),1,k);plot(UPNEWSPIKES(j).windows(:,:,k));
	%	hold on;

	%end

	%pause();


	[samples,trials,nchannels]=size(UPNEWSPIKES(j).windows);

	tmp=[];

	for k=1:nchannels
		tmp=[tmp;UPNEWSPIKES(j).windows(:,:,k)];
	end

	clusterspikewindowscell{j}=tmp;
	
	tmp=[];

	for k=1:nchannels
		tmp=[tmp;UPNEWSPIKES(j).oldwindows(:,:,k)];
	end

	storespikewindowscell{j}=tmp;
	trialscell{j}=ones(trials,1).*j;

end

clusterspikewindows=cat(2,clusterspikewindowscell{:});
storespikewindows=cat(2,storespikewindowscell{:});
trialnum=cat(1,trialscell{:});
spiketimes=cat(2,UPNEWSPIKES(:).times);

clearvars NEWSPIKES SPIKES UPNEWSPIKES clusterspikewindowscell storespikewindowscell;

[idx spikedata MODEL]=ephys_spike_gmmsort(clusterspikewindows,...
	'proc_fs',proc_fs,'fs',fs,'interpolated_fs',interpolated_fs,...
	'smem',smem,'garbage',garbage,'maxnoisetraces',maxnoisetraces,...
	'cluststart',cluststart,'pcs',pcs,'workers',workers,'modelselection',...
	modelselection);

features=size(spikedata,2); % what's the dimensionality of the data used for sorting?
clusters=unique(idx(idx>0)); % how many clusters?
nclust=length(clusters);

% number of spikes per cluster is simply the number of labels

nspikes=[];
for i=1:nclust
	nspikes(i)=sum(idx==clusters(i));
end

[val loc]=sort(nspikes,'descend');

% make the number contiguous and sort by number of spikes, descending

LABELS=zeros(size(idx));

for i=1:length(clusters)
	LABELS(idx==clusters(loc(i)))=i;	
end


%TRIALS=trialnum;

clusters=unique(LABELS(LABELS>0));

% now assess the cluster quality ,
% take each cluster and check the FP and FN rate

% use l ratio of .05

for i=1:nclust

	clusterlocs=find(LABELS==i);
	otherlocs=find((LABELS~=i)&(LABELS>0));

	% get the feature data

	clusterpoints=spikedata(clusterlocs,:);
	otherpoints=spikedata(otherlocs,:);

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
		STATS.lratio(i)=sum(1-chi2cdf(mahaldist,features))/nclustpoints;
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

[uniq_trial trial_boundary trial_group]=unique(trialnum);
trial_boundary=[0;trial_boundary];

OUTLIERS=storespikewindows(:,isnan(idx));

for i=1:nclust

	WINDOWS{i}=storespikewindows(:,LABELS==i);
	TIMES{i}=spiketimes(LABELS==i);
	SPIKEDATA{i}=spikedata(LABELS==i,:);

	spikeisitmp=[];
	trialtmp=[];

	for j=1:length(uniq_trial)
		
		% all spike times in this trial	

		currtrial=spiketimes(trial_boundary(j)+1:trial_boundary(j+1));

		% now all spike ids from this trial

		currlabels=LABELS(trialnum==uniq_trial(j));

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
	SPIKEDATA{end+1}=spikedata(isnan(idx),:);
end
% if we have any outliers, assign them to the final cluster



