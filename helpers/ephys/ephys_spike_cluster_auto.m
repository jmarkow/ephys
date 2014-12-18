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

nparams=length(varargin);

SPIKEDATA=[];
CLUSTERPOINTS=[];
LABELS=[];
TRIALS=[];
ISI=[];
WINDOWS=[];
OUTLIERS=[];
MODEL=[];

fs=25e3;
%interpolated_fs=200e3;
interpolate_fs=200e3;
proc_fs=25e3;
maxnoisetraces=1e6;
cluststart=10;
pcs=4;
workers=1;
garbage=1;
smem=1;
modelselection='icl';
align_method='min';
regularize=.01;
noisewhiten=1;

%outlier_cutoff=.05; % posterior probability cutoff for outliers (.6-.8 work well) [0-1, high=more aggresive]

nfeatures=10; % number of features to use, ranked by dimreduction technique

if mod(nparams,2)>0
	error('ephysPipeline:argChk','Parameters must be specified as parameter/value pairs!');
end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'fs'
			fs=varargin{i+1};
		case 'interpolate_fs'
			interpolate_fs=varargin{i+1};
		case 'proc_fs'
			proc_fs=varargin{i+1};
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
		case 'align_method'
			align_method=varargin{i+1};
		case 'noisewhiten'
			noisewhiten=varargin{i+1};
	end
end

% TODO: chi2 test for dealing with outliers (simply take residuals and use chi2 CDF)
% need to deal with cell input (multiple trials), convert input to big matrix
% and spit out trial number

downfact=interpolate_fs/fs;

if mod(downfact,1)~=0
	error('Downsample rate must be an integer');
end

spikewindows=[];
spiketimes=[];
spikeifr=[];
trialnum=[];
spike_data=[];

% use ifr as a clustering feature

% string the channels together for clustering
% get the covariance matrices for whitening

[nsamples,ntrials,nchannels]=size(SPIKES(1).windows);

if noisewhiten
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
end

for j=1:length(SPIKES)

	[samples,trials,nchannels]=size(SPIKES(j).windows);

	% store unwhitened times and use the unwhitened spikes for spike times

	SPIKES(j).storewindows=SPIKES(j).windows;

	% comment out the next three lines to not noise-whiten

	if noisewhiten
		for k=1:nchannels
			SPIKES(j).windows(:,:,k)=[SPIKES(j).windows(:,:,k)'*invR{k}]';
		end
	end

	% upsample and align, then downsample and whiten!!!

	% masking

	%spikemask=ones(size(SPIKES(j).windows));
	%spikemask([1:15 end-15:end],:,:)=0;
	%SPIKES(j).windows=SPIKES(j).windows.*spikemask;

	alignspikes=ephys_spike_upsample_align(SPIKES(j),'interpolate_fs',interpolate_fs,'align_method',align_method);	
	CLUSTSPIKES(j)=alignspikes;

	% cluster with the decimated spikes

	[~,trials,nchannels]=size(CLUSTSPIKES(j).windows);

	tmp=[];
	
	for k=1:nchannels
		tmp=[tmp;CLUSTSPIKES(j).windows(:,:,k)];
	end

	clusterspikewindowscell{j}=tmp;
	
	tmp=[];

	for k=1:nchannels
		tmp=[tmp;CLUSTSPIKES(j).storewindows(:,:,k)];
	end

	storespikewindowscell{j}=tmp;
	trialscell{j}=ones(trials,1).*j;

end

[nsamples,ntrials,nchannels]=size(CLUSTSPIKES(1).windows);

clusterspikewindows=cat(2,clusterspikewindowscell{:});
storespikewindows=cat(2,storespikewindowscell{:});
trialnum=cat(1,trialscell{:});
spiketimes=cat(2,CLUSTSPIKES(:).storetimes);

clearvars SPIKES CLUSTSPIKES clusterspikewindowscell storespikewindowscell;

[idx spikedata MODEL]=ephys_spike_gmmsort(clusterspikewindows,...
	'proc_fs',proc_fs,'fs',fs,'interpolate_fs',interpolate_fs,...
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
OUTLIERS=storespikewindows(:,LABELS==0);

% now assess the cluster quality ,
% take each cluster and check the FP and FN rate

% use l ratio of .05

[WINDOWS TIMES TRIALS SPIKEDATA ISI STATS]=...
	check_clusterquality(storespikewindows,spiketimes,spikedata,LABELS,trialnum,MODEL);
