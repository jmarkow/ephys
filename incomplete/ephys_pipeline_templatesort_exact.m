function [LABELS,clust]=ephys_pipeline_templatesort_exact(SPIKE_DATA,NOISE_DATA,varargin)
%
%
%
%
%
%

% calculate merge and accept thresholds based on estimate of the noise
% match the noise sampling rate to the spike sampling rate and 

alpha=.0001;
interpolate_fs=200e3;
fs=25e3;
proc_fs=25e3; % downsample the spikes by a factor of 4 and upsample the noise by
spikecut=50;
maxnoisetraces=100e3;

nparams=length(varargin);

if mod(nparams,2)>0
	error('ephysPipeline:argChk','Parameters must be specified as parameter/value pairs!');
end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'alpha'
			alpha=varargin{i+1};
		case 'interpolate_fs'
			interpolate_fs=varargin{i+1};
		case 'fs'
			fs=varargin{i+1};
		case 'proc_fs'
			proc_fs=varargin{i+1};
		case 'maxnoisetraces'
			maxnoisetraces=varargin{i+1};
	end
end

downfact=interpolate_fs/proc_fs;

if mod(downfact,1)~=0
	error('ephyspipeline:templatesortexact:baddownfact',...
		'Need to downsample by an integer factor');
end

SPIKE_DATA=downsample(SPIKE_DATA,downfact);

% consider setting higher depending on results

[nsamples,ntrials]=size(SPIKE_DATA);

% resample the noise data to match the spike FS?

noisetrials=floor(length(NOISE_DATA)/nsamples);

if noisetrials>maxnoisetraces
	noisetrials=maxnoisetraces;
end

noisematrix=zeros(noisetrials,nsamples);

counter=0;
for i=1:noisetrials
	noisematrix(i,:)=NOISE_DATA(counter+1:counter+nsamples);
	counter=counter+nsamples;
end

% get the covariance matrix of the noise

noisecov=cov(noisematrix);
meancov=mean(noisecov(:));
noisecov=noisecov+(1e-6*meancov).*randn(size(noisecov));

% cholesky factorization

R=chol(noisecov);

% inverse gives us the whitening matrix

invR=inv(R);

SPIKE_DATA=SPIKE_DATA'*invR;
invC=inv(cov(SPIKE_DATA));
SPIKE_DATA=SPIKE_DATA';

% threshold is the distance that include n% of points in the chi2 distribution

thresh_t=chi2inv(1-alpha,nsamples);

% number of standard deviations for merging clusters

thresh_m=3;

% initialize relevant variables...

LABELS=zeros(ntrials,1);

% create the first cluster

clust{1}=[SPIKE_DATA(:,1)];
clust_template(:,1)=mean(clust{1},2);
LABELS(1)=1;

nclust=1;
repair=0;

for i=2:ntrials

	i
	% calculate the residuals

	%spikeresidual=zeros(1,length(clust));

	spikeresidual=[];
	
	nclust=length(clust)

	spikeresidual=sum((repmat(SPIKE_DATA(:,i),[1 nclust])-clust_template).^2);
	% get the min

	addclust=1;
	[val loc]=min(spikeresidual);

	if val<=thresh_t

		% accept the spike waveform into the cluster

		clust{loc}=[clust{loc} SPIKE_DATA(:,i)];
		clust_template(:,loc)=mean(clust{loc},2);
		LABELS(i)=loc;
		addclust=0;

	end

	% merge any clusters?

	if addclust
		clust{end+1}=SPIKE_DATA(:,i);
		clust_template=[clust_template mean(clust{end},2)];
		LABELS(i)=length(clust);
		nclust=nclust+1;
		repair=1;
	end
	
	count=1;

	while count>0 && nclust>1

		% nchoosek is a time-suck, only recompute the pairs of we need them

		if repair
			clustpairs=nchoosek(1:nclust,2);
			repair=0;
		end
		
		clustdist=zeros(size(clustpairs,1),1);
		count=0;

		for j=1:size(clustpairs,1)

			c1=clustpairs(j,1);
			c2=clustpairs(j,2);

			m1=clust_template(:,c1);
			m2=clust_template(:,c2);	

			%dist=(m1-m2)'*invC*(m1-m2);
			dist=sum((m1-m2).^2);

			clustdist(j)=sqrt(dist);

		end

		count=sum(clustdist<=thresh_m);

		if count>0

			[val loc]=min(clustdist);

			merge1=clustpairs(loc(1),1);
			merge2=clustpairs(loc(1),2);

			fprintf('Merging clusters %g and %g\n',merge1,merge2);

			LABELS(LABELS==merge2)=merge1;
			clust{merge1}=[clust{merge1} clust{merge2}];
			clust_template(:,merge1)=mean(clust_template(:,[merge1 merge2]),2);
			clust(merge2)=[];
			clust_template(:,merge2)=[];
			nclust=nclust-1;
			repair=1;

			% ensure the labeling is contiguous again
		end
	end

	%to_del=[];
	%for j=1:length(clust)
	%	if isempty(clust{j})
	%		to_del=[to_del j];
	%	end
	%end

	%clust(to_del)=[];

end

to_del=[];

for i=1:length(clust)
	nspikes=size(clust{i},2);

	if nspikes<=spikecut
		to_del=[to_del i];
	end
end

clust(to_del)=[];

for i=1:length(to_del)
	LABELS(LABELS==to_del(i))=NaN;
end

