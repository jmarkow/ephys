function [LABELS SPIKE_DATA]=ephys_pipeline_templatesort_exact(SPIKE_DATA,NOISE_DATA,varargin)
%
%
%
%
%
%

% calculate merge and accept thresholds based on estimate of the noise
% match the noise sampling rate to the spike sampling rate and 

interpolate_fs=200e3;
fs=25e3;
proc_fs=25e3; % downsample the spikes by a factor of 4 and upsample the noise by
spikecut=50;
correction=2.15;
mergecorrection=2.15;
periodcheck=3e3; % should be large, e.g. 10e3
periodcut=2; % should be small, e.g. 5-10
weights=[.0002 .0002]; % set the edges to zero to avoid edge effects
maxnoisetraces=1e6;
method='exact';
exact_alpha=.999;
exact_merge=4;
mean_spikes=100;

nparams=length(varargin);

if mod(nparams,2)>0
	error('ephysPipeline:argChk','Parameters must be specified as parameter/value pairs!');
end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'maxclust'
			maxclust=varargin{i+1};
		case 'interpolate_fs'
			interpolate_fs=varargin{i+1};
		case 'fs'
			fs=varargin{i+1};
		case 'proc_fs'
			proc_fs=varargin{i+1};
		case 'spikecut'
			spikecut=varargin{i+1};
		case 'correction'
			correction=varargin{i+1};
		case 'mergecorrection'
			mergecorrection=varargin{i+1};
		case 'periodcheck'
			periodcheck=varargin{i+1};
		case 'periodcut'
			periodcut=varargin{i+1};
		case 'weights'
			weights=varargin{i+1};
		case 'method'
			method=varargin{i+1};
		case 'exact_alpha'
			exact_alpha=varargin{i+1};
		case 'exact_merge'
			exact_merge=varargin{i+1};
		case 'maxnoisetraces'
			maxnoisetaces=varargin{i+1};
	end
end

if isempty(mean_spikes)
	mean_spikes=inf;
end

% consider setting higher depending on results

[nsamples,ntrials]=size(SPIKE_DATA);

% resample the noise data to match the spike FS?
% get the covariance matrix of the noise

disp(['Sort fs ' num2str(proc_fs)]);
disp(['Period cut ' num2str(periodcut)]);
disp(['Period check ' num2str(periodcheck)]);
disp(['Spike cut ' num2str(spikecut)]);

switch lower(method(1))

	case 'e'

		% downsample spikes back to original FS

		downfact=interpolate_fs/fs;

		if mod(downfact,1)~=0
			error('ephyspipeline:templatesortexact:baddownfact',...
				'Need to downsample by an integer factor');
		end
	
		SPIKE_DATA=downsample(SPIKE_DATA,downfact);

		[nsamples,ntrials]=size(SPIKE_DATA);

		disp('Exact method');
		disp(['Chi2 alpha ' num2str(exact_alpha)]);
		disp(['Merge cutoff (std) ' num2str(exact_merge)]);
		disp(['Max noise traces ' num2str(maxnoisetraces)]);

		% do we need to resample the noise data to match the sorting fs?

		noisetrials=floor(length(NOISE_DATA)/nsamples);

		if noisetrials>maxnoisetraces
			noisetrials=maxnoisetraces;
		end

		disp(['Noise trials ' num2str(noisetrials)]);
		
		noisematrix=zeros(noisetrials,nsamples);

		counter=0;
		for i=1:noisetrials
			noisematrix(i,:)=NOISE_DATA(counter+1:counter+nsamples);
			counter=counter+nsamples;
		end

		noisematrix=noisematrix+.01.*randn(size(noisematrix));
		noisecov=cov(noisematrix);

		R=chol(noisecov);
		invR=inv(R);

		% UPSAMPLE AGAIN USING SPLINES IF PROC_FS>FS

		resample_factor=proc_fs/fs;

		if mod(resample_factor,1)~=0 | resample_factor<1
			error('ephyspipeline:templatesort:resamplenotinteger');
		end

		OLD_DATA=SPIKE_DATA;

		[coef score]=princomp(SPIKE_DATA');

		SPIKE_DATA=SPIKE_DATA'*invR;	
		SPIKE_DATA=SPIKE_DATA';	

		if resample_factor>1
			oldtimepoints=1:nsamples;
			newtimepoints=linspace(1,nsamples,nsamples*resample_factor);
			newspikedata=zeros(nsamples*resample_factor,ntrials);

			disp(['Upsampling ' num2str(resample_factor) ' fold']);

			for i=1:ntrials
				newspikedata(:,i)=spline(oldtimepoints,SPIKE_DATA(:,i),newtimepoints);
			end

			SPIKE_DATA=newspikedata;
			clear newspikedata;
			nsamples=nsamples*resample_factor;
		end

		thresh_t=chi2inv(exact_alpha,nsamples);

		disp(['Chi2 cutoff ' num2str(thresh_t)]);

		noisematrix=noisematrix*invR;	
		thresh_m=exact_merge;

	case 'a'


		downfact=interpolate_fs/proc_fs;

		if mod(downfact,1)~=0
			error('ephyspipeline:templatesortexact:baddownfact',...
				'Need to downsample by an integer factor');
		end

		SPIKE_DATA=downsample(SPIKE_DATA,downfact);

		[nsamples,ntrials]=size(SPIKE_DATA);

		disp('Approximate method');
		disp(['Threshold correction ' num2str(correction)]);
		disp(['Merge correction ' num2str(mergecorrection)]);

		noise_est=std(NOISE_DATA);

		disp(['Noise estimate ' num2str(noise_est)]);

		thresh=nsamples*(noise_est^2)
		thresh_t=thresh*correction;
		thresh_m=thresh*mergecorrection;

	otherwise

end

if ~isempty(weights) & length(weights)==2

	disp(['Weights (zero values from edges) ' num2str(weights)]);

	weight_samples=ceil(weights.*proc_fs);
	weights=ones(nsamples,1);
	weights(1:weight_samples(1))=0;
	weights(nsamples-weight_samples(2):nsamples)=0;

	SPIKE_DATA=SPIKE_DATA.*repmat(weights,[1 ntrials]);
	nsamples=sum(weights>0);

	if lower(method(1))=='a'

		thresh=sum(weights.*(noise_est^2))
		thresh_t=thresh*correction;
		thresh_m=thresh*mergecorrection;

	else

		thresh_t=chi2inv(exact_alpha,nsamples)

	end
end

% threshold is the distance that include n% of points in the chi2 distribution

LABELS=zeros(ntrials,1);

% create the first cluster

clust_template(:,1)=SPIKE_DATA(:,1);
LABELS(1)=1;
nclust=1;

repair=0;
[nblanks formatstring]=progressbar(100);

disp('Pre-computing combinations');

load('combos.mat','combos');
maxclust=length(combos);

fprintf(1,['Progress:  ' blanks(nblanks)]);

for i=2:ntrials

	% calculate the residuals

	fprintf(1,formatstring,round((i/ntrials)*100));

	count=0;

	% any noise clusters to delete?

	if ~isempty(periodcheck) & mod(i,periodcheck)==0

		nclust=size(clust_template,2);

		to_del=[];
		for j=1:nclust
			if sum(LABELS==j)<periodcut
				to_del=[to_del j];
			end
		end

		clust_template(:,to_del)=[];
		grps=unique(LABELS(LABELS>0));

		for j=1:length(to_del)
			LABELS(LABELS==grps(to_del(j)))=NaN;
		end

		grps=unique(LABELS(LABELS>0));
		idx=LABELS;

		for j=1:length(grps)
			idx(LABELS==grps(j))=j;
		end

		LABELS=idx;

	end

	% collect the residuals

	spikeresidual=[];
	
	nclust=size(clust_template,2);
	spikeresidual=sum((repmat(SPIKE_DATA(:,i),[1 nclust])-clust_template).^2);

	[val loc]=min(spikeresidual);

	if val<=thresh_t

		% accept the spike waveform into the cluster
		% just update with the current waveform!

		LABELS(i)=loc;
		hitidx=find(LABELS==loc);

		hitcut=min([length(hitidx) mean_spikes]);
		hitidx=hitidx(1:hitcut);

		clust_template(:,loc)=mean(SPIKE_DATA(:,hitidx),2);

		newdist=zeros(1,nclust-1);
		checktemplates=setdiff(1:nclust,loc);

		for j=1:nclust-1
			newdist(j)=sum((clust_template(:,checktemplates(j))-clust_template(:,loc)).^2);
		end


	else

		clust_template=[clust_template SPIKE_DATA(:,i)];
		LABELS(i)=nclust+1;
		nclust=nclust+1;
		loc=nclust;
		newdist=[];

	end

	% get distance between update waveform and other mean waveforms

	if nclust>maxclust
		warning('ephyspipeline:templatesort:maxclustexceeded','Maximum number of clusters exceeded');
		LABELS(1:end)=NaN;
		return	
	end

end

count=1;
nclust=size(clust_template,2);

for i=1:nclust
	clust_template(:,i)=mean(SPIKE_DATA(:,LABELS==i),2);
end

while count>0 && nclust>1 

	% nchoosek is a time-suck

	clustpairs=combos{nclust};
	clustdist=zeros(size(clustpairs,1),1);
	count=0;

	% get the squared euclidean distance

	for j=1:size(clustpairs,1)

		c1=clustpairs(j,1);
		c2=clustpairs(j,2);

		switch lower(method(1))
			case 'e'
				dist=norm(clust_template(:,c1)-clust_template(:,c2));
			case 'a'
				dist=sum((clust_template(:,c1)-clust_template(:,c2)).^2);
			otherwise
		end

		clustdist(j)=dist;

	end

	count=sum(clustdist<thresh_m);

	if count>0

		[val loc]=min(clustdist);

		merge1=clustpairs(loc(1),1);
		merge2=clustpairs(loc(1),2);

		%fprintf('Merging clusters %g and %g\n',merge1,merge2);

		LABELS(LABELS==merge2)=merge1;

		grps=unique(LABELS(LABELS>0));

		idx=LABELS;

		for j=1:length(grps)
			idx(LABELS==grps(j))=j;
		end

		LABELS=idx;

		clust_template(:,merge1)=mean(SPIKE_DATA(:,[(LABELS==merge1)]),2);
		clust_template(:,merge2)=[];

		nclust=nclust-1;

		% ensure the labeling is contiguous again
	end


end

% one final round, assign spikes to template or discard as noise

count=1;
nclust=size(clust_template,2);
fprintf(1,'\n');

% remove noise clusters

grps=unique(LABELS(~isnan(LABELS)));
to_del=[];

for i=1:nclust

	nspikes=sum(LABELS==i);

	if nspikes<=spikecut
		to_del=[to_del i];
	end
end

for i=1:length(to_del)
	LABELS(LABELS==grps(to_del(i)))=NaN;
end

clust_template(:,to_del)=[];
nclust=size(clust_template,2);

for i=1:ntrials

	if LABELS(i)==NaN,
		continue;
	end

	spikeresidual=sum((repmat(SPIKE_DATA(:,i),[1 nclust])-clust_template).^2);

	if all(spikeresidual>thresh_t)
		LABELS(i)=NaN;
		continue;
	end

	[val loc]=min(spikeresidual);
	LABELS(i)=loc;

end

disp(['N Clust:  ' num2str(length(unique(LABELS(LABELS>0))))]);
