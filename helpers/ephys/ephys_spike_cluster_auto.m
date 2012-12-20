function [LABELS TRIALS ISI WINDOWS CLUSTERPOINTS STATS OUTLIERS]=...
		ephys_spike_cluster_auto(SPIKEWINDOWS,SPIKETIMES,varargin)
%automated spike clustering using a GMM with EM
%
%
%



% spikewindows', rows x samples, each row is a windowed spike waveform

if nargin<2
	error('ephysPipeline:suavis:notenoughparams','Need 2 arguments to continue, see documentation');
end

if ~license('test','Statistics_Toolbox')
	error('ephysPipeline:toolboxChk','Need statistics toolbox for clustering.');
end

CLUSTERPOINTS=[];
LABELS=[];
TRIALS=[];
ISI=[];
WINDOWS=[];
OUTLIERS=[];

fs=25e3;
use_spiketime=0; % use spiketime as a clustering feature (usually helps if SNR is low)
nparams=length(varargin);
maxcoeffs=10; % number of wavelet coefficients to use (sorted by KS statistic)
wavelet_method='bi';
wavelet_mpca=1;
%outlier_cutoff=.05; % posterior probability cutoff for outliers (.6-.8 work well) [0-1, high=more aggresive]

outlier_cutoff=75; % l-infty cutoff for residual test
clust_choice='bic'; % knee is more tolerant, choice BIC (b) or AIC (a) for more sensitive clustering
nfeatures=10; % number of features to use, ranked by dimreduction technique
dim_reduce='bi';
red_cutoff=.5;
outlier_detect=1;
kfolds=2;
cv_sim=50;
merge=.4;

if mod(nparams,2)>0
	error('ephysPipeline:argChk','Parameters must be specified as parameter/value pairs!');
end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'wavelet_method'
			wavelet_method=varargin{i+1};
		case 'wavelet_mpca'
			wavelet_mpca=varargin{i+1};
		case 'fs'
			fs=varargin{i+1};
		case 'use_spiketime'
			use_spiketime=varaargin{i+1};
		case 'maxcoeffs'
			maxcoeffs=varargin{i+1};
		case 'clust_choice'
			clust_choice=varargin{i+1};
		case 'nfeatures'
			nfeatures=varargin{i+1};
		case 'dim_reduce'
			dim_reduce=varargin{i+1};
		case 'red_cutoff'
			red_cutoff=varargin{i+1};
		case 'outlier_detect'
			outlier_detect=varargin{i+1};
		case 'outlier_cutoff'
			outlier_cutoff=varargin{i+1};
		case 'kfolds'
			kfolds=varargin{i+1};
		case 'merge'
			merge=varargin{i+1};
	end
end

% TODO: chi2 test for dealing with outliers (simply take residuals and use chi2 CDF)
% need to deal with cell input (multiple trials), convert input to big matrix
% and spit out trial number

disp(['Wavelet dimensionality reduction method:  ' wavelet_method]);
disp(['N wavelet coefficients:  ' num2str(maxcoeffs)]);
disp(['MPCA ' num2str(wavelet_mpca)]);
disp(['GMM mode selection method:  ' clust_choice]);
disp(['All dimensionality reduction method:  ' dim_reduce]);
disp(['Redundancy cutoff:  ' num2str(red_cutoff)]);
disp(['Outlier detection:  ' num2str(outlier_detect)]);
disp(['Outlier cutoff:  ' num2str(outlier_cutoff)]);
disp(['Merge cutoff (Bhattacharyya distance):  ' num2str(merge)]);

spikewindows=[];
spiketimes=[];
spikeifr=[];
trialnum=[];
spike_data=[];

cv_flag=(lower(clust_choice(1))=='c');

% use ifr as a clustering feature

if iscell(SPIKEWINDOWS)
	for i=1:length(SPIKEWINDOWS)

		[samples,trials]=size(SPIKEWINDOWS{i});

		spikewindows=[spikewindows SPIKEWINDOWS{i}];
		spiketimes=[spiketimes SPIKETIMES{i}];

		ifr_tmp=[];
		padded_spikes=[-inf SPIKETIMES{i} inf];
		for j=2:length(padded_spikes)-1
			
			curr_spike=padded_spikes(j);
			next_spike=min(padded_spikes(padded_spikes>padded_spikes(j)))-curr_spike;
			prev_spike=curr_spike-max(padded_spikes(padded_spikes<padded_spikes(j)));
			ifr_tmp(j-1)=1/(min(next_spike,prev_spike)/fs); % isi is the time from the closest spike, before or after
		
		end

		spikeifr=[spikeifr;ifr_tmp(:)];
		[samples trials]=size(SPIKEWINDOWS{i});
		trialnum=[trialnum;repmat(i,trials,1)];

	end
else

	spikewindows=SPIKEWINDOWS;
	spiketimes=SPIKETIMES;

	ifr_tmp=[];
	padded_spikes=[-inf spiketimes inf];
	
	for j=2:length(padded_spikes)-1

		curr_spike=padded_spikes(j);
		next_spike=min(padded_spikes(padded_spikes>padded_spikes(j)))-curr_spike;
		prev_spike=curr_spike-max(padded_spikes(padded_spikes<padded_spikes(j)));
		ifr_tmp(j-1)=(1/min(next_spike,prev_spike)/fs); % isi is the min spike distance
	end

	spikeifr=ifr_tmp(:);

	[samples,trials]=size(SPIKEWINDOWS);
	trialnum=repmat(1,trials,1);
end


[samples,trials]=size(spikewindows);

% take the time vector and expand to interpolated fs

timepoints=[1:samples]';

% need to check for wavelet toolbox, otherwise we would need to use PCA or standard features
% include an option for mpca
% calculate other features and re-evaluate the correlation-coefficients and multimodality

if license('test','Wavelet_Toolbox')

	disp('Computing wavelet coefficients...');
	[coeffs]=get_wavelet_coefficients(spikewindows,maxcoeffs,...
		'method',wavelet_method,'mpca',wavelet_mpca); % pick multi-modal method and weighted PCA
	spike_data=[spike_data coeffs];
end

disp('Computing PCA...');
[coef score]=princomp(spikewindows');

% take first ten

spike_data=[spike_data score(:,1:10)];

% compute geometric measures (not working so well, leaving out for now) [9/20/12 JM]

%disp('Computing geometric features...');
%[coef]=get_geometric_coefficients(spikewindows);
%
%% remove features with NaN
%
%todel=[];
%for i=1:size(coef,2)
%	if sum(isnan(coef(:,i)))>0
%		todel(end+1)=i;
%	end
%end

%coef(:,todel)=[];
%spike_data=[spike_data coef];

nbins=sqrt(trials);

for i=1:size(spike_data,2)

	% center the data points

	testpoints=spike_data(:,i);
	samplemedian=median(testpoints);
	samplevar=median(abs(testpoints-samplemedian))./.6745;

	testpoints=testpoints-samplemedian;
	
	spike_data(:,i)=testpoints;

	samplemedian=median(testpoints);
	samplevar=median(abs(testpoints-samplemedian))./.6745;

	% don't need to remove outliers, simply control for cases with large skew

	%testpoints(testpoints>samplemedian+4*samplevar)=[];
	
	% get the empirical pdf, kernel density estimate

	x=linspace(min(testpoints),max(testpoints),nbins);

	count=ksdensity(testpoints,x);

	fx_pdf=count./sum(count);
	fx_cdf=cumsum(fx_pdf);

	% evaluate normpdf

	gx_pdf=normpdf(x,samplemedian,samplevar);
	gx_pdf=gx_pdf./sum(gx_pdf);
	gx_cdf=cumsum(gx_pdf);

	normentropy=-sum(gx_pdf.*log2(gx_pdf+eps));
	obsentropy=-sum(fx_pdf.*log2(fx_pdf+eps));

	negentropy(i)=[normentropy-obsentropy];
	coeffbimodal(i)=(1+skewness(testpoints)^2)/(kurtosis(testpoints)+((3*(trials-1)^2)/((trials-2)*(trials-3))));
	ksstat(i)=max(abs(fx_pdf-gx_pdf));

	% if the skewness is too high, reject
	% if the skewness exceeds roughly .4 this estimate becomes highly unstable

	if abs(skewness(testpoints)>.4)
		coeffbimodal(i)=0;
	end

end

% take the top six dimensions
% use multimodal test to choose best features

switch lower(dim_reduce(1))
	case 'n'
		[val loc]=sort(negentropy,'descend');
	case 'b'
		[val loc]=sort(coeffbimodal,'descend');
	case 'k'
		[val loc]=sort(ksstat,'descend');
end

% take 3 extra features

if nfeatures<length(loc)
	loc=loc(1:nfeatures);
end

spike_data=spike_data(:,loc);

% need better criteria for optimizing redundancy without losing informative dimensions 
% maybe MIFS or variant using coefficient of bimodality 

features=length(loc);

disp(['Using features ' num2str(loc) ' for clustering']);

% remove all correlated features

for i=2:features

	similarity=[];

	for j=1:i-1

		% get covariance matrix

		covmat=cov(spike_data(:,i),spike_data(:,j));

		% get corrcoef
		% take abs value, don't wait perfectly anti-correlated features...

		similarity(j)=abs(covmat(1,2)/(sqrt(covmat(1,1)*covmat(2,2))));
	end

	[redundancy(i-1) newloc(i-1)]=max(similarity);
end

redundant_dims=find(redundancy>red_cutoff);

disp(['Will remove the following redundant dimensions: ' num2str(redundant_dims+1)]);

if length(redundant_dims)==size(spike_data,2)
	redundant_dims=redundant_dims(1:size(spike_data,2)-1);
end

spike_data(:,redundant_dims+1)=[];

% TODO:  include projection pursuit as an option after it has been thoroughly vetted...

if use_spiketime
	spike_data=[spike_data spiketimes];
end

% now compute the combination of 6 features that is the least redundant
% clustering, fit a Gaussian mixture model, first find the appropriate number of modes using
% Akaike Information Criterion (AIC)

options=statset('Display','off');
obj_fcn=[];

clustnum=1:7;
[datapoints,features]=size(spike_data);
if datapoints<=features
	warning('ephysPipeline:spikesort:toofewspikes','Too few spikes to fit');
	return;
end

% gaussian mixture seems to work better than fcm
% may also want to check for stats toolbox, fall back on kmeans perhaps

kmeans_labels=zeros(datapoints,length(clustnum));

if matlabpool('size')>0
	pctRunOnAll warning('off','stats:gmdistribution:FailedToConverge');
	pctRunOnAll warning('off','stats:kmeans:FailedToConverge');
else
	warning('off','stats:gmdistribution:FailedToConverge');
	warning('off','stats:kmeans:FailedToConverge');
end

if cv_flag

	% randomize or not?

	%trialpool=[1:datapoints];

	foldpoints=floor(datapoints/kfolds);

	for i=1:cv_sim

		% divide set into k partitions, test on rest

		trialpool=randperm(datapoints);
		testset{i}=trialpool(1:foldpoints);
		trainset{i}=setdiff(trialpool,testset{i});

	end
end

CVNLOGL=zeros(length(clustnum),1);

if ~cv_flag
	parfor i=1:length(clustnum)

		% first use kmeans to get a good starting point


		kmeans_labels(:,i)=kmeans(spike_data,clustnum(i),'replicates',5);
		testobj=gmdistribution.fit(spike_data,clustnum(i),'Regularize',1,...
			'Options',options,'replicates',1,'start',kmeans_labels(:,i));	
		
		AIC(i)=testobj.AIC;
		logl(i)=testobj.NlogL;
		BIC(i)=testobj.BIC;

		disp([ num2str(clustnum(i)) ' clusters']);

		disp([ 'AIC ' num2str(testobj.AIC)]) % Akaike information criterion
		disp([ 'BIC ' num2str(testobj.BIC)]) % Bayes information criterion

		% grab the variance for each dimension in K components

		K=testobj.NComponents;
		M=testobj.NDimensions;

		variance=[];

		for j=1:K
			variance(j)=sqrt(det(testobj.Sigma(:,:,j)));
		end

		% fuzzy hypervolume

		FHV(i)=sum(variance);

		% evidence density

		ED(i)=-logl(i)/FHV(i);

		% minimum description length

		L=K*(1+M+(((M+1)*M)/2))-1; % nparameters
		MDL(i)=logl(i)+((L)/2)*log(trials*M);

		disp([ 'FHV ' num2str(FHV(i))]);	
		disp([ 'ED ' num2str(ED(i))]);
		disp([ 'MDL ' num2str(MDL(i))]);

	end
else

	% cross-validated likelihood, slower than BIC and generally worse...

	for i=1:length(clustnum)

		testnlogl=[];
		newtestnlogl=[];

		parfor j=1:cv_sim

			% k-fold validation
			
			[kmeanslabels]=kmeans(spike_data(trainset{j},:),clustnum(i),'replicates',2);
			testobj=gmdistribution.fit(spike_data(trainset{j},:),clustnum(i),...
				'replicates',1,'regularize',1,'Options',options,'start',kmeanslabels);

			% cluster the test set and get the test likelihood (need to verify this number)

			[tmpidx,testnlogl(j),p,logpdf]=cluster(testobj,spike_data(testset{j},:));

		end

		CVNLOGL(i)=mean(testnlogl);
		disp([ 'CVNLOGL ' num2str(CVNLOGL(i))]);

	end
end

if matlabpool('size')>0
	pctRunOnAll warning('on','stats:gmdistribution:FailedToConverge');
	pctRunOnAll warning('on','stats:kmeans:FailedToConverge');
else
 	warning('on','stats:gmdistribution:FailedToConverge');
	warning('on','stats:kmeans:FailedToConverge');
end


% AIC and BIC have worked miserably here, simply using the elbow of the log-likelihood
% decided to use fuzzy hypervolume or scaled likelihood, pretty stable (8/7/2012)

switch lower(clust_choice(1))

	case 'b'

		% BIC

		[val loc]=min(BIC);
		nclust=clustnum(loc);
	case 'a'

		% AIC

		[val loc]=min(AIC);
		nclust=clustnum(loc);
	case 'k'

		% logl knee

		secondderiv=diff(diff(logl));
		x=clustnum(3:end);
		[val,loc]=max(secondderiv); % maximum derivative in log-likelihood over k
		nclust=x(loc(1));

	case 'e'

		% ED Measure, log likelihood scaled by hypervolume

		[val loc]=min(ED);
		nclust=clustnum(loc);

	case 'f'

		% hypervolume alone

		[val loc]=min(FHV);
		nclust=clustnum(loc);

	case 'm'

		% minimum description length

		[val loc]=min(MDL);
		nclust=clustnum(loc);

	case 'c'

		[val loc]=min(CVNLOGL);
		nclust=clustnum(loc);

	otherwise
		error('ephysPipeline:autoclust:badnclustchoice','Did not understand nclust choice method!');
end

disp(['Will use ' num2str(nclust) ' clusters']);

warning('off','stats:gmdistribution:FailedToConverge');
warning('off','stats:kmeans:FailedToConverge');
[kmeanslabels]=kmeans(spike_data,nclust,'replicates',10);
testobj=gmdistribution.fit(spike_data,nclust,'Regularize',1,...
	'Options',options,'replicates',1,'start',kmeanslabels);
warning('on','stats:gmdistribution:FailedToConverge');
warning('on','stats:kmeans:FailedToConverge');

[idx,nlogl,P]=cluster(testobj,spike_data);

% assign output variables, garbage collection

%clear testobj;

% instead, get mean and covariances, project onto FLD, assess cluster quality and then 
% use outlier cutoff per Hill et al. (2011)

% for outliers take membership assignment and compute the integral
% of the residuals, this should follow a chi2 distribution N degrees of freedom

clusters=unique(idx);

% number of spikes per cluster is simply the number of labels

for i=1:length(clusters)
	nspikes(i)=sum(idx==clusters(i));
end

[val loc]=sort(nspikes,'descend');

% make the number contiguous and sort by number of spikes, descending

LABELS=zeros(size(idx));

for i=1:length(clusters)
	LABELS(idx==clusters(loc(i)))=i;	
end

if merge>0
	
	% use exp(-d), where d is Battacharyya distance to find clusters for potential merging

	count=1;

	while count>0 && length(clusters)>1

		clustpairs=nchoosek(1:length(clusters),2);
		clustdist=zeros(size(clustpairs,1),1);
		count=0;

		for i=1:size(clustpairs,1)

			c1=clustpairs(i,1);
			c2=clustpairs(i,2);

			m1=mean(spike_data(LABELS==c1,:))';
			m2=mean(spike_data(LABELS==c2,:))';
			cov1=cov(spike_data(LABELS==c1,:));
			cov2=cov(spike_data(LABELS==c2,:));

			mcov=(cov1+cov2)/2;

			dist=(1/8)*(m1-m2)'*inv(mcov)*(m1-m2)+...
				(1/2)*log((det(mcov))/(sqrt(det(cov1)*det(cov2))));

			clustdist(i)=exp(-dist);

			fprintf('Exp(-d) between %g and %g:\t%.3f\n',c1,c2,clustdist(i));

		end

		count=sum(clustdist>merge);

		if count>0
	
			[val loc]=max(clustdist);

			merge1=clustpairs(loc(1),1);
			merge2=clustpairs(loc(1),2);

			fprintf('Merging clusters %g and %g\n',merge1,merge2);
			
			if merge1<merge2
				LABELS(LABELS==merge2)=merge1;
			else
				LABELS(LABELS==merge1)=merge2;
			end

			% ensure the labeling is contiguous again

			idx=LABELS;
			clusters=unique(idx);
			LABELS=zeros(size(LABELS));

			for i=1:length(clusters)
				LABELS(idx==clusters(i))=i;
			end

			clusters=unique(LABELS);

		end

	end

end

% resort by number of spikes

idx=LABELS;
clusters=unique(idx);

% number of spikes per cluster is simply the number of labels

nspikes=[];
for i=1:length(clusters)
	nspikes(i)=sum(idx==clusters(i));
end

[val loc]=sort(nspikes,'descend');

% make the number contiguous and sort by number of spikes, descending

LABELS=zeros(size(idx));

for i=1:length(clusters)
	LABELS(idx==clusters(loc(i)))=i;	
end

% recluster with mean templates, use a residual test and dump outliers into
% a "junk" cluster

if outlier_detect

	% construct a template

	templates=zeros(samples,length(clusters));
	newlabels=zeros(size(LABELS));

	for i=1:length(clusters)
		templates(:,i)=mean(spikewindows(:,LABELS==i),2);
	end

	counter=zeros(1,length(clusters));

	for i=1:trials

		score=zeros(1,length(clusters));

		% l-infty norm

		for j=1:length(clusters)
			score(j)=max(abs(spikewindows(:,i)-templates(:,j)));
		end

		[val loc]=min(score);

		if val>outlier_cutoff
			newlabels(i)=NaN;
			counter(LABELS(i))=counter(LABELS(i))+1;
		else
			newlabels(i)=loc;
		end

	end

	for i=1:length(clusters)
		prcoutliers=100*(counter(i)./sum(LABELS==i));
		disp(['Cluster ' num2str(i) ' : ' num2str(prcoutliers) ' % outliers']);
	end

	LABELS=newlabels;
	clusters=unique(LABELS(~isnan(LABELS)));

	for i=1:length(clusters)
		clusterlocs=find(LABELS==i);
		CLUSTERPOINTS{i}=spike_data(clusterlocs,:);
	end

	clusterlocs=find(isnan(LABELS));
	OUTLIERS=spikewindows(:,clusterlocs);

	% any outliers were assigned the new cluster, don't want
	% to include in any further analysis (save for checking)

	idx=LABELS;
	clusters=unique(idx(~isnan(idx)));

	% in the odd case an entire cluster is wiped out... (should NOT happen)

	LABELS=zeros(size(idx));
	for i=1:length(clusters)
		LABELS(idx==clusters(i))=i;	
	end

	LABELS(isnan(idx))=NaN;

	%for i=1:length(clusters)

	%	clusterlocs=find(LABELS==i);

	%	% get spike coordinates for trialnumber and spikewindow
	%	% matrices

	%	CLUSTERPOINTS{i}=spike_data(clusterlocs,:);

	%	% find outliers

	%	% get the residuals

	%	[nclusttrials,dims]=size(CLUSTERPOINTS{i});
	%	residuals=mahal(CLUSTERPOINTS{i},CLUSTERPOINTS{i});

	%	% cutoff point

	%	cutoff=chi2inv(1-1/nclusttrials,dims);
	%	outliers=find(residuals>cutoff);

	%	% display percentage outliers per cluster

	%	disp(['Cluster ' num2str(i) ' % outliers: '...
	%		num2str(length(outliers)*100/nclusttrials)]);

	%	% store outlier windows in a separate cell array, or just 

	%	trialnum(clusterlocs(outliers))=[];

	%	CLUSTERPOINTS{i}(outliers,:)=[];
	%	OUTLIERS{i}=spikewindows(:,clusterlocs(outliers));
	%	spikewindows(:,clusterlocs(outliers))=[];
	%	spiketimes(clusterlocs(outliers))=[];

	%	% give outliers a new cluster identity, or delete

	%	LABELS(clusterlocs(outliers))=[];
	

	%end
end

TRIALS=trialnum;

% now assess the cluster quality ,
% take each cluster and check the FP and FN rate

% use l ratio of .05

for i=1:length(clusters)
	
	clusterlocs=find(LABELS==i);
	otherlocs=find(LABELS~=i);
	% get the feature data

	clusterpoints=spike_data(clusterlocs,:);
	otherpoints=spike_data(otherlocs,:);

	% l ratio is the sum inv chi2cdf of mahal distance of all other points over
	% n spikes

	nclustpoints=size(clusterpoints,1);

	mahaldist=mahal(otherpoints,clusterpoints);
	
	if length(otherpoints>0)
		STATS.lratio(i)=sum(1-chi2cdf(mahaldist.^2,features))/nclustpoints;
	else
		STATS.lratio(i)=NaN;
	end

	if length(mahaldist)>=nclustpoints
		sortmahal=sort(mahaldist.^2,'ascend');
		STATS.isod(i)=sortmahal(nclustpoints);
	else
		STATS.isod(i)=NaN;
	end

end

clear spike_data;

% resort labels by number of spikes, first vector is cluster id, and the second the 
% point where those labels end (e.g. 1 1 1 2 2 2 2 would return 1 2 and 3 7)

% return labels, and windows and ISI sorted by cluster IDX

[uniq_trial trial_boundary trial_group]=unique(trialnum);
trial_boundary=[0;trial_boundary];

for i=1:length(clusters)

	WINDOWS{i}=spikewindows(:,LABELS==i);
	spikeisitmp=[];

	for j=1:length(uniq_trial)
		
		% all spike times in this trial
		
		currtrial=spiketimes(trial_boundary(j)+1:trial_boundary(j+1));

		% now all spike ids from this trial

		currlabels=LABELS(trialnum==uniq_trial(j));

		% spike times for this cluster

		currtrial=currtrial(currlabels==clusters(i));
		currisi=(diff(currtrial));
		spikeisitmp=[spikeisitmp;currisi(:)];
	end

	ISI{i}=spikeisitmp; 
end
