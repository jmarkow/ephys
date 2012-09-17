function [LABELS TRIALS ISI WINDOWS]=ephys_spike_cluster_auto(SPIKEWINDOWS,SPIKETIMES,varargin)
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

LABELS=[];
TRIALS=[];
ISI=[];
WINDOWS=[];

fs=25e3;
use_spiketime=0; % use spiketime as a clustering feature (usually helps if SNR is low)
nparams=length(varargin);
maxcoeffs=10; % number of wavelet coefficients to use (sorted by KS statistic)
wavelet_method='bi';
wavelet_mpca=0;
outlier_cutoff=.05; % posterior probability cutoff for outliers (.6-.8 work well) [0-1, high=more aggresive]
clust_choice='ed'; % knee is more tolerant, choice BIC (b) or AIC (a) for more sensitive clustering
nfeatures=10; % number of features to use, ranked by dimreduction technique
dim_reduce='bi';
red_cutoff=.1;

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
	end
end

% need to deal with cell input (multiple trials), convert input to big matrix
% and spit out trial number

disp(['Wavelet dimensionality reduction method:  ' wavelet_method]);
disp(['N wavelet coefficients:  ' num2str(maxcoeffs)]);
disp(['MPCA ' num2str(wavelet_mpca)]);
disp(['GMM mode selection method:  ' clust_choice]);

spikewindows=[];
spiketimes=[];
spikeifr=[];
trialnum=[];
spike_data=[];

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

TRIALS=trialnum;

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

% compute geometric measures

disp('Computing geometric features...');
[coef]=get_geometric_coefficients(spikewindows);

% remove features with NaN

todel=[];
for i=1:size(coef,2)
	if sum(isnan(coef(:,i)))>0
		todel(end+1)=i;
	end
end

coef(:,todel)=[];

spike_data=[spike_data coef];
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

for i=2:features

	jsd=[];

	for j=1:i-1
		jsd(j)=sqrt(kld(spike_data(:,i),spike_data(:,j),'jsd',1));
	end

	[redundancy(i-1) newloc(i-1)]=min(jsd);
end

redundant_dims=find(redundancy<red_cutoff);

disp(['Will remove the following redundant dimensions: ' num2str(redundant_dims+1)]);

if length(redundant_dims)<size(spike_data,2)
	spike_data(:,redundant_dims+1)=[];
end

% TODO:  include projection pursuit as an option after it has been thoroughly vetted...

if use_spiketime
	spike_data=[spike_data spiketimes];
end

% now compute the combination of 6 features that is the least redundant
% clustering, fit a Gaussian mixture model, first find the appropriate number of modes using
% Akaike Information Criterion (AIC)

options=statset('Display','off');
obj_fcn=[];

clustnum=2:7;
[datapoints,features]=size(spike_data);
if datapoints<=features
	warning('ephysPipeline:spikesort:toofewspikes','Too few spikes to fit');
	return;
end

% gaussian mixture seems to work better than fcm
% may also want to check for stats toolbox, fall back on kmeans perhaps

parfor i=1:length(clustnum)

	warning('off','stats:gmdistribution:FailedToConverge');
	testobj=gmdistribution.fit(spike_data,clustnum(i),'Regularize',1,'Options',options,'replicates',10);
	warning('on','stats:gmdistribution:FailedToConverge');

	%[center,u,obj_fcn]=fcm(spike_data,clustnum(i),[NaN NaN NaN 0]);
	
	% compute the partition coefficient, simply all membership indices squared and summed

	%partition_coef(i)=sum(sum(u.^2))/datapoints;

	AIC(i)=testobj.AIC;
	logl(i)=testobj.NlogL;
	BIC(i)=testobj.BIC;
	disp([ num2str(clustnum(i)) ' clusters']);

	%disp(['Partition coefficient ' num2str(partition_coef(i))]);
	disp([ 'AIC ' num2str(testobj.AIC)]) % Akaike information criterion
	disp([ 'BIC ' num2str(testobj.BIC)]) % Bayes information criterion

	[m,n,K]=size(testobj.Sigma);

	% grab the variance for each dimension in K components

	variance=[];

	for j=1:K
		variance(j)=sqrt(det(testobj.Sigma(:,:,j)));
	end

	% the hypervolume is defined by the product of the variances for a given components,
	% so then take the summed volume of all components as a measure of model fit (FHV)
	% another method is to take the hypervolume as a penalty term for the log-likelihood (ED)

	FHV(i)=sum(variance);
	ED(i)=-logl(i)/FHV(i);

	d=testobj.NDimensions;

	L=d+d*((d+1)/2);
	MDL(i)=logl(i)+((K*L)/2)*log(trials);
	
	disp([ 'FHV ' num2str(FHV(i))]);	
	disp([ 'ED ' num2str(ED(i))]);
	disp([ 'MDL ' num2str(MDL(i))]);

end

secondderiv=diff(diff(logl));

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

	otherwise
		error('ephysPipeline:autoclust:badnclustchoice','Did not understand nclust choice method!')
end

disp(['Will use ' num2str(nclust) ' clusters']);

warning('off','stats:gmdistribution:FailedToConverge');
testobj=gmdistribution.fit(spike_data,nclust,'Regularize',1,'Options',options,'replicates',5);
warning('on','stats:gmdistribution:FailedToConverge');

[idx,nlogl,P]=cluster(testobj,spike_data);

%[center,u,obj_fcn]=fcm(spike_data,nclust,[NaN NaN NaN 0]);
% place all points with posterior probability <cutoff in the same junk cluster

counter=1;
for i=1:datapoints

	%[membership(i),idx(i)]=max(u(:,i)); % take posterior probability of the chosen cluster
						% given observation i as the measure of "membership"
	membership(i)=P(i,idx(i));

	if membership(i)<outlier_cutoff
		idx(i)=nclust+1; % assign new "junk cluster"
		counter=counter+1;
	end
end

disp([ num2str(counter) ' outliers']);
clusters=unique(idx);

% number of spikes per cluster is simply the number of labels

for i=1:length(clusters)
	nspikes(i)=sum(idx==clusters(i));
end

[val loc]=sort(nspikes,'descend');

% make the number contiguous and sort by number of spikes, descending

LABLES=zeros(size(idx));

for i=1:length(clusters)
	LABELS(idx==clusters(loc(i)))=i;	
end

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
