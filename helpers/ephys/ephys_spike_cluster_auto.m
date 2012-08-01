function [LABELS TRIALS ISI WINDOWS]=ephys_spike_cluster_auto(SPIKEWINDOWS,SPIKETIMES,varargin)
%automated spike clustering using a GMM with EM
%
%
%

% spikewindows', rows x samples, each row is a windowed spike waveform

if ~license('test','Statistics_Toolbox')
	error('Need statistics toolbox for clustering!');
end

LABELS=[];
TRIALS=[];
ISI=[];
WINDOWS=[];

sr=25e3;
interpolate=1;
interpolate_fs=50e3;
use_spiketime=0; % use spiketime as a clustering feature (usually helps if SNR is low)
nparams=length(varargin);
maxcoeffs=10; % number of wavelet coefficients to use (sorted by KS statistic)
outlier_cutoff=.5; % posterior probability cutoff for outliers (.6-.8 work well) [0-1, high=more aggresive]

if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs!');
end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'sr'
			sr=varargin{i+1};
		case 'interpolate'
			interpolate=varargin{i+1};
		case 'interpolate_fs'
			interpolate_fs=varargin{i+1};
		case 'use_spiketime'
			spiketime=varaargin{i+1};
		case 'maxcoeffs'
			maxcoeffs=varargin{i+1};
	end
end

% need to deal with cell input (multiple trials), convert input to big matrix
% and spit out trial number

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
			ifr_tmp(j-1)=1/(min(next_spike,prev_spike)/sr); % isi is the time from the closest spike, before or after
		
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
		ifr_tmp(j-1)=(1/min(next_spike,prev_spike)/sr); % isi is the min spike distance
	end

	spikeifr=ifr_tmp(:);

	[samples,trials]=size(SPIKEWINDOWS);
	trialnum=repmat(1,trials,1);
end

TRIALS=trialnum;

[samples,trials]=size(spikewindows);

expansion=interpolate_fs/sr;
interpspikes=zeros(samples*expansion,trials);

% take the time vector and expand to interpolated fs

timepoints=[1:samples]';
newtimepoints=linspace(1,samples,expansion*samples)';

% either sinc or spline interpolation

parfor i=1:trials
	%interpspikes(:,i)=sinc(newtimepoints(:,ones(size(timepoints)))-timepoints(:,ones(size(newtimepoints)))')*spikewindows(:,i);
	interpspikes(:,i)=spline(timepoints,spikewindows(:,i),newtimepoints);
end

% need to check for wavelet toolbox, otherwise we would need to use PCA or standard features

if license('test','Wavelet_Toolbox')
	[coeffs]=get_wavelet_coefficients(interpspikes,maxcoeffs,'method','bimodal'); % use bimodality coefficient to pick relevant dimensions
	spike_data=[spike_data coeffs];
else

	% fall back on PCA if necessary, just first two components seem to be most useful

	disp('Could not find the wavelet toolbox, falling back on PCA...');

	[coef score]=princomp(interpspikes');
	spike_data=[spike_data score(:,1:2)];

end


if use_spiketime
	spike_data=[spike_data spiketimes];
end

% clustering, fit a Gaussian mixture model, first find the appropriate number of modes using
% Akaike Information Criterion (AIC)

options=statset('Display','off');
obj_fcn=[];

clustnum=1:6;
[datapoints,features]=size(spike_data);
if datapoints<=features
	disp('Too few spikes to fit');
	return;
end

% gaussian mixture seems to work better than fcm
% may also want to check for stats toolbox, fall back on kmeans perhaps


parfor i=1:length(clustnum)

	testobj=gmdistribution.fit(spike_data,clustnum(i),'Regularize',1,'Options',options);
	
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

end

% what gives us the max AIC?

% max partition coefficient instead

%[val,loc]=max(AIC);
%[val,loc]=max(partition_coef);

% second derivative is defined over 2:length-1

for i=1:length(logl)-2
	secondderiv(i)=logl(i+2)+logl(i)-2*logl(i+1);
end

% AIC and BIC have worked miserably here, simply using the elbow of the log-likelihood

x=clustnum(2:end-1);
[val,loc]=max(secondderiv); % maximum derivative in log-likelihood over k
nclust=x(loc(1));

disp(['Will use ' num2str(nclust) ' clusters']);

testobj=gmdistribution.fit(spike_data,nclust,'Regularize',1,'Options',options);
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

for i=1:length(nclust)
	idx(idx==clusters(i))=i;	
end

% return labels, and windows and ISI sorted by cluster IDX

LABELS=idx;

[uniq_trial trial_boundary trial_group]=unique(trialnum);
trial_boundary=[1;trial_boundary];

for i=1:length(clusters)

	WINDOWS{i}=interpspikes(:,LABELS==i);

	spikeisitmp=[];

	for j=1:length(uniq_trial)
		
		% all spike times in this trial
		
		currtrial=spiketimes(trial_boundary(j):trial_boundary(j+1));

		% now all spike ids from this trial

		currlabels=LABELS(trialnum==uniq_trial(j));

		% spike times for this cluster

		currtrial=currtrial(currlabels==clusters(i));

		currisi=(diff(currtrial)); 
		spikeisitmp=[spikeisitmp;currisi(:)];
	end


	ISI{i}=spikeisitmp; 
end
