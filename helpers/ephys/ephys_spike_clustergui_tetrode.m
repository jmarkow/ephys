function [LABELS TRIALS ISI WINDOWS STATS]=ephys_spike_clustergui_tetrode(SPIKEWINDOWS,SPIKETIMES,SPIKELESS,varargin)
%GUI for spike cluster cutting
%
%

% spikewindows', rows x samples, each row is a windowed spike waveform

% The functions are all nested inside the main function but not each other
% thus by default all variables declared in the main function ARE GLOBAL

nparams=length(varargin);

TRIALS=[]; % by default we don't need this, unless input is over multiple trials
fs=25e3;
interpolate_fs=200e3;

% all features excluding IFR and spiketimes

features_all={'max','min','ne','^2','neo','wid','pgrad','ngrad','PC1','PC2','PC3','PC4'}; 
features={'PCA','pose','nege','posgrad','neggrad','min','width','ISI','Spiketimes','wavelets'}; 

% possible features include, min, max, PCA, width, energy and wavelet coefficients
 
channel_labels=[];
colors={'b','r','g','c','m','y','k','r','g','b'};
NDIMS=[];
cluster_data=[];
LEGEND_LABELS={};
LABELS=[];
ISI={};
WINDOWS={};
CLUSTERS=[];
wavelets=10; % chooses top N non-normal wavelets according to either negentropy or KS test
wavelet_method='bi';
wavelet_mpca=1;
template_cutoff=[]; % l-infty norm cutoff for template clustering (max(abs) residual)

% the template cutoff could be defined by the 95th prctile of the abs(noise) magnitude

if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs!');
end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'fs'
			fs=varargin{i+1};
		case 'interpolate_fs'
			interpolate_fs=varargin{i+1};
		case 'features'
			features=varargin{i+1};
		case 'channel_labels'
			channel_labels=varargin{i+1};
		case 'wavelets'
			wavelets=varargin{i+1};
		case 'wavelet_method'
			wavelet_method=varargin{i+1};
		case 'wavelet_mpca'
			wavelet_mpca=varargin{i+1};
		case 'template_cutoff'
			template_cutoff=varargin{i+1};

	end
end


% get the noise estimate
% get the number of total noise samples, if it's over 10e3 downsample


noisedata=[];

for i=1:length(SPIKELESS)
	noisedata=[noisedata;downsample(SPIKELESS{i}(:),5)];
end

SIGMA_EST=var(noisedata); % conservative noise estimate, sup(variance) over trials

if isempty(template_cutoff)	
	template_cutoff=prctile(abs(noisedata),90);
end

disp(['Template outlier cutoff ' num2str(template_cutoff)]);

spikewindows=[];
spiketimes=[];
spikeifr=[];
trialnum=[];
spike_data=[];

% reformat windows and spiketimes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% IFR and spike times

if iscell(SPIKEWINDOWS)
	for i=1:length(SPIKEWINDOWS)

		[samples,trials,channels]=size(SPIKEWINDOWS{i});

		spikewindows=[spikewindows SPIKEWINDOWS{i}];
		spiketimes=[spiketimes SPIKETIMES{i}];

		ifr_tmp=[];
		padded_spikes=[-inf SPIKETIMES{i} inf];
		
		for j=2:length(padded_spikes)-1
			
			curr_spike=padded_spikes(j);
			next_spike=min(padded_spikes(padded_spikes>padded_spikes(j)))-curr_spike;
			prev_spike=curr_spike-max(padded_spikes(padded_spikes<padded_spikes(j)));
			ifr_tmp(j-1)=1/(min(next_spike,prev_spike)/fs); % ifr is the time from the closest spike, before or after
		
		end

		spikeifr=[spikeifr;ifr_tmp(:)];
		[samples,trials,channels]=size(SPIKEWINDOWS{i});
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

	[samples,trials,channels]=size(SPIKEWINDOWS);
	trialnum=repmat(1,trials,1);
end

[samples,trials,channels]=size(spikewindows);

if isempty(channel_labels)
	channel_labels=[1:channels];
end

clear SPIKEWINDOWS;
TRIALS=trialnum;

timepoints=[1:samples]';

% expand to include features from all channels processed
% by convention let's keep each channel a separate column

for i=1:wavelets
	features_all{end+1}=['WC' num2str(i)];
end

features_status=zeros(length(features_all),1,'int8');

spike_data=[];

for i=1:channels

	% cheap to compute standard features


	geom_features=get_geometric_coefficients(spikewindows);

	%[max_value,max_sample]=max(spikewindows(:,:,i)',[],2); % will get first max if there are multiple peaks
	%[min_value,min_sample]=min(spikewindows(:,:,i)',[],2); % will get first min if there are multiple peaks (maybe interpolate?)
	%peak_to_peak=max_value-min_value;
	%width=abs(min_sample-max_sample);
	%energy=sum(spikewindows(:,:,i)'.^2,2);


	if any(strcmp('max',lower(features)))
		spike_data=[spike_data geom_features(:,1)];
		features_status(1)=1;
	end

	if any(strcmp('min',lower(features)))
		spike_data=[spike_data geom_features(:,2)];
		features_status(2)=1;
	end

	if any(strcmp('pose',lower(features)))
		spike_data=[spike_data geom_features(:,3)];
		features_status(3)=1;
	end

	if any(strcmp('nege',lower(features)))
		spike_data=[spike_data geom_features(:,4)];
		features_status(4)=1;
	end

	if any(strcmp('tote',lower(features)))
		spike_data=[spike_data geom_features(:,5)];
		features_status(5)=1;
	end

	if any(strcmp('neo',lower(features)))
		spike_data=[spike_data geom_features(:,6)];
		features_status(6)=1;
	end

	if any(strcmp('width',lower(features)))
		spike_data=[spike_data geom_features(:,7)];
		features_status(7)=1;
	end

	if any(strcmp('posgrad',lower(features)))
		spike_data=[spike_data geom_features(:,8)];
		features_status(8)=1;
	end

	if any(strcmp('neggrad',lower(features)))
		spike_data=[spike_data geom_features(:,9)];
		features_status(9)=1;
	end

	if any(strcmp('pca',lower(features)))
		[coef score variance t2]=princomp(spikewindows(:,:,i)');
		spike_data=[spike_data score(:,1:4)];
		features_status(10:13)=1;
	end

	if any(strcmp('wavelets',lower(features)))
		[coeffs]=get_wavelet_coefficients(spikewindows(:,:,i),wavelets,'method',wavelet_method,'mpca',wavelet_mpca);
		[m,ncoeffs]=size(coeffs);
		spike_data=[spike_data coeffs];
		features_status(13:end)=1;
	end
	
end

features_status=features_status(1:13+ncoeffs);
features_all=features_all(1:13+ncoeffs);

active_features=find(features_status);

property_names={};
for i=1:channels
	for j=1:length(active_features)
		property_names{end+1}=[ features_all{active_features(j)} ' CH ' num2str(channel_labels(i)) ];
	end
end

if any(strcmp('spiketimes',lower(features)))
	spike_data=[spike_data spiketimes(:)];
	property_names{end+1}='Spike times';
end

if any(strcmp('isi',lower(features)))
	spike_data=[spike_data log(spikeifr(:))];
	property_names{end+1}='Spike ISI';

end

spike_data(isnan(spike_data))=0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GUI setup

main_window=figure('Visible','off','Position',[360,500,700,600],'Name','Data Plotter','NumberTitle','off');
plot_axis=axes('Units','pixels','Position',[50,50,425,425]);

pop_up_x= uicontrol('Style','popupmenu',...
	'String',property_names,...
	'Position',[400,90,75,25],'call',@change_plot,...
	'Value',1);
pop_up_x_text= uicontrol('Style','text',...
	'String','X',...
	'Position',[405,130,50,45]);

pop_up_y= uicontrol('Style','popupmenu',...
	'String',property_names,...
	'Position',[495,90,75,25],'call',@change_plot,...
	'Value',2);
pop_up_y_text= uicontrol('Style','text',...
	'String','Y',...
	'Position',[500,130,50,45]);

pop_up_z= uicontrol('Style','popupmenu',...
	'String',property_names,...
	'Position',[595,90,75,25],'call',@change_plot,...
	'Value',3);
pop_up_z_text= uicontrol('Style','text',...
	'String','Z',...
	'Position',[600,130,50,45]);

pop_up_clusters= uicontrol('Style','popupmenu',...
	'String',{'auto','1','2','3','4','5','6','7','8','9'},...
	'Position',[475,210,75,25]);
pop_up_clusters_text= uicontrol('Style','text',...
	'String','Number of Clusters',...
	'Position',[500,250,100,45]);

push_replot_save= uicontrol('Style','pushbutton',...
	'String','Show cluster stats',...
	'Position',[500,40,100,25],'call',@show_stats);

push_draw_mode= uicontrol('Style','pushbutton',...
	'String','Draw mode (x and y only)',...
	'Position',[500,450,100,35],'value',0,...
	'Call',@change_cluster);

push_recluster= uicontrol('Style','pushbutton',...
	'String','Recluster',...
	'Position',[500,550,100,35],'value',0,...
	'Call',@change_cluster);

push_templatecluster= uicontrol('Style','pushbutton',...
	'String','Template cluster',...
	'Position',[500,350,100,35],'value',0,...
	'Call',@template_cluster);

rows=ceil(length(property_names)/5);

i=1;
while i<=length(property_names)
	row=ceil(i/7);
	column=mod(i,7);
	if column==0, column=7; end
	cluster_data_check{i}=uicontrol('Style','checkbox',...
		'String',property_names{i},...
		'Value',i==1,'Position',[5+column*60,600-row*35,70,25]);
	set(cluster_data_check{i},'Units','Normalized')
	i=i+1;
end

% now align everything and send the main_window handle to the output
% so we can use the gui with uiwait (requires the handle as a return value)

align([pop_up_clusters,pop_up_clusters_text,push_replot_save,push_templatecluster],'Center','None');
align([pop_up_x,pop_up_x_text],'Center','None');
align([pop_up_y,pop_up_y_text],'Center','None');
align([pop_up_z,pop_up_z_text],'Center','None');

change_cluster();
change_plot();

% run change_plot, which updates the plot according to the defaults

set([main_window,plot_axis,pop_up_x,pop_up_x_text,pop_up_y,pop_up_y_text,pop_up_z,...
	pop_up_z_text,pop_up_clusters,pop_up_clusters_text,...
	push_replot_save,push_draw_mode,push_recluster,push_templatecluster],'Units','Normalized');
movegui(main_window,'center')

set(main_window,'Visible','On');
uiwait(main_window);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Callbacks

% this callback changes the plot and returns the sum of the distances
% from the centroid for each point in a cluster

% change the plot if we change any of our dimensions, DO NOT RECLUSTER!


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot change

function change_plot(varargin)

% get the number of dimensions for the plot

% clear the axes and draw the points

cla;

viewdim(1)=get(pop_up_x,'value');
viewdim(2)=get(pop_up_y,'value');
viewdim(3)=get(pop_up_z,'value');

view_data=spike_data(:,viewdim);
CLUSTERS=length(unique(LABELS));

if NDIMS==2
	for i=1:CLUSTERS
		points=find(LABELS==i);
		h(:,i)=plot(view_data(points,1),view_data(points,2),...
			'o','markerfacecolor',colors{i},'markeredgecolor','none');hold on
	end

else
	for i=1:CLUSTERS
		points=find(LABELS==i);
		h(:,i)=plot3(view_data(points,1),view_data(points,2),view_data(points,3),...
			'o','markerfacecolor',colors{i},'markeredgecolor','none');hold on

	end
end

grid on
view(NDIMS)

xlabel(property_names{viewdim(1)});
ylabel(property_names{viewdim(2)});
zlabel(property_names{viewdim(3)});

L=legend(h,LEGEND_LABELS,'Location','NorthEastOutside');legend boxoff
set(L,'FontSize',20,'FontName','Helvetica')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Recluster

function change_cluster(varargin)

% label everything
% use the dimensions ticked in the top box for clustering

dim=[];
for i=1:length(cluster_data_check)

	value=get(cluster_data_check{i},'Value');

	if value
		dim=[dim i];
	end

end

draw_mode=get(push_draw_mode,'value');

clusterchoices=get(pop_up_clusters,'string');
clusterselection=get(pop_up_clusters,'value');

clusterchoice=clusterchoices{clusterselection};

% perform the kmeans analysis and return the labels, centroid coordinates,
% sum of all points in each cluster from their respective centroid and
% the distance of all points from all centroids

% start with one cluster, go up to 10 and check the within distance for all clusters

cluster_data=spike_data(:,dim);
[datapoints,features]=size(spike_data);


if ~draw_mode

	if matlabpool('size')>0
		pctRunOnAll warning('off','stats:gmdistribution:FailedToConverge');
		pctRunOnAll warning('off','stats:kmeans:FailedToConverge');
	else
		warning('off','stats:gmdistribution:FailedToConverge');
		warning('off','stats:kmeans:FailedToConverge');
	end

	options=statset('Display','off');

	clustnum=2:9;
	if datapoints<=features
		disp('Too few spikes to fit');
		return;
	end

	% gaussian mixture seems to work better than fcm

	if strcmp(lower(clusterchoice),'auto')


		parfor i=1:length(clustnum)

			kmeans_labels(:,i)=kmeans(cluster_data,clustnum(i),'replicates',5);
			testobj=gmdistribution.fit(cluster_data,clustnum(i),'Regularize',1,...
				'Options',options,'replicates',1,'start',kmeans_labels(:,i));	

			AIC(i)=testobj.AIC;
			BIC(i)=testobj.BIC;
			logl(i)=testobj.NlogL;

			disp([ num2str(clustnum(i)) ' clusters']);
			disp([ 'AIC ' num2str(testobj.AIC)]) % Akaike information criterion
			disp([ 'BIC ' num2str(testobj.BIC)]) % Bayes information criterion

		end

		[val loc]=min(BIC);
		nclust=clustnum(loc);
		CLUSTERS=nclust;

	else
		CLUSTERS=str2num(clusterchoice);

	end

	disp(['Will use ' num2str(CLUSTERS) ' clusters']);

	[kmeanslabels]=kmeans(cluster_data,CLUSTERS,'replicates',10);
	testobj=gmdistribution.fit(cluster_data,CLUSTERS,'Regularize',1,...
		'Options',options,'replicates',1,'start',kmeanslabels);

	if matlabpool('size')>0
		pctRunOnAll warning('on','stats:gmdistribution:FailedToConverge');
		pctRunOnAll warning('on','stats:kmeans:FailedToConverge');
	else
		warning('on','stats:gmdistribution:FailedToConverge');
		warning('on','stats:kmeans:FailedToConverge');
	end

	[idx,nlogl,P]=cluster(testobj,cluster_data);
	counter=1;

	clusterlabels=unique(idx);
	CLUSTERS=length(clusterlabels);

	for i=1:length(clusterlabels)
		idx(idx==clusterlabels(i))=i;	
	end

	% return labels, and windows and ISI sorted by cluster IDX

	LABELS=idx;
	NDIMS=3;

	% clear the plot axis
else

	viewdim(1)=get(pop_up_x,'value');
	viewdim(2)=get(pop_up_y,'value');

	view_data=spike_data(:,viewdim);
	cluster_data=view_data;

	cla;
	LABELS=ones(datapoints,1);
	counter=2;

	plot(view_data(:,1),view_data(:,2),'o','markerfacecolor',colors{1});view(2);
	hold on
	disp('Select the corners of the enclosing polygon then press RETURN to continue...');
	hold off;
	response=[];

	while isempty(response)
	
		response2=[];

		while isempty(response2)

			[xv,yv]=ginput;
			k=convhull(xv,yv);	
			hold on;
			plot(xv(k),yv(k),'b-','linewidth',1.25);
			choice=inpolygon(view_data(:,1),view_data(:,2),xv(k),yv(k));
			LABELS(choice==1)=counter;
			response2=input('(D)one drawing or (c)ontinue?  ','s');

			switch lower(response2)
				case 'd'
					break;
				case 'c'
					response2=[];
				otherwise
					response2=[];
			end


		end

		response=input('(D)one clustering or (c)ontinue?  ','s');

		switch lower(response)
			case 'd'
				break;
			case 'c'
				response=[];
			otherwise
				response=[];
		end

		counter=counter+1;
	
	end	
	
	CLUSTERS=counter;
	NDIMS=3;

end

% plot in either 2 or 3 dims
% turns out plot is MUCH faster than scatter, changed accordingly...

LEGEND_LABELS={};
for i=1:CLUSTERS
	LEGEND_LABELS{i}=num2str(i);
end

% compute any other stats we want, ISI, etc...

collect_stats();
change_plot();

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Recluster

function template_cluster(varargin)

% get the templates for each cluster by extracting mean
% consider removing outliers first

CLUSTERS=length(unique(LABELS));

for i=1:CLUSTERS
	template(:,i)=mean(WINDOWS{i},2);
	mvartemplate(:,i)=2.25.*iqr(WINDOWS{i},2);
	vartemplate(i)=var(template(:,i));
end

% could slide the window as well, probably not necessary with precise alignment



tmp=zeros(size(spikewindows,2),CLUSTERS);
%cutoff=chi2inv(.95,size(spikewindows,1)-1)
%cutoff=300;
samples=size(spikewindows,1);

for i=1:size(spikewindows,2)
	
	for j=1:CLUSTERS
		tmp(i,j)=max(abs(spikewindows(:,i,1)-template(:,j)));
	end

	[val loc]=min(tmp(i,:));

	% outlier check

	LABELS(i)=loc;

	% chi square test per Wang et al.

	%residual=var(spikewindows(:,i,1)-template(:,loc));
	%chisquare(i)=(residual*(samples-1))/SIGMA_EST;

	residual=tmp(i,loc);

	% simply use the L-infty norm here

	if residual>template_cutoff
		LABELS(i)=CLUSTERS+1;
	end

end

% new cluster assignment is based on template clusters

idx=LABELS;
clusterlabels=unique(idx);
CLUSTERS=length(clusterlabels);

for i=1:length(clusterlabels)
	idx(idx==clusterlabels(i))=i;	
end

% return labels, and windows and ISI sorted by cluster IDX

LABELS=idx;
clusterlabels=unique(LABELS);
outliers=[];

LABELS(outliers)=length(clusterlabels)+1;

idx=LABELS;
clusterlabels=unique(idx);
CLUSTERS=length(clusterlabels);

for i=1:length(clusterlabels)
	idx(idx==clusterlabels(i))=i;	
end

LABELS=idx;
CLUSTERS=length(unique(LABELS));

collect_stats();
change_plot();

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Stats collection

function collect_stats(varargin)

[uniq_trial trial_boundary trial_group]=unique(trialnum);
trial_boundary=[0;trial_boundary];

idx=LABELS;
clusterlabels=unique(idx);
CLUSTERS=length(clusterlabels);

for i=1:length(clusterlabels)
	idx(idx==clusterlabels(i))=i;	
end

LABELS=idx;
CLUSTERS=length(unique(LABELS));

WINDOWS={};
ISI={};

for i=1:CLUSTERS

	spikewintmp=spikewindows(:,LABELS==i,1);

	% need to collect isi within trial, don't count the first spike!
	
	spikeifrtmp=[];

	for j=1:length(uniq_trial)
		
		% all spike times in this trial
		
		currtrial=spiketimes(trial_boundary(j)+1:trial_boundary(j+1));

		% now all spike ids from this trial

		currlabels=LABELS(trialnum==uniq_trial(j));

		% spike times for this cluster

		currtrial=currtrial(currlabels==i);

		currisi=(diff(currtrial));

		spikeifrtmp=[spikeifrtmp;currisi(:)];
	end

	ISI{i}=spikeifrtmp;
	WINDOWS{i}=spikewintmp;

end

% compute remaining stats

STATS.lratio=[];
STATS.isod=[];

for i=1:CLUSTERS
	
	clusterlocs=find(LABELS==i);
	otherlocs=find(LABELS~=i);
	
	% get the feature data

	clusterpoints=cluster_data(clusterlocs,:);
	otherpoints=cluster_data(otherlocs,:);

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

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Show a window with spike stats

function show_stats(varargin)

% get the labels from the main_window

statfig=figure('Visible','off');

counter=1;
binedges=[0:.1:10];

ymin=inf;
ymax=-inf;

for i=1:CLUSTERS

	ax(i)=subplot(CLUSTERS,2,counter);
	[samples,trials]=size(WINDOWS{i});

	if trials>100
		plotwindow=WINDOWS{i}(:,1:100);
	else
		plotwindow=WINDOWS{i}(:,1:trials);
	end	

	plot(([1:samples]./interpolate_fs)*1e3,plotwindow);

	ylabel('microvolts');
	xlabel('msec');
	axis tight
	ylimits=ylim();
	
	if ylimits(1)<ymin
		ymin=ylimits(1);
	end
	
	if ylimits(2)>ymax
		ymax=ylimits(2);
	end

	%[f,xi]=ksdensity(spikeifr);
	
	xlimits=xlim();
	xlim([0 xlimits(2)]);

	counter=counter+1;
	ax2=subplot(CLUSTERS,2,counter);
	density=histc((ISI{i}/fs)*1e3,binedges);
	bar(binedges,density,'histc');
	xlabel('ISI (msec)');
	ylabel('N');
	xlim([binedges(1) binedges(end)]);

	counter=counter+1;
end

linkaxes(ax,'xy');
set(ax(1),'ylim',[ymin ymax]);
linkaxes(ax2,'xy');

set(statfig,'visible','on');

end

end
