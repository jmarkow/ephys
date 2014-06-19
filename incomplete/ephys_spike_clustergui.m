function [WINDOWS TIMES TRIALS ISI STATS OUTLIERS SPIKEDATA MODEL]=ephys_spike_clustergui_tetrode(SPIKES,NOISEDATA,PARAMETERS,varargin)
%GUI for spike cluster cutting
%
%

% spikewindows', rows x samples, each row is a windowed spike waveform

% The functions are all nested inside the main function but not each other
% thus by default all variables declared in the main function ARE GLOBAL

nparams=length(varargin);

% all features excluding IFR and spiketimes

features_all={'max','min','ne','^2','neo','wid','pgrad','ngrad','PC1','PC2','PC3','PC4'}; 
features={'PCA','pose','nege','posgrad','neggrad','min','width','ISI'}; 

% possible features include, min, max, PCA, width, energy and wavelet coefficients
 
channel_labels=[];
colors={'b','r','g','c','m','y','r','g','b'};
outliercolor='k';

NDIMS=2;

LEGEND_LABELS={};
LABELS=[];

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
		case 'pcs'
			pcs=varargin{i+1};
		case 'garbage'
			garbage=varargin{i+1};
		case 'proc_fs'
			proc_fs=varargin{i+1};
		case 'noisewhiten'
			noisewhiten=varargin{i+1};
		case 'align_method'
			align_method=varargin{i+1};
		case 'spike_window'
			spike_window=varargin{i+1};
	end
end


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

% expand to include features from all channels processed
% by convention let's keep each channel a separate column

spike_data=[];
property_names={};

% downsample spikes back to original FS

downfact=interpolate_fs/proc_fs;

if mod(downfact,1)~=0
	error('ephyspipeline:templatesortexact:baddownfact',...
		'Need to downsample by an integer factor');
end

clusterspikewindows=downsample(clusterspikewindows,downfact);


% cheap to compute standard features


if nchannels==1

	geom_features=get_geometric_coefficients(clusterspikewindows);
	
	if any(strcmp('max',lower(features)))
		spike_data=[spike_data geom_features(:,1)];
		property_names{end+1}='max';
	end

	if any(strcmp('min',lower(features)))
		spike_data=[spike_data geom_features(:,2)];
		property_names{end+1}='min';
	end

	if any(strcmp('pose',lower(features)))
		spike_data=[spike_data geom_features(:,3)];
		property_names{end+1}='pose';
	end

	if any(strcmp('nege',lower(features)))
		spike_data=[spike_data geom_features(:,4)];
		property_names{end+1}='nege';
	end

	if any(strcmp('tote',lower(features)))
		spike_data=[spike_data geom_features(:,5)];
		property_names{end+1}='tote';
	end

	if any(strcmp('neo',lower(features)))
		spike_data=[spike_data geom_features(:,6)];
		property_names{end+1}='neo';
	end

	if any(strcmp('width',lower(features)))
		spike_data=[spike_data geom_features(:,7)];
		property_names{end+1}='width';
	end

	if any(strcmp('posgrad',lower(features)))
		spike_data=[spike_data geom_features(:,8)];
		property_names{end+1}='posgrad';
	end

	if any(strcmp('neggrad',lower(features)))
		spike_data=[spike_data geom_features(:,9)];
		property_names{end+1}='neggrad';
	end
end

outlierpoints=[];
if any(strcmp('pca',lower(features)))
	newmodel=gmem(clusterspikewindows',[],1,'garbage',1,'merge',0,'debug',0);
	[v,d]=eigs(newmodel.sigma(:,:,1));
	newscore=-clusterspikewindows'*v;
	spike_data=[spike_data newscore(:,1:pcs)];

	% these comprise the outliers before the projection... set to >1 to include all (default for now)

	outlierpoints=newmodel.R(:,2)>=2;

	for i=1:pcs
		property_names{end+1}=['PC ' num2str(i)];
	end
end
	
spike_data(isnan(spike_data))=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GUI setup

main_window=figure('Visible','on','Position',[360,500,700,600],'Name','Data Plotter','NumberTitle','off');
plot_axis=axes('Units','pixels','Position',[50,50,425,425]);

pop_up_x= uicontrol('Style','popupmenu',...
	'String',property_names,...
	'Position',[425,90,75,25],'call',@change_plot,...
	'Value',1);
pop_up_x_text= uicontrol('Style','text',...
	'String','X',...
	'Position',[430,130,50,45]);

pop_up_y= uicontrol('Style','popupmenu',...
	'String',property_names,...
	'Position',[520,90,75,25],'call',@change_plot,...
	'Value',2);
pop_up_y_text= uicontrol('Style','text',...
	'String','Y',...
	'Position',[525,130,50,45]);

pop_up_z= uicontrol('Style','popupmenu',...
	'String',property_names,...
	'Position',[620,90,75,25],'call',@change_plot,...
	'Value',3);
pop_up_z_text= uicontrol('Style','text',...
	'String','Z',...
	'Position',[625,130,50,45]);

pop_up_clusters= uicontrol('Style','popupmenu',...
	'String',{'1','2','3','4','5','6','7','8','9'},...
	'Position',[500,210,75,25]);
pop_up_clusters_text= uicontrol('Style','text',...
	'String','Number of Clusters',...
	'Position',[525,250,100,45]);

push_replot_save= uicontrol('Style','pushbutton',...
	'String','Show cluster stats',...
	'Position',[500,40,100,25],'call',@show_stats);

push_recluster= uicontrol('Style','pushbutton',...
	'String','Recluster',...
	'Position',[500,550,100,35],'value',0,...
	'Call',@change_cluster);

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

align([pop_up_clusters,pop_up_clusters_text,push_replot_save],'Center','None');
align([pop_up_x,pop_up_x_text],'Center','None');
align([pop_up_y,pop_up_y_text],'Center','None');
align([pop_up_z,pop_up_z_text],'Center','None');

change_cluster();
change_plot();

% run change_plot, which updates the plot according to the defaults

set([main_window,plot_axis,pop_up_x,pop_up_x_text,pop_up_y,pop_up_y_text,pop_up_z,...
	pop_up_z_text,pop_up_clusters,pop_up_clusters_text,...
	push_replot_save,push_recluster],'Units','Normalized');
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
clusters=unique(LABELS(LABELS>0));

if NDIMS==2
	for i=1:length(clusters)
		points=find(LABELS==clusters(i));
		h(:,i)=plot(view_data(points,1),view_data(points,2),...
			'o','markerfacecolor',colors{i},'markeredgecolor','none');hold on
	end

	points=find(LABELS==0);
	if ~isempty(points)
		h(:,length(clusters)+1)=plot(view_data(points,1),view_data(points,2),...
			'o','markerfacecolor',outliercolor,'markeredgecolor','none');hold on
	end
else
	for i=1:length(clusters)
		points=find(LABELS==clusters(i));
		h(:,i)=plot3(view_data(points,1),view_data(points,2),view_data(points,3),...
			'o','markerfacecolor',colors{i},'markeredgecolor','none');hold on

	end

	points=find(LABELS==0);
	if ~isempty(points)
		h(:,length(clusters)+1)=plot3(view_data(points,1),view_data(points,2),view_data(points,3),...
			'o','markerfacecolor',outliercolor,'markeredgecolor','none');hold on
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

clusterchoices=get(pop_up_clusters,'string');
clusterselection=get(pop_up_clusters,'value');
clusterchoice=clusterchoices{clusterselection};

% perform the kmeans analysis and return the labels, centroid coordinates,
% sum of all points in each cluster from their respective centroid and
% the distance of all points from all centroids

% start with one cluster, go up to 10 and check the within distance for all clusters

cluster_data=spike_data(:,dim);
[datapoints,features]=size(spike_data);

options=statset('Display','off');

clustnum=2:9;
if datapoints<=features
	disp('Too few spikes to fit');
	return;
end

% gaussian mixture seems to work better than fcm


startmu=[];
startcov=[];
mixing=[];

nclust=str2num(clusterchoice);

startobj=struct('mu',startmu,'sigma',startcov,'mixing',mixing);
idx=kmeans(cluster_data,nclust,'replicates',5);

%% set up initial model

mu=[];
for j=1:nclust
	startmu(j,:)=mean(cluster_data(idx==j,:))';
	startcov(:,:,j)=diag(var(cluster_data));
end

startobj.mu=startmu;
startobj.sigma=startcov;

for j=1:nclust
	startobj.mixing(j)=sum(idx==j)/length(idx);
end

clustermodel=gmem(cluster_data,startobj,nclust,...
		'garbage',garbage,'merge',smem,'debug',0);

MODEL=clustermodel;

idx=[];
for i=1:size(clustermodel.R,1)
	posteriors=clustermodel.R;
	[~,idx(i)]=max(posteriors(i,:));
end

if garbage
	garbageidx=find(clustermodel.garbage);
	idx(idx==garbageidx)=NaN;
end

LABELS=zeros(size(cluster_data,1),1);

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


OUTLIERS=storespikewindows(:,isnan(idx));

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

for i=1:nclust
	LABELS(idx==clusters(loc(i)))=i;	
end

% return labels, and windows and ISI sorted by cluster IDX

% clear the plot axis

% plot in either 2 or 3 dims
% turns out plot is MUCH faster than scatter, changed accordingly...

LEGEND_LABELS={};
for i=1:nclust
	LEGEND_LABELS{i}=['Cluster ' num2str(i)];
end

if garbage & any(isnan(idx))
	LEGEND_LABELS{end+1}='Outliers';
end

% compute any other stats we want, ISI, etc...

[WINDOWS TIMES TRIALS SPIKEDATA ISI STATS]=...
	check_clusterquality(storespikewindows,spiketimes,cluster_data,LABELS,trialnum,clustermodel);
change_plot();

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Show a window with spike stats

function show_stats(varargin)

% get the labels from the main_window

cluster=[];
cluster=struct('windows',{WINDOWS},'times',{TIMES},'trials',{TRIALS},'spikedata',{SPIKEDATA},'stats',{STATS},...
	'isi',{ISI},'model',MODEL);
cluster.parameters=PARAMETERS;
ephys_cluster_autostats(cluster,NOISEDATA,'fig_num',[],'clust_plot',0,'savemode',0)

end

end
