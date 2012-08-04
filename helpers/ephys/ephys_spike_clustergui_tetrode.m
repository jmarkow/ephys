function [LABELS TRIALS ISI WINDOWS]=ephys_spike_clustergui_tetrode(SPIKEWINDOWS,SPIKETIMES,varargin)
%GUI for spike cluster cutting
%
%

% spikewindows', rows x samples, each row is a windowed spike waveform

TRIALS=[]; % by default we don't need this, unless input is over multiple trials
fs=25e3;
features={'P2P','width','energy','min','ISI','Spiketimes','PCA','wavelets'}; % possible features include, min, max, PCA, width
				     % energy and wavelet coefficients
outlier_cutoff=.5;				     
channel_labels=[];
nparams=length(varargin);
colors={'b','r','g','c','m','y','k','r','g','b'};
NDIMS=[];
LEGEND_LABELS={};
LABELS=[];
ISI={};
WINDOWS={};
CLUSTERS=[];
wavelets=10; % chooses top N non-normal wavelets according to either negentropy or KS test
%

wavelet_method='ks';
wavelet_mpca=1;


if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs!');
end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'fs'
			fs=varargin{i+1};
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

	end
end

% need to deal with cell input (multiple trials), convert input to big matrix
% and spit out trial number

% get the number of channels

spikewindows=[];
spiketimes=[];
spikeifr=[];
trialnum=[];
spike_data=[];

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

[samples,trials,channels]=size(spikewindows);


if isempty(channel_labels)
	channel_labels=[1:channels];
end

clear SPIKEWINDOWS;
TRIALS=trialnum;



% take the time vector and expand to 2*fs (from 25k to 50k)

timepoints=[1:samples]';

% expand to include features from all channels processed
% by convention let's keep each channel a separate column

features_all={'max','min','p2p','width','energy','PC1','PC2'}; % all features excluding IFR and spiketimes

for i=1:wavelets
	features_all{7+i}=['WC' num2str(i)];
end

features_status=zeros(7+wavelets,1,'int8');

spike_data=[];
for i=1:channels

	% cheap to compute standard features


	[max_value(:,i),max_sample]=max(spikewindows(:,:,i)',[],2); % will get first max if there are multiple peaks
	[min_value(:,i),min_sample]=min(spikewindows(:,:,i)',[],2); % will get first min if there are multiple peaks (maybe interpolate?)
	peak_to_peak(:,i)=max_value(:,i)-min_value(:,i);
	width(:,i)=abs(min_sample-max_sample);
	energy(:,i)=sum(spikewindows(:,:,i)'.^2,2);


	if any(strcmp('max',lower(features)))
		spike_data=[spike_data max_value];
		features_status(1)=1;
	end

	if any(strcmp('min',lower(features)))
		spike_data=[spike_data min_value];
		features_status(2)=1;
	end

	if any(strcmp('p2p',lower(features)))
		spike_data=[spike_data peak_to_peak];
		features_status(3)=1;
	end

	if any(strcmp('width',lower(features)))
		spike_data=[spike_data width];
		features_status(4)=1;
	end

	if any(strcmp('energy',lower(features)))
		spike_data=[spike_data energy];
		features_status(5)=1;
	end

	if any(strcmp('pca',lower(features)))
		[coef score variance t2]=princomp(spikewindows(:,:,i)');
		spike_data=[spike_data score(:,1:2)];
		features_status(6:7)=1;
	end

	if any(strcmp('wavelets',lower(features)))
		[coeffs]=get_wavelet_coefficients(spikewindows(:,:,i),wavelets,'method',wavelet_method,'mpca',wavelet_mpca);
		[m,ncoeffs]=size(coeffs);
		spike_data=[spike_data coeffs];
		features_status(8:end)=1;
	end
end

features_status=features_status(1:7+ncoeffs);
features_all=features_all(1:7+ncoeffs);

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
	spike_data=[spike_data spikeifr(:)];
	property_names{end+1}='Spike ISI';

end

main_window=figure('Visible','off','Position',[360,500,700,600],'Name','Data Plotter','NumberTitle','off');
plot_axis=axes('Units','pixels','Position',[50,50,400,400]);

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

rows=ceil(length(property_names)/5);

i=1;
while i<=length(property_names)
	row=ceil(i/5);
	column=mod(i,5);
	if column==0, column=5; end
	cluster_data_check{i}=uicontrol('Style','checkbox',...
		'String',property_names{i},...
		'Value',i==1,'Position',[50+column*60,600-row*35,70,25]);
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
	push_replot_save,push_draw_mode,push_recluster],'Units','Normalized');
movegui(main_window,'center')

set(main_window,'Visible','On');



uiwait(main_window);

%% Callbacks

% this callback changes the plot and returns the sum of the distances
% from the centroid for each point in a cluster

% change the plot if we change any of our dimensions, DO NOT RECLUSTER!

function change_plot(varargin)

% get the number of dimensions for the plot (number of principal components)

cla;

viewdim(1)=get(pop_up_x,'value');
viewdim(2)=get(pop_up_y,'value');
viewdim(3)=get(pop_up_z,'value');

view_data=spike_data(:,viewdim);

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
	options=statset('Display','off');

	clustnum=2:5;
	if datapoints<=features
		disp('Too few spikes to fit');
		return;
	end

	% gaussian mixture seems to work better than fcm

	if strcmp(lower(clusterchoice),'auto')

		parfor i=1:length(clustnum)

			testobj=gmdistribution.fit(cluster_data,clustnum(i),'Regularize',1,'Options',options);
			
			AIC(i)=testobj.AIC;
			logl(i)=testobj.NlogL;
			disp([ num2str(clustnum(i)) ' clusters']);

			%disp(['Partition coefficient ' num2str(partition_coef(i))]);
			disp([ 'AIC ' num2str(testobj.AIC)]) % Akaike information criterion
			%disp([ 'BIC ' num2str(testobj.BIC)]) % Bayes information criterion

		end

		[val,loc]=max(diff(diff(logl))); % maximum derivative in log-likelihood over k
		nclust=clustnum(loc);
		CLUSTERS=nclust;

	else
		CLUSTERS=str2num(clusterchoice);

	end

	disp(['Will use ' num2str(CLUSTERS) ' clusters']);

	testobj=gmdistribution.fit(cluster_data,CLUSTERS,'Regularize',1,'Options',options);
	[idx,nlogl,P]=cluster(testobj,cluster_data);
	%[center,u,obj_fcn]=fcm(cluster_data,CLUSTERS,[NaN NaN NaN 0]);

	counter=1;

	for i=1:datapoints

		%[membership(i),idx(i)]=max(u(:,i)); % take posterior probability of the chosen cluster
							% given observation i as the measure of "membership"
		membership(i)=P(i,idx(i));

		if membership(i)<outlier_cutoff
			idx(i)=CLUSTERS+1; % assign new "junk cluster"
			counter=counter+1;
		end
	end

	disp([ num2str(counter) ' outliers']);

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

	cla;
	LABELS=ones(datapoints,1);
	response=[];
	counter=2;

	plot(view_data(:,1),view_data(:,2),'o','markerfacecolor',colors{1});view(2);
	hold on
	disp('Select the corners of the enclosing polygon then press RETURN to continue...');
	hold off;

	while isempty(response)

		[xv,yv]=ginput;
		k=convhull(xv,yv);	
		hold on;
		plot(xv(k),yv(k),'b-','linewidth',1.25);
		choice=inpolygon(view_data(:,1),view_data(:,2),xv(k),yv(k));
		LABELS(choice==1)=counter;
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

[uniq_trial trial_boundary trial_group]=unique(trialnum);
trial_boundary=[1;trial_boundary];

for i=1:CLUSTERS

	spikewintmp=spikewindows(:,LABELS==i,1);

	% need to collect isi within trial, don't count the first spike!
	
	spikeifrtmp=[];

	for j=1:length(uniq_trial)
		
		% all spike times in this trial
		
		currtrial=spiketimes(trial_boundary(j):trial_boundary(j+1));

		% now all spike ids from this trial

		currlabels=LABELS(trialnum==uniq_trial(j));

		% spike times for this cluster

		currtrial=currtrial(currlabels==i);

		currisi=(diff(currtrial)); % isi in msec
		spikeifrtmp=[spikeifrtmp;currisi(:)];
	end

	ISI{i}=spikeifrtmp;
	WINDOWS{i}=spikewintmp;

end

change_plot();

end	

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

	plot(([1:samples]./fs)*1e3,WINDOWS{i});

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
