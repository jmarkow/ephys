function ephys_visual_sua(EPHYS_DATA,HISTOGRAM,CHANNELS,varargin)
%generates song-aligned single-unit rasters
%
%	[MUA TIME LABEL HISTOGRAM]=ephys_visual_sua(MIC_DATA,EPHYS_DATA,CHANNELS,varargin)
%
%	MIC_DATA
%	aligned microphrone traces from extracted_data.mat (should be the variable mic_data)
%
%	HISTOGRAM
%	contour histogram returned by ephys_visual_histogram.m (or loaded from histogram.mat)
%
%	CHANNELS
%	channel labels (i.e. the channel that corresponds to a given element in the cell array
%	ephys_data) from extracted_data.mat
%
%	the following may be specified as parameter/value pairs:
%
%		tetrode_channels
%		channels may act as a putative n-trode (can specify an arbitrary number of channels),
%		defaults to manual cluster cutting
%
%		car_exclude
%		carelectrodes to exclude from noise estimate
%
%		SR
%		data sampling rate (default: 25e3)
%		
%		noise
%		noise rejection method ('car' for common average 'nn' for nearest neighbor, or 'none',
%		default: 'none')
%
%		freq_range
%		vector with two elements to specify the frequency range (one element specifies low pass, default: [500 8e3])
%
%		savedir
%		directory to store results (default: pwd)
%
%		min_f
%		lowermost frequency to display for contour histogram (default: 1e3)
%
%		max_f
%		uppermost frequency to display for contour histogram (default: 10e3)
%
%		auto_clust
%		perform automated cluster cutting via fitting a GMM through EM (default: 0)
%
%		sigma_t
%		multiple of variance estimate for automatic threshold setting (uses the Quiroga formulation, default: 4)
%
%		singletrials
%		number of singletrials to plot
%
%		savedir
%		directory to store results (default: pwd)
%
%		subtrials
%		vector of trials to include in the analysis (default: all trials)
%
%
% see also ephys_visual_mua.m,ephys_visual_lfp_amp.m,ephys_visual_lfp_tf.m,ephys_spike_cluster_auto.m,ephys_spike_clustergui_tetrode.m,ephys_spike_detect.m


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETER COLLECTION %%%%%%%%%%%%%%%%%

if nargin<3
	error('ephysPipeline:suavis:notenoughparams','Need 3 arguments to continue, see documentation');
end

nparams=length(varargin);

if mod(nparams,2)>0
	error('ephysPipeline:argChk','Parameters must be specified as parameter/value pairs!');
end

SR=25e3;
noise='none'; % common-average reference for noise removal
car_exclude=[];
savedir=pwd;
min_f=1;
max_f=10e3;
hist_colors='jet';
figtitle='';
freq_range=[500 8000]; % frequency range for filtering
sort=1; % do we want to sort?
auto_clust=0;
tetrode_channels=[];
sigma_t=4;
jitter=4;
singletrials=10; % number of random single trials to plot per cluster
subtrials=[];
align='com'; % how to align spike waveforms can be min for minimum peak or COM
	     % COM seems a bit sloppy, may move to MIN
% parameter for smooth density estimate

channels=CHANNELS;
smooth_rate=1e3;
sigma=.0025;

colors={'b','r','g','c','m','y','k','r','g','b'};

for i=1:2:nparams
	switch lower(varargin{i})
		case 'sr'
			SR=varargin{i+1};
		case 'noise'
			noise=varargin{i+1};
		case 'sigma_t'
			sigma_t=varargin{i+1};
		case 'filtering'
			filtering=varargin{i+1};
		case 'savedir'
			savedir=varargin{i+1};
		case 'min_f'
			min_f=varargin{i+1};
		case 'max_f'
			max_f=varargin{i+1};
		case 'hist_colors'
			hist_colors=varargin{i+1};
		case 'car_exclude'
			car_exclude=varargin{i+1};
		case 'figtitle'
			figtitle=varargin{i+1};
		case 'freq_range'
			freq_range=varargin{i+1};
		case 'sort'
			sort=varargin{i+1};
		case 'channels'
			channels=varargin{i+1};
		case 'tetrode_channels'
			tetrode_channels=varargin{i+1};
		case 'auto_clust'
			auto_clust=varargin{i+1};
		case 'sigma_t'
			sigma_t=varargin{i+1};
		case 'subtrials'
			subtrials=varargin{i+1};
		case 'align'
			align=varargin{i+1};
		case 'jitter'
			jitter=varargin{i+1};
	end
end

startidx=max([find(HISTOGRAM.f<=min_f)]);

if isempty(startidx)
	startidx=1;
end

stopidx=min([find(HISTOGRAM.f>=max_f)]);

if isempty(stopidx)
	stopidx=length(HISTOGRAM.f);
end

[samples,ntrials,ncarelectrodes]=size(EPHYS_DATA);
proc_data=zeros(samples,ntrials,length(channels));

if isempty(subtrials)
	subtrials=1:ntrials;
end

TIME=[1:samples]./SR; % time vector for plotting

% kernel density estimate of smooth firing rate

kernedges=[-3*sigma:1/smooth_rate:3*sigma];
binedges=[0:(1/smooth_rate):samples/SR];

smooth_spikes.time=binedges;

if exist('normpdf')>0
	kernel=normpdf(kernedges,0,sigma);
else
	kernel=(1/(sigma*sqrt(2*pi)))*exp((-(kernedges-0).^2)./(2*sigma^2));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SIGNAL CONDITIONING %%%%%%%%%%%%%%%%

proc_data=ephys_denoise_signal(EPHYS_DATA,CHANNELS,channels,'method',noise,'car_exclude',car_exclude);
proc_data=ephys_condition_signal(proc_data,'s','freq_range',freq_range);

if ~isempty(tetrode_channels)
	tetrode_data=ephys_denoise_signal(EPHYS_DATA,CHANNELS,tetrode_channels,'method',noise,'car_exclude',car_exclude);
	tetrode_data=ephys_condition_signal(tetrode_data,'s','freq_range',freq_range);
else
	tetrode_data=[];
end

clear EPHYS_DATA;

proc_data=proc_data(:,subtrials,:);

if ~isempty(tetrode_channels)
	tetrode_data=tetrode_data(:,subtrials,:);
end
[samples,ntrials,newchannels]=size(proc_data);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SPIKE DETECTION %%%%%%%%%%%%%%%%%%%%

disp('Entering spike detection...');
disp(['Alignment method:  ' align]);

% need a noise cutoff...let's start with 3*std or Quiroga's measure

for i=1:length(channels)

	disp(['Electrode ' num2str(channels(i))])

	if ~isempty(tetrode_channels)
		disp(['Will use tetrode channels ' num2str(tetrode_channels) ' for sorting']);
	end

	% collect spikes

	tmp_spikewindows={};
	tmp_spiketimes={};

	sort_data=cat(3,proc_data(:,:,i),tetrode_data);

	parfor j=1:ntrials
		

		spike_pp=[];
		spikeless_data=[];

		threshold=sigma_t*median(abs(proc_data(:,j,i))/.6745)
		%disp([num2str(threshold)]);
		spike_pp=ephys_spike_detect(squeeze(sort_data(:,j,:)),threshold,'sr',SR,'visualize','n','align',align,'jitter',jitter);

		% after spike detect also collect trace without spikes

		[winsamples,nspikes,tetrodes]=size(spike_pp.abs.windows);

		if nspikes<2
			spikeless_data=NaN;
			continue;
		end

		% if tetrode then check each channel for significant spike

		window=floor((winsamples-1)/2); % how many samples to the left and right of the spike

		% should we delete?

		adjust=0;

		% compute the peak to peak of the mean waveform
	
		spikeless_data=proc_data(:,j,i);

		% eliminate each spike from the data

		for k=1:nspikes

			spike_point=spike_pp.abs.times(k)-adjust;

			if spike_point-window<1
				spikeless_data=NaN;
				break;
			end

			spikeless_data(spike_point-window:spike_point+window)=[];
			adjust=adjust+((window*2)+1);
		end

		tmp_spiketimes{j}=spike_pp.abs.times;
		tmp_spikewindows{j}=spike_pp.abs.windows;
		tmp_threshold(j)=threshold;
		tmp_spikeless{j}=spikeless_data;
	
	end

	clear sort_data;

	spikewindows{i}=tmp_spikewindows;
	spikeless{i}=tmp_spikeless;
	spiketimes{i}=tmp_spiketimes;
	threshold(i,:)=tmp_threshold;
	
	[winsamples,trials,tetrodes]=size(spikewindows{1}{1});

	if tetrodes>1
		tetrode_preview=figure('Name','Tetrode Preview','Visible','off');
		subplot(tetrodes,1,1);
		plot(spikewindows{1}{1}(:,:,1));
		title(['Channel ' num2str(CHANNELS(channels(i)))]);
		axis tight;
		for k=1:tetrodes-1
			subplot(tetrodes,1,k+1);
			plot(spikewindows{1}{1}(:,:,k+1));
			title(['Channel ' num2str(tetrode_channels(k))]);
			axis tight
		end
		set(tetrode_preview,'Visible','on');
		pause();
		close([tetrode_preview]);
	end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% simplest way to draw raster...plot([sample;sample]
% vector of times [120 130 200;120 130 200], then next vector specifies bottom and top of tickmark
% e.g. [1.5 1.5 1.5;.5 .5 .5] for trial 1, etc.

% collapse spike windows into a single matrix, get cluster ids back and then 
% change color according to ID

disp('Generating figures...');
disp(['Will save to directory:  ' savedir]);

% scale pixels by time
% now we have clusterid and the trial for all spikes, we need to color them appropriately

[path,name,ext]=fileparts(fullfile(savedir));
savefilename=[ name '_sua_freqrange_' num2str(freq_range) '_electrode_' ];

savedir=fullfile(savedir,'sua');

if ~exist(savedir,'dir')
	mkdir(fullfile(savedir));
end





for i=1:length(channels)

	% [time;time] [trial+.5;trial-.5]

	

	disp(['Channel ' num2str(channels(i))]);

	if sort
		if auto_clust
			[clusterid clustertrial clusterisi clusterwindows]=ephys_spike_cluster_auto(spikewindows{i},spiketimes{i},'sr',SR);
		else
			[clusterid clustertrial clusterisi clusterwindows]=ephys_spike_clustergui_tetrode(spikewindows{i},spiketimes{i},'sr',SR);
		end

		if isempty(clusterid)
			disp('No labels returned, bailing...');
			return;
		end

	else
		clusterid=[];
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% IFR, SMOOTH RATE %%%%%%%%%%%%%%%%%%%


	uniq_clusters=unique(clusterid);

	if ~isempty(clusterid)
	
		% cycle through each cluster id

		for j=1:length(uniq_clusters)

			IFR{i}{j}=zeros(ntrials,samples);
			spike_vec{j}=[];
			trial_vec{j}=[];

			for k=1:ntrials

				% grab the spike labels for this electrode and trial

				currids=clusterid(find(clustertrial==k));

				% get the spike times (in seconds)

				currspikes=spiketimes{i}{k};
			
				% now grab the spikes for THIS CLUSTER only

				clusterspikes=currspikes(find(currids==uniq_clusters(j)));
				
				% IFR will use the current sampling rate
				
				IFR{i}{j}(k,:)=ephys_ifr(round(clusterspikes),samples,SR);

				clusterspikes=clusterspikes./SR;
				
				tmp=repmat(clusterspikes,2,1);
				spike_vec{j}=[spike_vec{j} tmp];

				tmp=[k+.3;k-.3];
				tmp=repmat(tmp,1,length(clusterspikes));
				trial_vec{j}=[trial_vec{j} tmp];

				% for storage resort the spike vector by cluster

				clust_spike_vec{i}{j}{k}=clusterspikes;

				binned_spikes=histc(clusterspikes,binedges);
				smooth_spikes.density{j}(k,:)=conv(binned_spikes,kernel,'same');


			end
		end
	else

		clust_spike_vec=[];
		spike_vec{1}=[];
		trial_vec{1}=[];

		for j=1:ntrials

			% now for each trial, we need to label the spikes appropriately,
			% first sort by id

			% find command is already sorted, ascending

			tmp=repmat(spiketimes{i}{j},2,1)./SR;
			spike_vec{1}=[spike_vec{1} tmp];

			tmp=[j+.3;j-.3];
			tmp=repmat(tmp,1,length(spiketimes{i}{j}));
			trial_vec{1}=[trial_vec{1} tmp];

			binned_spikes=histc(spiketimes{i}{j},binedges);
			smooth_spikes{1}(j,:)=conv(binned_spikes,kernel,'same');

		end

	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


	if ~isempty(uniq_clusters)
		nplots=length(uniq_clusters);
	else
		nplots=1;
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% STATS PLOTTING %%%%%%%%%%%%%%%%%%%%%


	% plot the spike stats

	savefilename_stats=[ savefilename num2str(channels(i))...
		       	'_stats_cluster_'];

	% delete all old files

	delete(fullfile(savedir,[ '*_sua_freqrange*electrode_' num2str(channels(i)) '*']));
	
	if ~isempty(clusterid)

		% delete old stats figures if they exist

		delete(fullfile(savedir,[savefilename_stats '*.png']));
		delete(fullfile(savedir,[savefilename_stats '*.eps']));

		for j=1:length(uniq_clusters)

			% estimate SNR from 6*std of spikeless trace

			ntrials=length(spikeless{i});

			for k=1:ntrials
				noise_p2p(k)=6*std(spikeless{i}{k});
			end

			noise_p2p(isnan(noise_p2p))=[];
			mean_noise_p2p=mean(noise_p2p);
			stats_fig=ephys_visual_spikestats(clusterwindows{uniq_clusters(j)},clusterisi{uniq_clusters(j)},...
				'noise_p2p',mean_noise_p2p,'sr',25e3,'spike_sr',50e3);

			multi_fig_save(stats_fig,savedir,...
				[ savefilename_stats num2str(uniq_clusters(j))],'eps,png');
			close([stats_fig]);

		end
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RASTER PLOTTING %%%%%%%%%%%%%%%%%%%%%

	savefilename_raster=[ savefilename num2str(channels(i))...
		       	'_raster_cluster_'];

	ax=[];

	% delete any old single trial plots

	if exist(fullfile(savedir,'singletrials',['ch' num2str(channels(i)) ]),'dir')
		rmdir(fullfile(savedir,'singletrials',['ch' num2str(channels(i)) ]),'s');
	end

	for j=1:length(uniq_clusters)

		raster_fig=figure('visible','off','Units','Pixels','Position',[0 0 600 800]);
		ax(1)=subaxis(6,1,1,1,1,2,'margin',.1,'spacingvert',0);
		imagesc(HISTOGRAM.t,HISTOGRAM.f(startidx:stopidx),HISTOGRAM.imask(startidx:stopidx,:));
		colormap(hist_colors);
		freezeColors;
		set(gca,'ydir','normal','tickdir','out','xtick',[],'ytick',[min_f max_f],'linewidth',1.5,'ticklength',[.025 .025],...
			'FontSize',11,'FontName','Helvetica');
		ylim([min_f max_f]);
		box off;
		title({[ figtitle ];[ 'Channel ' num2str(channels(i))]},'FontSize',18,'FontName','Helvetica');

		ax(2)=subaxis(6,1,1,3,1,1,'spacingvert',0,'margin',0.1,'paddingbottom',.025);
		plot(TIME,HISTOGRAM.mean_osc,'-k');
		ylabel('Osc.','FontSize',13,'FontName','Helvetica');
		axis tight;
		set(gca,'tickdir','out','xtick',[],'ytick',[]);

		ax(3)=subaxis(6,1,1,4,1,3,'spacingvert',0.025,'margin',0.1,'paddingbottom',0);

		plot(spike_vec{j},trial_vec{j},'-','color','k');	
		axis([0 TIME(end) 0 ntrials+1]);
		xlabel('Time (in s)','FontSize',13,'FontName','Helvetica');
		ylabel('Trial','FontSize',13,'FontName','Helvetica');
		box off
		set(gca,'tickdir','out','linewidth',1.5,'ticklength',[.025 .025],'FontSize',11,'FontName','Helvetica','ydir','rev');

		linkaxes(ax,'x');

		if ~isempty(figtitle)
			name=figtitle;
		end

		set(raster_fig,'PaperPositionMode','auto');
		multi_fig_save(raster_fig,savedir,...
			[ savefilename_raster num2str(uniq_clusters(j)) ],'eps,png');
		close([raster_fig]);

		% generate single trial plots, simply the spikes marked along with the filtered traces

		randpopulation=randperm(ntrials);

		if singletrials<=ntrials
			randtrials=randpopulation(1:singletrials);
		else
			randtrials=randpopulation(1:ntrials);
			singletrials=ntrials;
		end

		singletrialdir=fullfile(savedir,'singletrials',['ch' num2str(channels(i)) ],[ 'clust' num2str(j)]);

		if ~isempty(clusterid)

			if ~exist(fullfile(singletrialdir),'dir');
				mkdir(singletrialdir);
			end

			disp('Plotting single trials');

			for k=1:singletrials

				singletrialfig=figure('Visible','off');
				idx=randtrials(k);
				currspikes=clust_spike_vec{i}{j}{idx};
				plot(TIME,proc_data(:,idx,i),'k-');
				hold on
				plot(currspikes,proc_data(round(currspikes*SR),idx,i),'r*');
				set(gca,'TickDir','out','TickLength',[.02 .02],'FontName','Helvetica','FontSize',11);
				xlabel('Time (s)','FontName','Helvetica','FontSize',13,'interpreter','latex');
				ylabel('Voltage ($\mu$V)','FontName','Helvetica','FontSize',13,'interpreter','latex');
				box off
				axis tight
				multi_fig_save(singletrialfig,singletrialdir,['trial' num2str(idx)],'eps,png');

				close([singletrialfig]);

			end
		end
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

if sort 
	save(fullfile(savedir,['sua_channels ' num2str(channels) '.mat']),'clusterid','clustertrial','clusterisi','clusterwindows',...
		'threshold','CHANNELS','channels','TIME','proc_data','freq_range','clust_spike_vec','smooth_spikes','spikeless','IFR','subtrials');
else
	save(fullfile(savedir,['sua_channels ' num2str(channels) '.mat']),'threshold','CHANNELS','channels','TIME','proc_data','freq_range');
end

