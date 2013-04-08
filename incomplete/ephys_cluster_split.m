function ephys_cluster_split(CHANNEL,CLUSTERNUM,HISTOGRAM,varargin)
%
%
%
%
%
%
%
nparams=length(varargin);

if mod(nparams,2)>0
	error('ephysPipeline:argChk','Parameters must be specified as parameter/value pairs!');
end

fs=25e3;
noise='none'; % none, nn for nearest neighbor, or car for common average
car_exclude=[];
savedir=pwd;
min_f=1; % min frequency to show for song histogram
max_f=10e3; % max frequency
hist_colors='jet'; % colormap for histogram
figtitle='';
freq_range=[500]; % bandpassing <10e3 distorted results, reasoning that >800 Hz is fine for spikes < 1ms long
filt_type='high'; % high,low or bandpass
filt_order=3;
spikesort=1; % do we want to sort?
auto_clust=0; % 0 for manual cluster cutting (GUI), 1 for automated clustering
tetrode_channels=[];
sigma_t=4; % multiple of noise estimate for spike threshold (generally 3-4, using Quiroga's method)
jitter=4; % max jitter in samples for spike re-alignment (4 is reasonable
singletrials=0; % number of random single trials to plot per cluster
subtrials=[];
align='min'; % how to align spike waveforms can be min, max or com for center-of-mass
interpolate_fs=200e3; % 200 has worked best here
smooth_rate=1e3;
sigma=.0025;
wavelet_method='ks'; % ks or bimodal have been sucessful
wavelet_mpca=0; % mpca seems to help...
wavelet_coeffs=20; % 10 has worked well (3-5 for mpca)
clust_choice='bic'; % with merging, BIC works just fine...
red_cutoff=.6;
outlier_detect=1;
outlier_cutoff=150;
nfeatures=10;
merge=.6; % exp(-d), lower to merge more aggressively (see Hennig "Methods for merging...", 2009)
savename=''; % add if doing multiple manual sorts, will append a name to the filename
isi_cutoff=.05; % percentage of ISI values <.001
lratio_cutoff=.2; % l-ratio cutoff from Reddish
isod_cutoff=20; % isolation distance defined by Harris et al.
snr_cutoff=1.1; % SNR definition from Ludwig et al. 2009 (J. Neurophys)
spike_window=[.0004 .0004];
trial_timestamps=[];
auto_clust=1;

for i=1:2:nparams
	switch lower(varargin{i})
		case 'fs'
			fs=varargin{i+1};
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
		case 'outlier_cutoff'
			outlier_cutoff=varargin{i+1};
		case 'freq_range'
			freq_range=varargin{i+1};
		case 'auto_clust'
			auto_clust=varargin{i+1};
		case 'clust_choice'
			clust_choice=varargin{i+1};
		case 'wavelet_coeffs'
			wavelet_coeffs=varargin{i+1};
		case 'wavelet_mpca'
			wavelet_mpca=varargin{i+1};
		case 'wavelet_method'
			wavelet_method=varargin{i+1};
		case 'interpolate_fs'
			interpolate_fs=varargin{i+1};
		case 'outlier_detect'
			outlier_detect=varargin{i+1};
		case 'red_cutoff'
			red_cutoff=varargin{i+1};
		case 'isi_cutoff'
			isi_cutoff=varargin{i+1};
		case 'isod_cutoff'
			isod_cutoff=varargin{i+1};
		case 'lratio_cutoff'
			lratio_cutoff=varargin{i+1};
		case 'snr_cutoff'
			snr_cutoff=varargin{i+1};
		case 'nfeatures'
			nfeatures=varargin{i+1};
		case 'merge'
			merge=varargin{i+1};
	end
end

if ~isempty(HISTOGRAM)
	startidx=max([find(HISTOGRAM.f<=min_f)]);

	if isempty(startidx)
		startidx=1;
	end

	stopidx=min([find(HISTOGRAM.f>=max_f)]);

	if isempty(stopidx)
		stopidx=length(HISTOGRAM.f);
	end
end

% start by checking the the right file and cluster
% assume we're already in the correct directory

load(fullfile(pwd,['sua_channels ' num2str(CHANNEL) '.mat']),'threshold','CHANNELS','channels',...
	'TIME','proc_data','freq_range','subtrials','spikeless','trial_timestamps',...
	'clust_spike_vec','clusterwindows','clusterpoints','proc_data','clusterid','clustertrial',...
	'clusterstats','spikeless','clusterisi');

name='ephys';
savefilename=[ name '_sua_freqrange_' num2str(freq_range) '_electrode_' ];
[nsamples,ntrials,nchannels]=size(proc_data);

savedir=pwd;

% retain data from old clustering, delete cluster to split and append results

splittimes=clust_spike_vec{1}{CLUSTERNUM};

for i=1:length(splittimes)
	splittimes{i}=round(splittimes{i}.*fs);
end
	
splitclustertrials=clustertrial(find(clusterid==CLUSTERNUM));
uniqtrials=unique(splitclustertrials);

for i=1:length(uniqtrials)
	splitwindows{i}=clusterwindows{CLUSTERNUM}(:,splitclustertrials==uniqtrials(i));
end

% arrange windows into trial cell array

if auto_clust
	[newclusterid newclustertrial newclusterisi newclusterwindows newclusterpoints newclusterstats newoutliers]=...
		ephys_spike_cluster_auto(splitwindows,splittimes,...
		'fs',fs,'wavelet_method',wavelet_method,'wavelet_mpca',wavelet_mpca,...
		'clust_choice',clust_choice,'maxcoeffs',wavelet_coeffs,'outlier_detect',outlier_detect,...
		'red_cutoff',red_cutoff,'nfeatures',nfeatures,'merge',merge,'outlier_cutoff',outlier_cutoff,...
		'dim_reduce','ks');
end





% reshape clust_spike_vec and all

clusterwindows(CLUSTERNUM)=[];
clusterpoints(CLUSTERNUM)=[];
clusterisi(CLUSTERNUM)=[];

clusterstats.lratio(CLUSTERNUM)=[];
clusterstats.isod(CLUSTERNUM)=[];
clusterstats.snr(CLUSTERNUM)=[];

clust_spike_vec{1}(CLUSTERNUM)=[];
old_uniq_clusters=unique(clusterid);
clusterid(clusterid==CLUSTERNUM)=newclusterid+length(old_uniq_clusters);
newclusters=1:length(newclusterwindows);

noldclust=length(clust_spike_vec{1});

for i=1:length(newclusters)

	clusterwindows{end+1}=newclusterwindows{i};
	clusterpoints{end+1}=newclusterpoints{i};
	clusterisi{end+1}=newclusterisi{i};
	clusterstats.lratio(end+1)=newclusterstats.lratio(i);
	clusterstats.isod(end+1)=newclusterstats.isod(i);

	for j=1:ntrials
		currids=newclusterid(find(newclustertrial==j));
		currspikes=splittimes{j};
		clusterspikes=currspikes(find(currids==newclusters(i)));
		clust_spike_vec{1}{noldclust+i}{j}=clusterspikes;
	end
end

uniq_clusters=1:length(clusterwindows);

% kernel density estimate of smooth firing rate

kernedges=[-3*sigma:1/smooth_rate:3*sigma];
binedges=[0:(1/smooth_rate):nsamples/fs];

smooth_spikes.time=binedges;

if exist('normpdf')>0
	kernel=normpdf(kernedges,0,sigma);
else
	kernel=(1/(sigma*sqrt(2*pi)))*exp((-(kernedges-0).^2)./(2*sigma^2));
end


% cycle through each cluster id

for j=1:length(uniq_clusters)

	IFR{1}{j}=zeros(ntrials,nsamples);
	spike_vec{j}=[];
	trial_vec{j}=[];

	for k=1:ntrials

		clusterspikes=clust_spike_vec{1}{j}{k};

		IFR{1}{j}(k,:)=ephys_ifr(round(clusterspikes),nsamples,fs);

		clusterspikes=clusterspikes./fs;

		tmp=repmat(clusterspikes,2,1);
		spike_vec{j}=[spike_vec{j} tmp];

		tmp=[k+.3;k-.3];
		tmp=repmat(tmp,1,length(clusterspikes));
		trial_vec{j}=[trial_vec{j} tmp];

		% for storage resort the spike vector by cluster

		binned_spikes=histc(clusterspikes,binedges);
		smooth_spikes.density{j}(k,:)=conv(binned_spikes,kernel,'same');

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

savefilename_stats=[ savefilename num2str(CHANNEL)...
	'_stats_cluster_'];

% delete all old files

delete(fullfile(savedir,[ '*_sua_freqrange*electrode_' num2str(CHANNEL) '*']));

% delete old candidate files

delete(fullfile(savedir,[ 'candidate_unit_ch' num2str(CHANNEL) '*']));
delete(fullfile(savedir,[ '.sua_channels ' num2str(CHANNEL) '*']));

if ~isempty(clusterid)

	% delete old stats figures if they exist

	delete(fullfile(savedir,[savefilename_stats '*.png']));
	delete(fullfile(savedir,[savefilename_stats '*.eps']));

	for j=1:nplots

		% estimate SNR from 6*std of spikeless trace

		ntrials=length(spikeless{1});

		for k=1:ntrials
			noise_p2p(k)=6*std(spikeless{1}{k});
		end

		noise_p2p(isnan(noise_p2p))=[];
		mean_noise_p2p=mean(noise_p2p);
		mean_wave=mean(clusterwindows{uniq_clusters(j)},2);
		max_wave=max(mean_wave);
		min_wave=min(mean_wave);

		clusterstats.snr(uniq_clusters(j))=[ abs(max_wave-min_wave)/mean_noise_p2p ];

		stats_fig=figure('Visible','off');

		note=[];

		note=['L-ratio ' sprintf('%.2f',clusterstats.lratio(uniq_clusters(j))) ...
			' IsoD ' sprintf('%.2f',clusterstats.isod(uniq_clusters(j))) ];

		stats_fig=ephys_visual_spikestats(clusterwindows{uniq_clusters(j)},clusterisi{uniq_clusters(j)},...
			'noise_p2p',mean_noise_p2p,'fs',fs,'spike_fs',interpolate_fs,'fig_num',stats_fig,'note',note);

		% TODO:  function for cluster stats (Fisher test, etc.)

		set(stats_fig,'Position',[0 0 450 600]);
		set(stats_fig,'PaperPositionMode','auto');

		% label candidate units if they meet our criteria
		% isi intervals < absolute refractory period

		isi_violations=sum((clusterisi{uniq_clusters(j)}./fs)<.001);
		isi_violations=isi_violations/length(clusterisi{uniq_clusters(j)});

		if clusterstats.snr(uniq_clusters(j))>=snr_cutoff && ...
				(isnan(clusterstats.lratio(uniq_clusters(j))) || ...
				clusterstats.lratio(uniq_clusters(j))<=lratio_cutoff) &&  ...
				isi_violations<isi_cutoff 
			if isnan(clusterstats.isod(uniq_clusters(j))) || clusterstats.isod(uniq_clusters(j))>=isod_cutoff
				fid=fopen(fullfile(savedir,['candidate_unit_ch' num2str(CHANNEL) ... 
					'_cl' num2str(uniq_clusters(j))]),'w');
				fclose(fid);
			end
		end

		multi_fig_save(stats_fig,savedir,...
			[ savefilename_stats num2str(uniq_clusters(j))],'png','res',200);
		close([stats_fig]);

	end

	stats_fig=figure('Visible','off','renderer','painters');

	stats=[];

	if auto_clust
		stats=clusterstats;
	end

	stats_fig=ephys_visual_cluststats(clusterwindows,clust_spike_vec{1},...
		'spike_fs',interpolate_fs,'fig_num',stats_fig,'stats',stats);

	set(stats_fig,'Position',[0 0 250+500*nplots 250+500*nplots]);
	set(stats_fig,'PaperPositionMode','auto');

	multi_fig_save(stats_fig,savedir,...
		[ savefilename_stats 'clstats' ],'eps,png','res',150,'renderer','painters');

	close([stats_fig]);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RASTER PLOTTING %%%%%%%%%%%%%%%%%%%%%

savefilename_raster=[ savefilename num2str(CHANNEL)...
	'_raster_cluster_'];

ax=[];

% delete any old single trial plots

if exist(fullfile(savedir,'singletrials',['ch' num2str(CHANNEL) ]),'dir')
	rmdir(fullfile(savedir,'singletrials',['ch' num2str(CHANNEL) ]),'s');
end

for j=1:length(uniq_clusters)

	raster_fig=figure('visible','off','Units','Pixels','Position',[0 0 600 800]);

	if ~isempty(HISTOGRAM)
		ax(1)=subaxis(6,1,1,1,1,2,'margin',.1,'spacingvert',0);
		imagesc(HISTOGRAM.t,HISTOGRAM.f(startidx:stopidx),HISTOGRAM.imask(startidx:stopidx,:));
		colormap(hist_colors);
		freezeColors;
		set(gca,'ydir','normal','tickdir','out','xtick',[],'ytick',[min_f max_f],'linewidth',1.5,'ticklength',[.025 .025],...
			'FontSize',11,'FontName','Helvetica');
		ylim([min_f max_f]);
		box off;
		title({[ figtitle ];[ 'Channel ' num2str(CHANNEL)]},'FontSize',18,'FontName','Helvetica');

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

	else
		plot(spike_vec{j},trial_vec{j},'-','color','k');
		axis([0 TIME(end) 0 ntrials+1]);

		xlabel('Time (in s)','FontSize',13,'FontName','Helvetica');
		ylabel('Trial','FontSize',13,'FontName','Helvetica');
		box off
		set(gca,'tickdir','out','linewidth',1.5,'ticklength',[.025 .025],'FontSize',11,'FontName','Helvetica','ydir','rev');

	end

	if ~isempty(figtitle)
		name=figtitle;
	end

	set(raster_fig,'PaperPositionMode','auto');

	savedir
	savefilename_raster

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

	singletrialdir=fullfile(savedir,'singletrials',['ch' num2str(CHANNEL) ],[ 'clust' num2str(j)]);

	if ~isempty(clusterid)

		if ~exist(fullfile(singletrialdir),'dir');
			mkdir(singletrialdir);
		end

		disp('Plotting single trials');

		for k=1:singletrials

			singletrialfig=figure('Visible','off');
			idx=randtrials(k);
			currspikes=clust_spike_vec{1}{j}{idx};
			plot(TIME,proc_data(:,idx,1),'k-');
			hold on
			plot(currspikes,proc_data(round(currspikes*fs),idx,i),'r*');
			set(gca,'TickDir','out','TickLength',[.02 .02],'FontName','Helvetica','FontSize',11);
			xlabel('Time (s)','FontName','Helvetica','FontSize',13,'interpreter','latex');
			ylabel('Voltage ($\mu$V)','FontName','Helvetica','FontSize',13,'interpreter','latex');
			box off
			axis tight
			multi_fig_save(singletrialfig,singletrialdir,['trial' num2str(idx)],'png','res',150);

			close([singletrialfig]);

		end
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%save(fullfile(savedir,['sua_channels ' num2str(channels) '.mat']),...
%	'clusterid','clustertrial','clusterisi','clusterwindows','fs','interpolate_fs',...
%	'threshold','CHANNELS','channels','TIME','proc_data','freq_range','clust_spike_vec','smooth_spikes',...
%	'spikeless','IFR','subtrials','clusterpoints','clusterstats','trial_timestamps');
