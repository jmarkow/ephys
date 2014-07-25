function [LFP_RASTER TIME LABEL HISTOGRAM]=ephys_lfp_amp(EPHYS,HISTOGRAM,varargin)
%generates song-aligned LFP amplitude rasters
%
%	[LFP_RASTER TIME LABEL HISTOGRAM]=ephys_lfp_amp(EPHYS.data,HISTOGRAM,EPHYS.labels,varargin)
%
%	EPHYS.data
%	aligned ephys data generated by ephys_cluster or in extracted_data/aggregated_data.mat
%	(should be the variable ephys_data), the data should be a matrix of doubles that is 
%	samples x trials x channels
%
%	HISTOGRAM
%	contour histogram returned by ephys_visual_histogram.m (or loaded from histogram.mat)
%
%	EPHYS.labels
%	channel labels (i.e. the channel that corresponds to a given element in the cell array
%	ephys_data) from ephys_cluster.m or extracted_data/aggregated_data.mat
%
%	the following may be specified as parameter/value pairs:
%
%		car_exclude
%		electrodes to exclude from noise estimate
%
%		fs
%		data sampling rate (default: 25e3)
%
%		noise
%		noise rejection method ('car' for common average 'nn' for nearest neighbor, or 'none',
%		default: 'none')
%
%		freq_range
%		vector with two elements to specify the frequency range (one element specifies low pass)
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
%		mua_colors
%		colormap (string) for raster data (default: hot)
%
%		hist_colors
%		colormap(string) for histogram (default: jet)
%
%		medfilt_scale
%		timescale for median filter (in ms, leave blank to skip, default: 1.5)
%
%		hampel
%		a simple hampel or median filter is used to throw out noise trials for display
%		(all trials are stored regardless), basically the number specifies how many MADs
%		the rms on a particular trial must be from the median computed across trials (30 trial window)
%to throw out 3-4 is standard, with lower being more aggresive (leave blank to skip, default: 3)
%
% see also ephys_visual_sua.m,ephys_visual_lfp_tf.m,ephys_visual_mua.m

if nargin<3
	error('ephysPipeline:lfpampvis:notenoughparams','Need 3 arguments to continue, see documentation');
end


nparams=length(varargin);

if mod(nparams,2)>0
	error('ephysPipeline:argChk','Parameters must be specified as parameter/value pairs!');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETER COLLECTION %%%%%%%%%%%%%%%%%

noise='none'; 
car_exclude=[];
savedir=pwd;
min_f=1;
max_f=10e3;
hist_colors='jet';
mua_colors='hot';
mua_colors_phase='hsv';
figtitle='';
freq_range=[10 80]; % frequency range for filtering
proc_fs=1.25e3;
channels=EPHYS.labels;
medfilt_scale=1.5; % median filter scale (in ms)
hampel=3;

%%%

for i=1:2:nparams
	switch lower(varargin{i})
		case 'noise'
			noise=varargin{i+1};
		case 'savedir'
			savedir=varargin{i+1};
		case 'min_f'
			min_f=varargin{i+1};
		case 'max_f'
			max_f=varargin{i+1};
		case 'mua_colors'
			mua_colors=varargin{i+1};
		case 'hist_colors'
			hist_colors=varargin{i+1};
		case 'car_exclude'
			car_exclude=varargin{i+1};
		case 'figtitle'
			figtitle=varargin{i+1};
		case 'freq_range'
			freq_range=varargin{i+1};
		case 'downsampling'
			downsampling=varargin{i+1};
		case 'channels'
			channels=varargin{i+1};
		case 'medfilt_scale'
			medfilt_scale=varargin{i+1};
		case 'hampel'
			hampel=varargin{i+1};
		case 'proc_fs'
			proc_fs=varargin{i+1};
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% denoise and condition signal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SIGNAL CONDITIONING %%%%%%%%%%%%%%%%

fs=EPHYS.fs;


proc_data=ephys_denoise_signal(EPHYS.data,EPHYS.labels,channels,'method',noise,'car_exclude',car_exclude);
clear EPHYS.data;
proc_data=single(ephys_condition_signal(proc_data,'l','freq_range',freq_range,'medfilt_scale',medfilt_scale));

downfact=fs/proc_fs;

if mod(downfact,1)>0
	error('ephysPipeline:spectcoherence:downsamplenotinteger','Need to downsample by integer');
end

disp(['Downsampling to ' num2str(proc_fs) ]);
proc_data=downsample(proc_data,downfact);

[nsamples,ntrials,nchannels]=size(proc_data);
LFP_RASTER.t=[1:nsamples]./proc_fs;
LFP_RASTER.trials=[1:ntrials];

% are we downsampling
% downsampling is done by skipping samples, straightforward

LFP_RASTER.image_amp=zeros(ntrials,nsamples,nchannels);
LFP_RASTER.image_phase=zeros(size(LFP_RASTER.image_amp));

for i=1:length(channels)
	for j=1:ntrials
		currdata=hilbert(proc_data(:,j,i));
		LFP_RASTER.image_phase(j,:,i)=mod(unwrap(angle(currdata)),2*pi);
		LFP_RASTER.image_amp(j,:,i)=proc_data(:,j,i);
	end
end

clear proc_data;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTTING CODE %%%%%%%%%%%%%%%%%%%%%%

disp('Generating figures...');
disp(['Will save to directory:  ' savedir]);

% set up directories for output

if length(freq_range)==2
	subdir=[ num2str(freq_range(1)) '_' num2str(freq_range(2))];
else
	subdir=num2str(freq_range(1));
end


[path,name,ext]=fileparts(savedir);

savedir=fullfile(savedir,'lfp_amp');

if ~exist(fullfile(savedir,subdir,'amp'),'dir') || ~exist(fullfile(savedir,subdir,'phase'),'dir')
	mkdir(fullfile(savedir,subdir,'amp'));
	mkdir(fullfile(savedir,subdir,'phase'));
end

savefilename=[ name '_lfp_freqrange_' num2str(freq_range) '_electrode_'];

% delete any old rasters

delete(fullfile(savedir,subdir,[savefilename '*.png']));
delete(fullfile(savedir,subdir,[savefilename '*.eps']));

delete(fullfile(savedir,subdir,'amp',[savefilename '*.png']));
delete(fullfile(savedir,subdir,'amp',[savefilename '*.eps']));

delete(fullfile(savedir,subdir,'phase',[savefilename '*.png']));
delete(fullfile(savedir,subdir,'phase',[savefilename '*.eps']));
goodtrials=[];

for i=1:length(channels)

	ax=[];
	raster_fig=figure('visible','off','Units','Pixels','Position',[0 0 700 1e3]);

	reject=[];
	if ~isempty(hampel)	
		reject=hampel_filter(LFP_RASTER.image_amp(:,:,i)','hampel_factor',hampel);
	end
	goodtrials=setdiff(LFP_RASTER.trials,reject);

	PLOTLFP=LFP_RASTER;
	PLOTLFP.trials=PLOTLFP.trials(goodtrials);
	PLOTLFP.image=PLOTLFP.image_amp(goodtrials,:,i);

	if ~isempty(HISTOGRAM)
		multi_unit_raster(HISTOGRAM,PLOTLFP,'fs',fs,...
			'fig_num',raster_fig,'fig_title',{[figtitle];[ 'Channel ' num2str(channels(i))]},...
			'min_f',min_f,'max_f',max_f,'raster_colors',mua_colors,'hist_colors',hist_colors);
	else
		imagesc(PLOTLFP.t,PLOTLFP.trials,PLOTLFP.image);
		colormap(mua_colors);
		axis xy;
		xlabel('Time (in s)','FontSize',13,'FontName','Helvetica');
		ylabel('Trial','FontSize',13,'FontName','Helvetica');
		box off
		set(gca,'tickdir','out','linewidth',1.5,'ticklength',[.025 .025],...
			'FontSize',11,'FontName','Helvetica');
	end

	set(raster_fig,'PaperPositionMode','auto')
	multi_fig_save(raster_fig,fullfile(savedir,subdir,'amp'),...
		[ savefilename num2str(channels(i)) '_amp' ],'eps,png');
	
	close([raster_fig]);

	raster_fig=figure('visible','off','Units','Pixels','Position',[0 0 700 1e3]);

	PLOTLFP.image=PLOTLFP.image_phase(goodtrials,:,i);

	if ~isempty(HISTOGRAM)
		multi_unit_raster(HISTOGRAM,PLOTLFP,'fs',fs,...
			'fig_num',raster_fig,'fig_title',{[figtitle];[ 'Channel ' num2str(channels(i))]},...
			'min_f',min_f,'max_f',max_f,'raster_colors',mua_colors_phase,'hist_colors',hist_colors,...
			'show_colorbar',1,'scalelabel','Phase');
	else
		imagesc(PLOTLFP.t,PLOTLFP.trials,PLOTLFP.image);
		colormap(mua_colors);
		axis xy;
		xlabel('Time (in s)','FontSize',13,'FontName','Helvetica');
		ylabel('Trial','FontSize',13,'FontName','Helvetica');
		box off
		set(gca,'tickdir','out','linewidth',1.5,'ticklength',[.025 .025],...
			'FontSize',11,'FontName','Helvetica');
	end

	set(raster_fig,'PaperPositionMode','auto')
	multi_fig_save(raster_fig,fullfile(savedir,subdir,'phase'),...
		[ savefilename num2str(channels(i)) '_phase' ],'eps,png');

	close([raster_fig]);

end

save(fullfile(savedir,subdir,['lfp_amp_freqrange_' num2str(freq_range) '.mat']),...
	'LFP_RASTER','channels','freq_range','goodtrials','reject');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

