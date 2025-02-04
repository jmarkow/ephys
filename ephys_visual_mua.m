function [MUA TIME LABEL HISTOGRAM]=ephys_visual_mua(EPHYS,HISTOGRAM,varargin)
%generates song-aligned mult-unit rasters
%
%	[MUA TIME LABEL HISTOGRAM]=ephys_visual_mua(EPHYS.data,HISTOGRAM,EPHYS.labels,varargin)
%
%	EPHYS.data
%	sound-aligned voltage traces from extracted_data.mat (should be the variable ephys_data)
%
%	HISTOGRAM
%	contour histogram returned by ephys_visual_histogram.m (or loaded from histogram.mat)
%
%	EPHYS.labels
%	channel labels (i.e. the channel that corresponds to a given element in the cell array
%	ephys_data) from extracted_data.mat%
%	the following may be specified as parameter/value pairs:
%
%		exclude
%		electrodes to exclude from noise estimate
%
%		smooth_window
%		smoothing_window for multi-unit in seconds (default .005)
%
%		fs
%		data sampling rate (default: 25e3)
%
%		noise
%		noise rejection method (default: 'car', or common-average-rejection)
%
%		filtering
%		define as a two-element vector with lower and upper corner frequencies 
%		to filtering multi-unit traces (default: none)
%
%		savedir
%		directory to store results (default: pwd)
%
%		min_f
%		lowermost frequency to display for contour histogram
%
%		max_f
%		uppermost frequency to display for contour histogram
%
%		colors
%		colormap for contour histogram and multi-unit data
%
% see also ephys_visual_sua.m,ephys_visual_lfp_amp.m,ephys_visual_lfp_tf.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETER COLLECTION %%%%%%%%%%%%%%%%%

if nargin<2
	HISTOGRAM=[];
end

if nargin<1
	error('ephysPipeline:muavis:notenoughparams','Need 1 argument to continue, see documentation');
end

nparams=length(varargin);

if mod(nparams,2)>0
	error('ephysPipeline:argChk','Parameters must be specified as parameter/value pairs!');
end

%%%

winsigma=.0025; % smoothing window in secs
noise='none'; % common-average reference for noise removal, none to skip digital
	      % re-referencing
car_exclude=[];
savedir=pwd;
min_f=1;
max_f=10e3;
hist_colors='jet';
mua_colors='hot';
figtitle='';
freq_range=[500 4.5e3]; % frequency range for filtering
proc_fs=10e3;
downsampling=2;
channels=EPHYS.labels;
hampel=6;
attenuation=40;
ripple=.2;
car_trim=40;

for i=1:2:nparams
	switch lower(varargin{i})
		case 'winsigma'
			smooth_window=varargin{i+1};
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
		case 'hampel'
			hampel=varargin{i+1};
		case 'attenuation'
			attenuation=varargin{i+1};
		case 'ripple'
			ripple=varargin{i+1};
		case 'proc_fs'
			proc_fs=varargin{i+1};
		case 'car_trim'
			car_trim=varargin{i+1};
	end
end

% perhaps add a median filter here to remove spikes from the multi-unit trace

% intan nearest neighbor mapping

fs=EPHYS.fs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SIGNAL CONDITIONING %%%%%%%%%%%%%%%%

downfact=fs/proc_fs;

if mod(downfact,1)>0
	error('ephysPipeline:spectcoherence:downsamplenotinteger','Need to downsample by integer');
end

proc_data=ephys_denoise_signal(EPHYS.data,EPHYS.labels,channels,'method',noise,'car_exclude',car_exclude,'car_trim',car_trim);
clear EPHYS.data;

disp('Anti-alias filtering and downsampling');

proc_data=ephys_condition_signal(proc_data,'l','freq_range',[proc_fs/2],'filt_order',2,'filt_name','b',...
		'medfilt_scale',[],'fs',fs,'notch',0); 

disp(['Downsampling to ' num2str(proc_fs) ]);
proc_data=downsample(proc_data,downfact);
proc_data=single(ephys_condition_signal(proc_data,'m','freq_range',freq_range,'winsigma',winsigma,'fs',proc_fs));

[nsamples,ntrials,nchannels]=size(proc_data);
TIME=[1:nsamples]./proc_fs;

% are we downsampling

if ~isempty(downsampling)
	disp(['Downsampling by factor of ' num2str(downsampling)]);
	MUA.t=downsample(TIME,downsampling);
	MUA.image=zeros(ntrials,length(MUA.t),length(channels),'single');

	for i=1:nchannels
		MUA.image(:,:,i)=downsample(proc_data(:,:,i),downsampling)';
	end

else
	MUA.t=TIME;
	MUA.image=zeros(ntrials,length(MUA.t),length(channels),'single');

	for i=1:nchannels
		MUA.image(:,:,i)=proc_data(:,:,i)';
	end
end

MUA.channels=channels;
MUA.trials=1:ntrials;
clear proc_data;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Generating figures...');
disp(['Will save to directory:  ' savedir]);

% scale pixels by time

[path,name,ext]=fileparts(savedir);

savedir=fullfile(savedir,'mua');

if ~exist(savedir,'dir')
	mkdir(savedir);
end

savefilename=[ name '_mua_freqrange_' num2str(freq_range) '_electrode_' ];

% delete any old rasters

delete(fullfile(savedir,[savefilename '*.png']));
delete(fullfile(savedir,[savefilename '*.eps']));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTTING CODE %%%%%%%%%%%%%%%%%%%%%%

for i=1:length(channels)

	ax=[];
	raster_fig=figure('visible','off','Units','Pixels','Position',[0 0 700 1e3]);
	
	% hampel filter slides across trials and checks for extreme deviations
	% where the rms is > hampel*MAD

	reject=[];
	if ~isempty(hampel)
		reject=hampel_filter(MUA.image(:,:,i)','hampel_factor',hampel);
	end

	goodtrials=setdiff(MUA.trials,reject);

	PLOTMUA=MUA;
	PLOTMUA.trials=PLOTMUA.trials(goodtrials);
	PLOTMUA.image=PLOTMUA.image(goodtrials,:,i);

	if ~isempty(HISTOGRAM)
		multi_unit_raster(HISTOGRAM,PLOTMUA,'fs',fs,...
			'fig_num',raster_fig,'fig_title',{[figtitle];[ 'Channel ' num2str(channels(i))]},...
			'min_f',min_f,'max_f',max_f,'raster_colors',mua_colors,'hist_colors',hist_colors);
	else
		imagesc(PLOTMUA.t,PLOTMUA.trials,PLOTMUA.image);
		colormap(mua_colors);
		xlabel('Time (in s)','FontSize',13,'FontName','Helvetica');
		ylabel('Trial','FontSize',13,'FontName','Helvetica');
		box off
		set(gca,'tickdir','out','linewidth',1.5,'ticklength',[.025 .025],...
			'FontSize',11,'FontName','Helvetica');
	end

	set(raster_fig,'PaperPositionMode','auto')

	multi_fig_save(raster_fig,savedir,...
		[ savefilename num2str(channels(i)) ],'eps,png');
	
	close([raster_fig]);

end

save(fullfile(savedir,'mua.mat'),'MUA','freq_range','channels');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
