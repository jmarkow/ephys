function [fr lfp,playback]=ephys_ratetracker(DIR,varargin)
%generates song-aligned single-unit rasters
%
%	ephys_visual_sua(EPHYS.data,HISTOGRAM,EPHYS.labels,varargin)
%
%	EPHYS
%	structure with the following fields
%
%	EPHYS.data
%	sound-aligned voltage traces from extracted_data.mat (should be the variable ephys_data)
%
%	EPHYS.labels
%	channel labels (i.e. the channel that corresponds to a given element in the cell array
%	ephys_data) from extracted_data.mat
%
%	HISTOGRAM
%	contour histogram returned by ephys_visual_histogram.m (or loaded from histogram.mat)
%
%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETER COLLECTION %%%%%%%%%%%%%%%%%

nparams=length(varargin);

if mod(nparams,2)>0
	error('ephysPipeline:argChk','Parameters must be specified as parameter/value pairs!');
end

noise='none'; % none, nn for nearest neighbor, or car for common average
car_exclude=[];

lfp_bands=[1 15;20 40;40 60;60 100;100 200;200 400;500 1e3];

sigma_t=3; % multiple of noise estimate for spike threshold (generally 3-4, using Quiroga's method)
subtrials=[];
channels=[];
proc_fs=500;
fr_sigma=.005;

% remove eps generation, too slow here...

for i=1:2:nparams
	switch lower(varargin{i})
		case 'noise'
			noise=varargin{i+1};
		case 'savedir'
			savedir=varargin{i+1};
		case 'car_exclude'
			car_exclude=varargin{i+1};
		case 'figtitle'
			figtitle=varargin{i+1};
		case 'filt_type'
			filt_type=varargin{i+1};
		case 'filt_name'
			filt_name=varargin{i+1};
		case 'freq_range'
			freq_range=varargin{i+1};
		case 'channels'
			channels=varargin{i+1};
		case 'sigma_t'
			fr_sigma_t=varargin{i+1};
		case 'subtrials'
			subtrials=varargin{i+1};
		case 'trial_timestamps'
			trial_timestamps=varargin{i+1};
		case 'method'
			method=varargin{i+1};
		case 'filt_order'
			filt_order=varargin{i+1};
		case 'car_trim'
			car_trim=varargin{i+1};
		case 'lfp_bands'
			lfp_bands=varargin{i+1};

	end
end


[status,result]=unix(['find ' pwd ' -type f -name "aggregated_data.mat"']);
ephys_file=regexp(result,'\n','split');
ephys_file(end)=[];
load(ephys_file{1},'agg_ephys','agg_file_datenum');

if isempty(channels)
	channels=agg_ephys.labels;
end

for i=2:length(ephys_file)
	tmp=load(ephys_file{i},'agg_ephys','agg_file_datenum');
	agg_ephys.data=[ agg_ephys.data tmp.agg_ephys.data ];
	agg_file_datenum=[ agg_file_datenum tmp.agg_file_datenum ];
end

% loop over channels
%


for j=1:length(channels)

	ifr=ephys_murate(agg_ephys,'channels',channels(j),'sigma_t',sigma_t,'noise','car');

	kernedges=[-3*fr_sigma:1/agg_ephys.fs:3*fr_sigma];
	kernel=normpdf(kernedges,0,fr_sigma);
	kernel=kernel./sum(kernel);

	downfact=agg_ephys.fs/proc_fs;

	smooth_ifr=filter(kernel,1,ifr');

	fr(j).time_series=[ downsample(smooth_ifr,downfact) ];
	%fr.channels=[ agg_ephys.labels(j).*ones(1,size(smooth_ifr,2)) ];
	fr(j).datenum=[ agg_file_datenum ];
	fr(j).channels=channels(j);
	fr(j).fs=proc_fs;

end

lfp=[];

%nchannels=length(agg_ephys.labels);

%for j=find(agg_ephys.labels==20)
%
%	agg_ephys.labels(j)
%	ifr=ephys_murate(agg_ephys,'channels',agg_ephys.labels(j),'sigma_t',sigma_t,'noise','car');
%	
%	kernedges=[-3*fr_sigma:1/agg_ephys.fs:3*fr_sigma];
%	kernel=normpdf(kernedges,0,fr_sigma);
%	kernel=kernel./sum(kernel);

%	downfact=agg_ephys.fs/proc_fs;

%	smooth_ifr=filter(kernel,1,ifr');

%	fr.time_series=[ fr.time_series downsample(smooth_ifr,downfact) ];
%	fr.channels=[ fr.channels agg_ephys.labels(j).*ones(1,size(smooth_ifr,2)) ];
%	fr.datenum=[ fr.datenum agg_file_datenum ];

%	[lfptmp,lfpt]=ephys_bandpower(agg_ephys,'channels',agg_ephys.labels(j),'freq_range',lfp_bands);

%	% lfp rates etc. etc.

%	lfp.compval=[ lfp.compval lfptmp ];	
%	lfp.t = [ lfp.t lfpt ];

%end

