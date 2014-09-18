function [fr_smoothed]=ephys_ratetracker_plot(DATA,DATENUM,varargin)
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

% TODO:  remove multi-channel support

nparams=length(varargin);

if mod(nparams,2)>0
	error('ephysPipeline:argChk','Parameters must be specified as parameter/value pairs!');
end

channelboundary=[];
noise='none'; % none, nn for nearest neighbor, or car for common average
car_exclude=[];
savedir=pwd;

min_f=1; % min frequency to show for song histogram
max_f=10e3; % max frequency
hist_colors='jet'; % colormap for histogram

figtitle='';

lfp_bands=[1 15;20 40;40 60;60 100;100 200;200 400;500 1e3];

sigma_t=4; % multiple of noise estimate for spike threshold (generally 3-4, using Quiroga's method)
subtrials=[];
channels=[];
fr_sigma=.005;
bound=.5;
% remove eps generation, too slow here...
fbstart='';
fr_smoothing=50;
lfp_smoothing=30;
fs=500;
fig_num=[];

for i=1:2:nparams
	switch lower(varargin{i})
		case 'fr_smoothing'
			fr_smoothing=varargin{i+1};
		case 'lfp_smoothing'
			lfp_smoothing=varargin{i+1};
		case 'bound'
			bound=varargin{i+1};
		case 'fbstart'
			fbstart=varargin{i+1};
		case 'fs'
			fs=varargin{i+1};
		case 'fig_num'
			fig_num=varargin{i+1};

	end
end

if isempty(fig_num)
	fig_num=figure();
end

%%%% firing rate plot
%

time_elapsed=zeros(length(DATENUM),1);
days_since=zeros(size(time_elapsed));

if ~isempty(fbstart)	
	fbstart_datevec=datevec(fbstart,'mm-dd-yyyy HH:MM:SS');
else
	fbstart_datevec=datevec(DATENUM(1));
end

fbstart_day=fbstart_datevec(3);
% time since feedback started (in second)

for i=1:length(time_elapsed)
	tmp=datevec(DATENUM(i));
	time_elapsed(i)=etime(tmp,fbstart_datevec);
	days_since(i)=tmp(3)-fbstart_day;
end

% for each day, plot points and running average

time_elapsed=time_elapsed/(3600*24);

%bound=round(fs*bound)
%DATA=mean(TIMESERIES(bound:end-bound,:));

% interpolate to evenly spaced data, then smooth
%
fr_smoothed.y.raw=smooth(time_elapsed,DATA,fr_smoothing,'lowess');
fr_smoothed.x=time_elapsed;
fr_smoothed.timescale=fr_smoothing;
fr_smoothed.y.norm=(fr_smoothed.y.raw-min(fr_smoothed.y.raw))./(max(fr_smoothed.y.raw)-min(fr_smoothed.y.raw));

%%%% plot pre first
%
bins=unique(days_since);

if any(time_elapsed<0)
	plot(time_elapsed(time_elapsed<0),DATA(time_elapsed<0),'ko');
end

for i=1:length(bins)
	hold on;
	plot(time_elapsed(time_elapsed>=0&days_since==bins(i)),DATA(time_elapsed>=0&days_since==bins(i)),...
		'ro','markerfacecolor','r','markersize',4);
	plot(time_elapsed(time_elapsed>=0&days_since==bins(i)),fr_smoothed.y.raw(time_elapsed>=0&days_since==bins(i)),'k-','linewidth',2);
end

set(gca,'TickDir','out','FontSize',9,'FontName','Helvetica');
xlim([floor(min(time_elapsed)*1e2)/1e2 ceil(max(time_elapsed)*1e2)/1e2]);
