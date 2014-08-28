function [fr lfp]=ephys_ratetracker_plot(TIMESERIES,DATENUM,varargin)
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
proc_fs=500;
fr_sigma=.005;
bound=.25;
% remove eps generation, too slow here...
fbstart='';
fr_smoothing=50;
lfp_smoothing=30;
fs=10;


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

	end
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
% time since feedback started (in second

for i=1:length(time_elapsed)
	tmp=datevec(DATENUM(i));
	time_elapsed(i)=etime(tmp,fbstart_datevec);
	days_since(i)=tmp(3)-fbstart_day;
end

% for each day, plot points and running average

time_elapsed=time_elapsed/(3600*24);
bound=round(fs*bound);
fr_mu=mean(TIMESERIES(bound:end-bound,:));

% interpolate to evenly spaced data, then smooth
%
fr_smoothed=smooth(time_elapsed,fr_mu,fr_smoothing,'sgolay',2);

%%%% plot pre first
%
bins=unique(days_since);

if any(time_elapsed<0)
	plot(time_elapsed(time_elapsed<0),fr_mu(time_elapsed<0),'ko');
end

for i=1:length(bins)
	hold on;
	plot(time_elapsed(time_elapsed>=0&days_since==bins(i)),fr_mu(time_elapsed>=0&days_since==bins(i)),...
		'ro','markerfacecolor','r','markersize',4);
	plot(time_elapsed(time_elapsed>=0&days_since==bins(i)),fr_smoothed(time_elapsed>=0&days_since==bins(i)),'k-','linewidth',2);
end

%%%% bin time

if sum(time_elapsed<0)>0
	binx(1)=-1;
	biny{1}=fr_mu(time_elapsed<0);
else
	binx=[];
	biny={};
end


for i=1:length(bins)
	binx(end+1)=bins(i);
	biny{end+1}=fr_mu(days_since==bins(i));
end

mu=cellfun(@median,biny)
stdev=cellfun(@std,biny)
len=cellfun(@length,biny)
sem=stdev./sqrt(len)



