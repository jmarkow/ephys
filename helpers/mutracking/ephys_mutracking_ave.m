function [fig]=ephys_pipeline_mutracking_barplot(DATA,varargin)
%
%
%
%
%
%

% creates a scatter plot of multi-unit activity against time with a smooth average

% take the data, square and take the average

stats=[];
intan_fs=25e3;
freq_range=[500 3e3];
plotcolor=[.7 .7 .7];
windowsize=15;
mads=30;
grpid=[];
grpcolors=[0 0 0;...
	0 0 1;...
	1 0 0];
startid=1; % which datenum indicates feedback start?

nparams=length(varargin);

if mod(nparams,2)>0
	error('ephysPipeline:argChk','Parameters must be specified as parameter/value pairs!');
end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'grpid'
			grpid=varargin{i+1};
		case 'startid'
			startid=varargin{i+1};
		case 'edges'
			edges=varargin{i+1};
	end
end

% bandpass filter to recover the multi-unit data


% for each group get the first and last 100-200 trials, histogram and run simple statistics (ranksum maybe)

groups=unique(grpid);

for i=1:length(groups)
	groupdata{i}=[];
end

% take the beginning and end of the session for plotting

for i=2:length(edges)

	currgrp=edges(i-1):edges(i);
	grpidx=grpid(edges(i-1));

	% take the first and last 100, or percentage

	grpidx
	len=length(edges(i-1):edges(i));

	% take a time average of the filtered LFP trace, plot 

end

