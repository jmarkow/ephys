function [TODEL]=ephys_pipeline_mutracking_dataclean(DATA,varargin)
%
%
%
%
%
%

% creates a scatter plot of multi-unit activity against time with a smooth average

% take the data, square and take the average

windowsize=15;
mads=20;
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
	end
end


todel=[];
for i=1:length(DATA)-windowsize

	winmed=median(DATA(i:i+windowsize-1));
	winmad=mad(DATA(i:i+windowsize-1));
	todel=[todel find((DATA>winmed+mads*winmad)|(DATA<winmed-mads*winmad))];
end

