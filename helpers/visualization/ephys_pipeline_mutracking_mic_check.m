function [fig1,fig2]=ephys_pipeline_mutracking_miccheck(MICDATA,EPHYSDATA,GRPID,varargin)
%
%
%
%
%
%

% creates a scatter plot of multi-unit activity against time with a smooth average

% take the data, square and take the average

stats=[];
fs=25e3;

nparams=length(varargin);

for i=1:2:nparams
	switch lower(varargin{i})
		case 'fs'
			fs=varargin{i+1};
	end
end

% bandpass filter to recover the multi-unit data

timevec=[1:size(MICDATA,1)]/fs;

fig1=figure();
plot(timevec,sqrt(mean(MICDATA(:,GRPID==2).^2,2)),'b');
hold on;
plot(timevec,sqrt(mean(MICDATA(:,GRPID==3).^2,2)),'r');

%size(timevec)
%size(EPHYSDATA)
%
%fig2=figure();
%subplot(2,1,1);
%plot(timevec,zscore(mean(MICDATA(:,GRPID==2).^2,2)));
%hold on;
%plot(timevec,zscore(mean(EPHYSDATA(:,GRPID==2),2)));
%
%subplot(2,1,2);
%plot(timevec,zscore(mean(MICDATA(:,GRPID==3).^2,2)));
%hold on;
%plot(timevec,zscore(mean(EPHYSDATA(:,GRPID==3),2)));


