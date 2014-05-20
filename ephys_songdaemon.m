function ephys_songdaemon(DIR,varargin)
%runs Intan song detection indefinitely
%
%	ephys_songdaemon(DIR);
%
%	DIR
%	directory that contains files for processing
%
%
%Before running, make sure you have created the appropriate directory structure
%using ephys_pipeline_mkdirs.m.  Also, be sure to use the following file naming
%convention:
%
%BIRDID_RECORDINGID_mic#
%
%e.g. lpur35_hvc_mic12
%
%the Intan demo software automatically appends a timestamp.  Run the daemon
%in the "unprocessed" directory, where the Intan files should be saved.  As they
%are processed, they will be moved to the "processed" directory (NOTE:  FILES
%ARE NOT AUTOMATICALLY DELETED UNLESS SPECIFIED, THIS IS CONSIDERED A FEATURE 
%AND NOT A BUG AT THE MOMENT).  Song detections will be sorted by date, bird,
%and recording ID and placed in the "intan_data" folder. For detailed processing 
%options be sure to check intan_songdet_intmic.m
%
%see also ephys_pipeline_mkdir.m,intan_songdet_intmic.m
%

%
%
%

% simply loops intan_songdet_intmic in the current directory

if nargin<1
	DIR=pwd;
end

while 1==1
	intan_songdet_intmic(DIR,varargin{:});
	disp('Pausing for 10 seconds');
	pause(10);
end
