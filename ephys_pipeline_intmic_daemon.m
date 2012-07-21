function intan_songdet_intmic_daemon(DIR,varargin)
%runs Intan song detection indefinitely
%
%	intan_songdet_intmic_daemon(DIR);
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
%the Intan demo software automatically appends a timestamp
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
	intan_songdet_intmic(DIR);
	disp('Pausing for 10 seconds');
	pause(10);
end
