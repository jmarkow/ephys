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
%see also frontend_mkdirs.m,frontend_main.m
%

%
%
%

% simply loops intan_songdet_intmic in the current directory

email_monitor=0;

if nargin<1
	DIR=pwd;
end

nparams=length(varargin);

if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs!');
end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'email_monitor'
			email_monitor=varargin{i+1};
	end
end


if email_monitor
	disp(['Email monitoring enabled']);
end

mail_flag=0;

while 1==1

	tmp=dir(pwd);
	lastmod=datevec(tmp(1).datenum);

	fileopen_elapsed=etime(clock,lastmod)/60; % get time since directory last modified in minutes

	if email_monitor>0 & mail_flag==0 & fileopen_elapsed>email_monitor
		gmail_send(['An Intan file has not been created in ' num2str(fileopen_elapsed) ' minutes.']);
		mail_flag=1; % don't send another e-mail!
	end

	frontend_main(DIR,varargin{:});
	pause(10);

end
