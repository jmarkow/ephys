function frontend_sleepdata(DATA,FILENAME,SLEEP_WINDOW,SLEEP_SEGMENT,SLEEP_FILEINTERVAL,SLEEP_FOLDER,BASE_DIR,FOLDER_FORMAT,DELIM,PARSE_STRING)
%frontend_sleepdata determines if data should be added to sleep collection
%
%
%
%

newstring=repmat('?',[1 length(PARSE_STRING)+1]);
datesym=strfind(PARSE_STRING,'d')+1;
newstring(datesym)='d';

fs=DATA.fs;
ephys_labels=DATA.ephys_labels;
audio_extraction=[];
file_datenum=DATA.datenum;

% convert the sleep window times to datenum

[year,month,day,hour]=datevec(DATA.datenum);

% compare hour, are we in the window?

if hour<=SLEEP_WINDOW(2)
	new_datenum=addtodate(DATA.datenum,-1,'day');
else
	new_datenum=DATA.datenum;
end

sleep_foldername=fullfile(BASE_DIR,datestr(new_datenum,FOLDER_FORMAT));
sleep_dir=fullfile(sleep_foldername,SLEEP_FOLDER);
sleep_listing=dir(fullfile(sleep_dir,'*.mat'));
parameters=DATA.parameters;

% get the date number of the last saved file

if ~isempty(sleep_listing)
	[~,~,~,~,last_datenum]=frontend_fileparse(sleep_listing(end).name,DELIM,newstring);
	time_elapsed=etime(datevec(DATA.datenum),datevec(last_datenum));
else
	time_elapsed=(SLEEP_FILEINTERVAL*60)+1;
end

% is it greater than the proposed file interval?

if time_elapsed>=SLEEP_FILEINTERVAL*60

	if ~exist(sleep_dir,'dir')
		mkdir(sleep_dir);
	end

	disp('Processing sleep file...');

	% how much data to keep?

	stopsample=round(SLEEP_SEGMENT*fs);

	disp(['Keeping ' num2str(SLEEP_SEGMENT) ' seconds of data']);

	if length(DATA.conditioned_data)<stopsample
		disp('File too short to keep, skipping...');
		return;
	end

	if ~isempty(DATA.conditioned_data)
		audio_extraction=DATA.conditioned_data(1:stopsample);
	end

	ephys_extraction=DATA.data(1:stopsample,DATA.ephys_channels);

	if ~isempty(DATA.ttl_data)
	    ttl_extraction=DATA.ttl_data(1:stopsample);
	else
	    ttl_extraction=[];
	end
	
	save(fullfile(sleep_dir,['sleepdata1_' FILENAME '.mat']),...
		'ephys_extraction','audio_extraction','ttl_extraction','fs','ephys_labels','file_datenum',...
		'parameters','-v7.3');

end


