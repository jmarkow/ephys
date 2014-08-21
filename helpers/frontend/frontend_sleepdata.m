function frontend_sleepdata(DATA,FILENAME,SLEEP_WINDOW,SLEEP_SEGMENT,SLEEP_FILEINTERVAL,SLEEP_FOLDER,BASE_DIR,FOLDER_FORMAT,DELIM,PARSE_STRING)
%frontend_sleepdata determines if data should be added to sleep collection
%
%
%
%

newstring=repmat('?',[1 length(PARSE_STRING)+1]);
datesym=strfind(PARSE_STRING,'d')+1;
newstring(datesym)='d';

fs=DATA.audio.fs;
ephys_labels=DATA.ephys.labels;
audio_extraction=[];
file_datenum=DATA.file_datenum;

% convert the sleep window times to datenum

[year,month,day,hour]=datevec(file_datenum);

% compare hour, are we in the window?

if hour<=SLEEP_WINDOW(2)
	new_datenum=addtodate(file_datenum,-1,'day');
else
	new_datenum=file_datenum;
end

sleep_foldername=fullfile(BASE_DIR,datestr(new_datenum,FOLDER_FORMAT));
sleep_dir=fullfile(sleep_foldername,SLEEP_FOLDER);
sleep_listing=dir(fullfile(sleep_dir,'*.mat'));
parameters=DATA.parameters;

% get the date number of the last saved file

if ~isempty(sleep_listing)
    tokens=regexp(sleep_listing(end).name,DELIM,'split');
    datetoken=length(tokens)-1:length(tokens);
    last_datenum=datenum([tokens{datetoken}],'yymmddHHMMSS');
	disp(['Last extraction:  ' datestr(last_datenum)]);
	time_elapsed=etime(datevec(DATA.file_datenum),datevec(last_datenum))
else
	time_elapsed=(SLEEP_FILEINTERVAL*60)+1;
end

data_types={'ephys','ttl','digout','digin','adc','aux','playback'};

savefun=@(filename,datastruct) save(filename,'-struct','datastruct','-v7.3');

% is it greater than the proposed file interval?

if time_elapsed>=SLEEP_FILEINTERVAL*60

	if ~exist(sleep_dir,'dir')
		mkdir(sleep_dir);
	end

	disp('Processing sleep file...');

	% how much data to keep?

	stopsample=round(SLEEP_SEGMENT*fs);

	disp(['Keeping ' num2str(SLEEP_SEGMENT) ' seconds of data']);

	if length(DATA.audio.data)<stopsample
		disp('File too short to keep, skipping...');
		return;
	end
	
	for j=1:length(data_types)
		if ~isempty(DATA.(data_types{j}).data)
			DATA.(data_types{j}).data=DATA.(data_types{j}).data(1:stopsample,:);
			DATA.(data_types{j}).t=DATA.(data_types{j}).t(1:stopsample);
		end
	end

	DATA.audio=rmfield(DATA.audio,'norm_data');

	DATA.audio.data=DATA.audio.data(1:stopsample);

	savefun(fullfile(sleep_dir,['sleepdata1_' FILENAME '.mat']),DATA);

end


