function frontend_sleepdata(DATA,CONDITIONED_DATA,EPHYS_LABELS,FS,FILE_DATENUM,SLEEP_SEGMENT,SLEEP_FILEINTERVAL,BASE_DIR,FOLDER_FORMAT)
%frontend_sleepdata determines if data should be added to sleep collection
%
%
%
%

fs=FS;
ephys_labels=EPHYS_LABELS;
audio_extraction=[];

% convert the sleep window times to datenum

[year,month,day,hour]=datevec(FILE_DATENUM);

% compare hour, are we in the window?

if hour>=SLEEP_WINDOW(1) | hour<=SLEEP_WINDOW(2)

	disp(['Processing sleep data for file ' FILENAME]);

	if hour<=SLEEP_WINDOW(2)
		new_datenum=addtodate(FILE_DATENUM,-1,'day');
	else
		new_datenum=FILE_DATENUM;
	end

	sleep_foldername=fullfile(BASE_DIR,datestr(new_datenum,FOLDER_FORMAT));

	sleep_foldername=fullfile(ROOT_DIR,BIRDID,RECID,...
		datestr(new_datenum,FOLDER_FORMAT));
	sleep_dir=fullfile(sleep_foldername,sleep_pre);

	sleep_listing=dir(fullfile(sleep_dir,'*.mat'));

	% get the date number of the last saved file

	if ~isempty(sleep_listing)
		[~,~,~,~,last_datenum]=frontend_fileparse(sleep_listing(end).name,'????dd');
		time_elapsed=etime(datevec(FILE_DATENUM),datevec(last_datenum));
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

		stopsample=round(sleep_segment*intan_fs);

		disp(['Keeping ' num2str(sleep_segment) ' seconds of data']);

		if length(conditioned_data)<stopsample
			disp('File too short to keep, skipping...');
			continue;
		end

		if ~isempty(CONDITION_DATA)
			audio_extraction=CONDITIONED_DATA(1:stopsample);
		end

		ephys_extraction=DATA(1:stopsample,ephys_channels);

		if ~isempty(ttl_trace)
		    ttl_extraction=ttl_data(1:stopsample);
		else
		    ttl_extraction=[];
		end
		
		save(fullfile(sleep_dir,['sleepdata1_' name '.mat']),...
			'ephys_extraction','audio_extraction','ttl_extraction','fs','ephys_labels','FILE_DATENUM');

	end

end

