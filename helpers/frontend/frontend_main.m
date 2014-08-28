function intan_songdet_intmic(DIR,varargin)
%intan_songdet_intmic.m is the core script for processing Intan files
%on the fly.  
%
%	intan_songdet_intmic(DIR,varargin)
%
%	DIR
%	directory to process 
%
%	the following may be specified as parameter/value pairs
%
%		mic_pre
%		mic channel (only need to specify if nosort==1)
%
%		minfs
%		minimum fs used for song detection (default: 2e3)
%
%		maxfs
%		maximum fs used for song detection (default: 6e3)
%
%		ratio_thresh
%		ratio between song frequencies and non-song frequencies for song detection (default: 4)
%
%		window
%		spectrogram window for song detection (default: 250 samples)
%		
%		noverlap
%		window overlap for song detection (default: 0)
%
%		song_thresh
%		song threshold (default: .27)
%	
%		songduration
%		song duration for song detection in secs (default: .8 seconds)
%
%		low
%		parameter for spectrogram display (default: 5), lower if spectrogram are dim		
%
%		high
%		parameter for spectrogram display (default: 10)
%
%		colors
%		spectrogram colormap (default: hot)		
%
%		disp_minfs
%		minimum fs for spectrograms (default: 1e3)		
%
%		disp_maxfs
%		maximum fs for spectrograms (default: 7e3)		
%
%		filtering
%		high pass corner for mic trace (default: 700 Hz)
%
%		fs
%		Intan sampling rate (default: 25e3)
%
%		audio_pad
%		extra data to left and right of extraction points to extract (default: .2 secs)
%
%		folder_format
%		folder format (date string) (default: yyyy-mm-dd)
%
%		image_pre
%		image sub directory (default: 'gif')
%	
%		wav_pre
%		wav sub directory (default: 'wav')
%
%		data_pre
%		data sub directory (default: 'mat')
%	
%		delimiter
%		delimiter for filename parsing (default: '\_', or underscore)
%
%		nosort
%		set to 1 to not parse filename (data put into separate folder) (default: 0)
%
%		subdir
%		subdir if nosort==1 (default: 'pretty bird')
%
% see also ephys_pipeline_intmic_daemon.m, song_det.m, im_reformat.m, ephys_pipeline_mkdirs.m
%
%
% To run this in daemon mode, run ephys_pipeline_intmic_daemon.m in the directory with unprocessed Intan
% files.  Be sure to create the appropriate directory structure using epys_pipeline_mkdirs.m first.

% while running the daemon this can be changed 

mic_pre=12;
minfs=2e3; % the song 'band'
maxfs=6e3; % the song 'band'
ratio_thresh=2; % power ratio between song and non-song band
window=250; % window to calculate ratio in (samples)
noverlap=0; % just do no overlap, faster
song_thresh=.15; % between .2 and .3 seems to work best (higher is more exlusive)
songduration=.8; % moving average of ratio
low=5;
high=10;
colors='hot';
disp_minfs=1;
disp_maxfs=10e3;
filtering=300; % changed to 100 from 700 as a more sensible default, leave empty to filter later
fs=25e3;
audio_pad=5; % pad on either side of the extraction
error_buffer=5; % if we can't load a file, how many days old before deleting
ext='rhd';

% parameters for folder creation

folder_format='yyyy-mm-dd';
parse_string='bimpdd'; % how to parse filenames, b=birdid, i=recid, m=micid, t=ttlid, d=date
		       % character position indicates which token (after delim split) contains the info
date_string='yymmddHHMMSS'; % parse date using datestr format

% directory names

image_pre='gif';
wav_pre='wav';
data_pre='mat';
sleep_pre='sleep';

delimiter='\_'; % delimiter for splitting fields in filename
bird_delimiter='\|'; % delimiter for splitting multiple birds
nosort=0;
subdir='pretty_bird';
sleep_window=[ 22 7 ]; % times for keeping track of sleep data (24 hr time, start and stop)
auto_delete_int=1; % delete data n days old
sleep_fileinterval=10; % specify file interval (in minutes) 
sleep_segment=5; % how much data to keep (in seconds)
ttl_skip=1; % skip song detection if TTL detected?
email_monitor=0; % monitor file creation, email if no files created in email_monitor minutes
file_check=.2; % how long to wait between file reads to check if file is no longer being written (in seconds)

mfile_path = mfilename('fullpath');
[script_path,~,~]=fileparts(mfile_path);

% where to place the parsed files

root_dir=fullfile(pwd,'..','..','data','intan_data'); % where will the detected files go
proc_dir=fullfile(pwd,'..','processed'); % where do we put the files after processing, maybe auto-delete
					 % after we're confident in the operation of the pipeline
unorganized_dir=fullfile(pwd,'..','unorganized');

hline=repmat('#',[1 80]);

if ~exist(root_dir,'dir')
	mkdir(root_dir);
end

if ~exist(proc_dir,'dir')
	mkdir(proc_dir);
end

% directory for files that have not been recognized

if ~exist(unorganized_dir,'dir');
	mkdir(unorganized_dir);
end

% we should write out a log file with filtering parameters, when we started, whether song was
% detected in certain files, etc.

nparams=length(varargin);

if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs!');
end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'auto_delete_int'
			auto_delete_int=varargin{i+1};
		case 'sleep_window'
			sleep_window=varargin{i+1};
		case 'sleep_fileinterval'
			sleep_fileinterval=varargin{i+1};
		case 'sleep_segment'
			sleep_segment=varargin{i+1};
		case 'filtering'
			filtering=varargin{i+1};
		case 'audio_pad'
			audio_pad=varargin{i+1};
		case 'song_thresh'
			song_thresh=varargin{i+1};
		case 'error_buffer'
			error_buffer=varargin{i+1};
		case 'mic_pre'
			mic_pre=varargin{i+1};
		case 'colors'
			colors=varargin{i+1};
		case 'folder_format'
			folder_format=varargin{i+1};
		case 'delimiter'
			delimiter=varargin{i+1};
		case 'nosort'
			nosort=varargin{i+1};
		case 'subdir'
			subdir=varargin{i+1};
		case 'ttl_skip'
			ttl_skip=varargin{i+1};
		case 'ext'
			ext=varargin{i+1};
		case 'parse_string'
			parse_string=varargin{i+1};
		case 'email_monitor'
			email_monitor=varargin{i+1};
	end
end

if nargin<1
	DIR=pwd;
end

% read in int or rhd files

filelisting=dir(fullfile(DIR));

% delete directories

isdir=cat(1,filelisting(:).isdir);
filelisting(isdir)=[];

% read in appropriate suffixes 

filenames={filelisting(:).name};
hits=regexp(filenames,'\.(rhd|int)','match');
hits=cellfun(@length,hits)>0;

filenames(~hits)=[];

proc_files={};
for i=1:length(filenames)
	proc_files{i}=fullfile(DIR,filenames{i});
end

clear filenames;

% check all files in proc directory and delete anything older than 
% auto-delete days

if ~isempty(auto_delete_int)
	auto_delete(proc_dir,auto_delete_int,ext);    
end

fileopen_time1=clock; % get the current time to track file creation
mail_flag=0;

for i=1:length(proc_files)


	fclose('all'); % seems to be necessary

	% read in the data

	% parse for the bird name,zone and date
	% new folder format, yyyy-mm-dd for easy sorting (on Unix systems at least)	

	disp([repmat(hline,[2 1])]);
	disp(['Processing: ' proc_files{i}]);

	% try reading the file, if we fail, skip

	%%% check if file is still being written to, check byte change within N msec

	dir1=dir(proc_files{i});
	pause(file_check);
	dir2=dir(proc_files{i});

	bytedif=dir1.bytes-dir2.bytes;

	% if we haven't written any new data in the past (file_check) seconds, assume
	% file has been written

	if bytedif==0

		try

			datastruct=frontend_readdata(proc_files{i});
			datastruct.original_filename=proc_files{i};

			fileopen_time2=clock;
			fileopen_elapsed=etime(fileopen_time2,fileopen_time1)/60; % elapsed time in minutes

			disp(['Time since last file successfully opened (mins):  ' num2str(fileopen_elapsed)]);

			if email_monitor>0 & mail_flag==0
				if fileopen_elapsed>email_monitor
					gmail_send(['An Intan file has not been created in ' num2str(fileopen_elapsed) ' minutes.']);
					mail_flag=1; % don't send another e-mail!
				end

			end

			fileopen_time1=fileopen_time2;

		catch err

			file_age=daysdif(file_datenum,datenum(now));

			if file_age>error_buffer
				disp(['File too old and cannot process, deleting ' proc_files{i}]);
				delete(proc_files{i});
				continue;
			end

			disp([err])
			disp('Could not read file, continuing...');
			fclose('all'); % read_intan does not properly close file if it bails
			continue;
		end
	else
		disp('File still being written, continuing...');
		continue;
	end


	% if we're successful reading, then move the file to a processed directory

	[path,name,ext]=fileparts(proc_files{i});

	% if user passes multiple birds, they are split by bird_delimiter, parsing is done
	% independently for each bird

	bird_split=regexp(name,bird_delimiter,'split');

	tokens=regexp(bird_split{end},delimiter,'split');

	% get the date tokens from the last bird, append to all others

	%datetokens=find(parse_string=='d');
	datetokens=[length(tokens)-1 length(tokens)];
	datestring='';

	for j=1:length(datetokens)
		datestring=[ datestring delimiter(end) tokens{datetokens(j)} ];
	end

	nbirds=length(bird_split);

	% clear out all extraction variables to be safe

	found_ports=unique(datastruct.ephys.ports); % which ports are currently being used?

	disp(['Found ports:  ' found_ports]);

	for j=1:nbirds

		sleep_flag=0;
		song_bin=[];

		norm_data=[];
		conditioned_data=[];
		ttl_data=[];

		norm_extraction=[];
		audio_extraction=[];
		ephys_extraction=[];
		ttl_extraction=[];
		sonogram_im=[];
		chunk_sonogram_im=[];

		if nosort
			foldername=fullfile(unorganized_dir,subdir);
			mic_trace=mic_pre;
			file_datenum=[];
		else

			% parse the file using the format string

			if j<nbirds
				bird_split{j}=[bird_split{j} datestring];
			end

			[birdid,recid,mic_trace,mic_source,mic_port,ports,ttl_trace,ttl_source,...
				playback_trace,playback_source,file_datenum]=...
				frontend_fileparse(bird_split{j},delimiter,parse_string,date_string);

			disp(['Processing bird ' num2str(j) ' of ' num2str(nbirds) ]);
			disp(['File date: ' datestr(file_datenum)]);
			disp(['Bird ID:  ' birdid]);
			disp(['Rec ID:  ' recid]);
			disp(['Mic ch:  ' num2str(mic_trace)]);
			disp(['Mic source:  ' mic_source]);
			disp(['Mic port:  ' mic_port]);
			disp(['Playback ch:  ' num2str(playback_trace)]);
			disp(['Playback source:  ' playback_source]);
			disp(['Data ports:  ' ports]);

			% now create the folder it doesn't exist already

			foldername=fullfile(root_dir,birdid,recid,datestr(file_datenum,folder_format));	

			% create the bird directory

			if ~exist(fullfile(root_dir,birdid),'dir')
				mkdir(fullfile(root_dir,birdid));
			end

			% create the template directory and a little readme

			if ~exist(fullfile(root_dir,birdid,'templates'),'dir')
				mkdir(fullfile(root_dir,birdid,'templates'));
				copyfile(fullfile(script_path,'template_readme.txt'),...
					fullfile(root_dir,birdid,'templates','README.txt'));
			end
		end

		if ~isempty(ports)
		
			include_ports=[];

			for k=1:length(ports)

				if any(ismember(lower(found_ports(:)),lower(ports(k))))
					include_ports=[include_ports ports(k)];
				end
			end
		else
			include_ports=found_ports;
		end

		include_ports=upper(include_ports);

		disp(['Will extract from ports: ' include_ports]);

		include_ephys=[];
		include_aux=[];
		include_id='';

		for k=1:length(include_ports)
			
			len=length(find(datastruct.ephys.ports==include_ports(k)));

			include_ephys=[include_ephys find(datastruct.ephys.ports==include_ports(k))];
			include_aux=[include_aux find(datastruct.aux.ports==include_ports(k))];
			include_id=[include_id repmat(include_ports(k),[1 len])];

		end

		fprintf(1,'Raw channel mapping for port: ');
		
		for k=1:length(include_ephys)
			fprintf(1,'%i(%s) ',include_ephys(k),include_id(k));	
		end

		fprintf(1,'\n');

		% map to a new structure with the appropriate ports

		datastruct.file_datenum=file_datenum;

		birdstruct=datastruct;

		birdstruct.ephys.labels=birdstruct.ephys.labels(include_ephys);
		birdstruct.ephys.ports=birdstruct.ephys.ports(include_ephys);
		birdstruct.ephys.data=birdstruct.ephys.data(:,include_ephys);

		birdstruct.aux.labels=birdstruct.aux.labels(include_aux);
		birdstruct.aux.ports=birdstruct.aux.ports(include_aux);
		birdstruct.aux.data=birdstruct.aux.data(:,include_aux);

		% if file contains sampling rate, overwrite and use file's fs

		if ~exist(foldername,'dir')
			mkdir(foldername);
		end

		% standard song detection

		ismic=~isempty(mic_trace);
		isttl=~isempty(ttl_trace);
		isplayback=~isempty(playback_trace);

		disp(['Flags: mic ' num2str(ismic) ' ttl ' num2str(isttl) ' playback ' num2str(isplayback)]);

		% if we use a ttl trigger, assume the source is digital

		if isplayback
			switch(lower(playback_source(1)))

				case 'c'

					playback_channel=find(playback_trace==birdstruct.adc.labels);

					birdstruct.playback.data=birdstruct.adc.data(:,playback_channel);
					birdstruct.playback.fs=birdstruct.adc.fs;
					birdstruct.playback.t=birdstruct.adc.t;

					birdstruct.adc.data(:,playback_channel)=[];
					birdstruct.adc.labels(playback_channel)=[];

					if isempty(birdstruct.adc.data)
						birdstruct.adc.t=[];
					end


				case 'd'

					playback_channel=find(playback_trace==birdstruct.digin.labels);

					birdstruct.playback.data=birdstruct.digin.data(:,playback_channel);
					birdstruct.playback.fs=birdstruct.digin.fs;
					birdstruct.playback.t=birdstruct.digin.t;

					birdstruct.digin.data(:,playback_channel)=[];
					birdstruct.digin.labels(playback_channel)=[];

					if isempty(birdstruct.digin.data)
						birdstruct.digin.t=[];
					end


				end
			else
				birdstruct.playback.data=[];
			end

			if isttl

				switch lower(ttl_source(1))

				case 'c'

					ttl_channel=find(ttl_trace==birdstruct.adc.labels);

					birdstruct.ttl.data=birdstruct.adc.data(:,ttl_channel);
					birdstruct.ttl.fs=birdstruct.adc.fs;
                    			birdstruct.ttl.t=birdstruct.adc.t;
	
			   	    	birdstruct.adc.data(:,ttl_channel)=[];
					birdstruct.adc.labels(ttl_channel)=[];

					if isempty(birdstruct.adc.data)
						birdstruct.adc.t=[];
					end

                    
				case 'd'

					ttl_channel=find(ttl_trace==birdstruct.digin.labels);

					birdstruct.ttl.data=birdstruct.digin.data(:,ttl_channel);
					birdstruct.ttl.fs=birdstruct.digin.fs;
                    			birdstruct.ttl.t=birdstruct.digin.t;

					birdstruct.digin.data(:,ttl_channel)=[];
					birdstruct.digin.labels(ttl_channel)=[];
					
					if isempty(birdstruct.digin.data)
						birdstruct.digin.t=[];
					end


			end
		else
				birdstruct.ttl.data=[];
		end

		if ismic		

			% (m)ain channels (i.e. electrode channel), (a)ux or a(d)c?

			switch lower(mic_source(1))

				case 'm'

					mic_channel=find(birdstruct.ephys.labels==mic_trace&birdstruct.ephys.ports==mic_port);

					% take out the mic channel from the ephys labels

					birdstruct.audio.data=birdstruct.ephys.data(:,mic_channel);
					birdstruct.audio.fs=birdstruct.ephys.fs;
					birdstruct.audio.t=birdstruct.ephys.t;

					birdstruct.ephys.data(:,mic_channels)=[];
					birdstruct.ephys.labels(mic_channel)=[];

					if isempty(birdstruct.ephys.data)
						birdstruct.ephys.t=[];
					end

				case 'a'

					mic_channel=find(birdstruct.aux.labels==mic_trace&birdstruct.aux.ports==mic_port);

					birdstruct.audio.data=birdstruct.aux.data(:,mic_channel);
					birdstruct.audio.fs=birdstruct.aux.fs;
					birdstruct.audio.t=birdstruct.aux.t;

					birdstruct.aux.data(:,mic_channel)=[];
					birdstruct.aux.labels(mic_channel)=[];

					if isempty(birdstruct.aux.data)
						birdstruct.aux.t=[];
					end

				case 'c'

					mic_channel=find(birdstruct.adc.labels==mic_trace);

					birdstruct.audio.data=birdstruct.adc.data(:,mic_channel);
					birdstruct.audio.fs=birdstruct.adc.fs;
					birdstruct.audio.t=birdstruct.adc.t;

					birdstruct.adc.data(:,mic_channel)=[];
					birdstruct.adc.labels(mic_channel)=[];

					if isempty(birdstruct.adc.data)
						birdstruct.adc.t=[];
					end

			end

			% set up high-pass for mic data if indicated by the user

			if ~isempty(filtering)
				[b,a]=butter(5,[filtering/(fs/2)],'high'); % don't need a sharp cutoff, butterworth should be fine
			else
				b=[];
				a=[];
			end


			if ~isempty(filtering)
				birdstruct.audio.norm_data=filtfilt(b,a,birdstruct.audio.data);
			else
				birdstruct.audio.norm_data=birdstruct.audio.data;
			end

			birdstruct.audio.norm_data=birdstruct.audio.norm_data./max(abs(birdstruct.audio.norm_data));
		else
			birdstruct.audio.data=[];
			birdstruct.audio.norm_data=[];

		end


		if ~isempty(file_datenum) & length(sleep_window)==2

			% convert the sleep window times to datenum

			[~,~,~,hour]=datevec(file_datenum);

			% compare hour, are we in the window?

			if hour>=sleep_window(1) | hour<=sleep_window(2)

				disp(['Processing sleep data for file ' proc_files{i}]);

				frontend_sleepdata(birdstruct,bird_split{j},sleep_window,sleep_segment,sleep_fileinterval,sleep_pre,...
					fullfile(root_dir,birdid,recid),folder_format,delimiter,parse_string);	

				sleep_flag=1;

				% TODO: skip song detection?

			end
		end

		frontend_extract_mkdirs(foldername,image_pre,wav_pre,data_pre,isttl);

		image_dir=fullfile(foldername,image_pre);
		wav_dir=fullfile(foldername,wav_pre);
		data_dir=fullfile(foldername,data_pre);

		image_dir_ttl=fullfile(foldername,[image_pre '_ttl']);
		wav_dir_ttl=fullfile(foldername,[wav_pre '_ttl']);
		data_dir_ttl=fullfile(foldername,[data_pre '_ttl']);

		if ~ismic & ~isttl & ~sleep_flag

			save(fullfile(data_dir,['songdet1_' bird_split{j} '.mat']),'-struct','birdstruct','-v7.3');
			clearvars birdstruct;

			continue;

		end

		% if we have a TTL trace, extract using the TTL

		dirstructttl=struct('image',image_dir_ttl,'wav',wav_dir_ttl,'data',data_dir_ttl);
		dirstruct=struct('image',image_dir,'wav',wav_dir,'data',data_dir);

		if isttl

			ttl_pts=find(birdstruct.ttl.data(:)>.5)';

			if ~isempty(ttl_pts)

				ttl_idx=[0 find(diff(ttl_pts)>audio_pad*2*fs) length(ttl_pts)];

				idx=1:length(ttl_idx)-1;

				startpoints=floor(ttl_pts(ttl_idx(idx)+1)-audio_pad*fs);
				stoppoints=ceil(ttl_pts(ttl_idx(idx+1))+audio_pad*fs);

				ext_pts=[startpoints(:) stoppoints(:)];

				if ~isempty(ext_pts)
					disp(['TTL detected in file:  ' proc_files{i}]);
				end

				frontend_dataextract(bird_split{j},birdstruct,dirstructttl,ext_pts,disp_minfs,disp_maxfs,1,colors);

				% if we found TTL pulses and ttl_skip is on, skip song detection and move on to next file

				if ttl_skip
					disp('Skipping song detection...');
					continue;
				end	

			end
		end

		% did we detect song?

		if ismic

			try
				disp('Entering song detection...');
				[song_bin]=song_det(birdstruct.audio.norm_data,fs,minfs,maxfs,window,...
					noverlap,songduration,ratio_thresh,song_thresh);
			catch err
				disp([err]);
				disp('Song detection failed, continuing...');
				fclose('all');
				continue;
			end


			song_pts=find(song_bin>0);

			if isempty(song_pts)
				disp(['No song detected in file:  ' proc_files{i}]);
				continue;
			else
				disp(['Song detected in file:  ' proc_files{i}]);
			end

			% if we're here, we've detected song
			% factor to move from sonogram coordinates to raw audio data coordinates

			son_to_vec=(length(birdstruct.audio.norm_data)-noverlap)/(length(song_bin));

			% use diff to find non_continguous song bouts separated by the audio pad + 1 second

			song_idx=[0 find(diff(song_pts*son_to_vec)>fs+audio_pad*2*fs) length(song_pts)];

			idx=1:length(song_idx)-1;

			startpoints=floor(song_pts(song_idx(idx)+1)*son_to_vec-audio_pad*fs);
			stoppoints=ceil(song_pts(song_idx(idx+1))*son_to_vec+audio_pad*fs);

			ext_pts=[startpoints(:) stoppoints(:)];

			frontend_dataextract(bird_split{j},birdstruct,dirstruct,ext_pts,disp_minfs,disp_maxfs,0,colors);

		end

		% clear the datastructure for this bird

		clear birdstruct;

	end

	% if there is neither a mic nor a TTL signal, store everything

	clearvars datastruct dirstruct dirstructttl;

	try
		movefile(proc_files{i},proc_dir);
	catch
		disp(['Could not move file ' proc_files{i}]);
		fclose('all');
		continue;
	end

end
