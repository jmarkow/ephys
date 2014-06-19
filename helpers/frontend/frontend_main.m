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
ratio_thresh=4; % power ratio between song and non-song band
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
ext='int';

% parameters for folder creation

folder_format='yyyy-mm-dd';
parse_string='bimtdd'; % how to parse filenames, b=birdid, i=recid, m=micid, t=ttlid, d=date
		       % character position indicates which token (after delim split) contains the info
date_string='yymmddHHMMSS'; % parse date using datestr format

% directory names

image_pre='gif';
wav_pre='wav';
data_pre='mat';
sleep_pre='sleep';

delimiter='\_';
nosort=0;
subdir='pretty_bird';
sleep_window=[ 22 7 ]; % times for keeping track of sleep data (24 hr time, start and stop)
auto_delete_int=1; % delete data n days old
sleep_fileinterval=3; % specify file interval (in minutes) 
sleep_segment=5; % how much data to keep (in seconds)
ttl_skip=1; % skip song detection if TTL detected?

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

for i=1:length(proc_files)

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

	% read in the data

	% parse for the bird name,zone and date
	% new folder format, yyyy-mm-dd for easy sorting (on Unix systems at least)	

	if nosort
		foldername=fullfile(unorganized_dir,subdir);
		mic_trace=mic_pre;
		file_datenum=[];
	else

		% parse the file using the format string

		[birdid,recid,mic_trace,ttl_trace,file_datenum]=frontend_fileparse(proc_files{i},delimiter,parse_string,date_string);

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

	disp([repmat(hline,[2 1])]);
	disp(['Processing: ' proc_files{i}]);
	disp(['File date: ' datestr(file_datenum)]);
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

			[t,amps,data,amps_aux,aux,parameters,dig,adc]=frontend_readdata(proc_files{i});

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

	% if file contains sampling rate, overwrite and use file's fs

	if isfield(parameters,'amplifier_sample_rate')
		fs=parameters.amplifier_sample_rate;
	end

	% set up high-pass for mic data if indicated by the user
	
	if ~isempty(filtering)
	    [b,a]=butter(5,[filtering/(fs/2)],'high'); % don't need a sharp cutoff, butterworth should be fine
	else
	    b=[];
	    a=[];
	end

	% if we're successful reading, then move the file to a processed directory

	[path,name,ext]=fileparts(proc_files{i});

	try
		movefile(proc_files{i},proc_dir);
	catch
		disp(['Could not move file ' proc_files{i}]);
		fclose('all');
		continue;
	end

	if ~exist(foldername,'dir')
		mkdir(foldername);
	end

	% standard song detection

	ismic=~isempty(mic_trace);
	isttl=~isempty(ttl_trace);

	disp(['Flags: mic ' num2str(ismic) ' ttl ' num2str(isttl)]);


	if isttl
		ttl_data=aux(:,find(ttl_trace==amps_aux));
	end

	if ismic		

		mic_channel=find(amps==mic_trace);
		ephys_labels=setdiff(amps,mic_trace);
		
		conditioned_data=data(:,mic_channel);

		if ~isempty(filtering)
			norm_data=filtfilt(b,a,conditioned_data);
		else
			norm_data=conditioned_data;
		end

		norm_data=norm_data./max(abs(norm_data));
	else
		ephys_labels=amps;
	end

	ephys_channels=[];
	for j=1:length(ephys_labels)
		ephys_channels(j)=find(amps==ephys_labels(j));
	end

	% prepare data as a structure to pass to processing functions

	datastruct=struct('data',data,'norm_data',norm_data,'conditioned_data',conditioned_data,...
		'ttl_data',ttl_data,'ephys_labels',ephys_labels,'ephys_channels',ephys_channels,...
		'datenum',file_datenum,'fs',fs,'aux',aux,'amps_aux',amps_aux,'dig',dig,'adc',adc,...
		'parameters',parameters);

	clearvars data norm_data conditioned_data ttl_data;

	if ~isempty(file_datenum) & length(sleep_window)==2

		% convert the sleep window times to datenum

		[~,~,~,hour]=datevec(file_datenum);

		% compare hour, are we in the window?

		if hour>=sleep_window(1) | hour<=sleep_window(2)

			disp(['Processing sleep data for file ' proc_files{i}]);

			frontend_sleepdata(datastruct,name,sleep_window,sleep_segment,sleep_fileinterval,sleep_pre,...
				fullfile(root_dir,birdid,recid),folder_format,delimiter,parse_string);	
			
			sleep_flag=1;

			% skip song detection?

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

		ephys_extraction=datastruct.data;
		audio_extraction=[];
		ttl_extraction=[];

		save(fullfile(data_dir,['songdet1_' name '.mat']),...
			'ephys_extraction','ttl_extraction','audio_extraction','fs','ephys_labels','file_datenum',...
			'parameters','dig','adc','-v7.3');

		continue;

	end

	% if we have a TTL trace, extract using the TTL
	
	dirstructttl=struct('image',image_dir_ttl,'wav',wav_dir_ttl,'data',data_dir_ttl);
	dirstruct=struct('image',image_dir,'wav',wav_dir,'data',data_dir);

	if isttl

		ttl_pts=find(ttl_data>.5)';

		if ~isempty(ttl_pts)

			ttl_idx=[0 find(diff(ttl_pts)>audio_pad*2*fs) length(ttl_pts)];

			idx=1:length(ttl_idx)-1;

			startpoints=floor(ttl_pts(ttl_idx(idx)+1)-audio_pad*fs);
			stoppoints=ceil(ttl_pts(ttl_idx(idx+1))+audio_pad*fs);

			ext_pts=[startpoints(:) stoppoints(:)];

			disp(['TTL detected in file:  ' proc_files{i}]);

			frontend_dataextract(name,datastruct,dirstructttl,ext_pts,disp_minfs,disp_maxfs,1,colors);

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
			[song_bin]=song_det(datastruct.norm_data,fs,minfs,maxfs,window,...
				noverlap,songduration,ratio_thresh,song_thresh);
		catch err
			disp([err]);
			disp('Song detection failed, continuing...');
			fclose('all');
			continue;
		end


		song_pts=find(song_bin>0);

		if isempty(song_pts)
			continue;
		else
			disp(['Song detected in file:  ' proc_files{i}]);
		end

		% if we're here, we've detected song
		% factor to move from sonogram coordinates to raw audio data coordinates

		son_to_vec=(length(datastruct.norm_data)-noverlap)/(length(song_bin));

		% use diff to find non_continguous song bouts separated by the audio pad + 1 second

		song_idx=[0 find(diff(song_pts*son_to_vec)>fs+audio_pad*2*fs) length(song_pts)];
		
		idx=1:length(song_idx)-1;

		startpoints=floor(song_pts(song_idx(idx)+1)*son_to_vec-audio_pad*fs);
		stoppoints=ceil(song_pts(song_idx(idx+1))*son_to_vec+audio_pad*fs);

		ext_pts=[startpoints(:) stoppoints(:)];

		frontend_dataextract(name,datastruct,dirstruct,ext_pts,disp_minfs,disp_maxfs,0,colors);

	end

	% if there is neither a mic nor a TTL signal, store everything

	clearvars datastruct dirstruct dirstructttl;

end
