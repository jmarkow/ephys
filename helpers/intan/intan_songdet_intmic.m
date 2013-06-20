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
%		intan_fs
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
colors=hot;
disp_minfs=1;
disp_maxfs=10e3;
filtering=300; % changed to 100 from 700 as a more sensible default, leave empty to filter later
intan_fs=25e3;
audio_pad=5; % pad on either side of the extraction
error_buffer=5; % if we can't load a file, how many days old before deleting

% parameters for folder creation

folder_format='yyyy-mm-dd';
image_pre='gif';
wav_pre='wav';
data_pre='mat';
sleep_pre='sleep';

delimiter='\_';
nosort=0;
subdir='pretty_bird';
sleep_window=[ 22 7 ]; % times for keeping track of sleep data (24 hr time, start and stop)
auto_delete=1; % delete data n days old
sleep_fileinterval=3; % specify file interval (in minutes) 
sleep_segment=5; % how much data to keep (in seconds)
sleep_timestamp=[];
sleep_birdid={};
sleep_recid={};

% where to place the parsed files

root_dir=fullfile(pwd,'..','..','data','intan_data'); % where will the detected files go
proc_dir=fullfile(pwd,'..','processed'); % where do we put the files after processing, maybe auto-delete
					 % after we're confident in the operation of the pipeline
unorganized_dir=fullfile(pwd,'..','unorganized');

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
		case 'auto_delete'
			auto_delete=varargin{i+1};
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
	end
end

if nargin<1
	DIR=pwd;
end

% list the files to process

if ~isempty(filtering)
	[b,a]=butter(3,[filtering/(intan_fs/2)],'high');
else
    b=[];
    a=[];
end

intlisting=dir(fullfile(DIR,'*.int'));

proc_files={};
for i=1:length(intlisting)
	proc_files{i}=fullfile(DIR,intlisting(i).name);
end

% check all files in proc directory and delete anything older than 
% auto-delete days

if ~isempty(auto_delete)

	% get the proc_dir listing
	
	oldproc_listing=dir(fullfile(proc_dir,'*.int'));
	oldproc_dates={oldproc_listing.date};
	oldproc_filenames={oldproc_listing.name};

	% make sure recycle is set to off

	oldstate=recycle;
	newstate=recycle('off');

	for i=1:length(oldproc_dates)

		% for each date check to see if it is > than
		% the threshold

		% days elapsed

		delapsed=daysdif(datenum(oldproc_dates{i}),datenum(now));

		if delapsed>auto_delete
			disp(['Deleting ' fullfile(proc_dir,oldproc_filenames{i})]);
			delete(fullfile(proc_dir,oldproc_filenames{i}));
		end

	end

	newstate=recycle(oldstate);
    
end

for i=1:length(proc_files)

	t=[];
	amps=[];
	data=[];
	aux=[];
	song_bin=[];
	mic_trace=[];

	% read in the data

	% parse for the bird name,zone and date
	% new folder format, yyyy-mm-dd for easy sorting (on Unix systems at least)	

	if nosort
		foldername=fullfile(unorganized_dir,subdir);
		mic_trace=mic_pre;
		file_datenum=[];
	else
		tokens=regexpi(intlisting(i).name,delimiter,'split');

		% first token should be bird number

		birdid=tokens{1};

		% second should be recording tag (normally nucleus)

		recid=tokens{2};

		% third should be mic trace

		mictokens=regexpi(tokens{3},'\d+','match');

		if ~isempty(mictokens)
			mic_trace=str2num(mictokens{1});
		else
			mic_trace=mic_pre;
		end

		if length(tokens)>3
			ttl_trace=regexpi(tokens{4},'\d+','match');
		else
			ttl_trace=[];
		end

		% fourth is date

		file_datenum=datenum([tokens{4} tokens{5}(1:end-4)],'yymmddHHMMSS');

		% now create the folder it doesn't exist already

		foldername=fullfile(root_dir,birdid,recid,datestr(file_datenum,folder_format));	

		% create the bird directory

		if ~exist(fullfile(root_dir,birdid),'dir')
			mkdir(fullfile(root_dir,birdid));
		end

		% create the template directory and a little readme

		if ~exist(fullfile(root_dir,birdid,'templates'),'dir')
			mkdir(fullfile(root_dir,birdid,'templates'));
			fid=fopen(fullfile(root_dir,birdid,'templates','README.txt'),'w');
			fprintf(fid,'Templates are stored in this directory, follow these steps:\n\n');
			fprintf(fid,['1) Create a directory within this directory (templates) with' ... 
				'the name of the template (e.g. motif1)\n']);
			fprintf(fid,['2) Place template_data.mat,classify_data.mat and template.png' ... 
				'from a folder created by ephys_cluster in the subdirectory\n']);
			fprintf(fid,['3) When the pipeline is active, it will automatically' ...
				'extract examples of the template from ALL the data for a given bird']);
			fclose(fid);
		end
	end

	disp([proc_files{i}]);

	try
		[t,amps,data,aux]=read_intan_data_cli(proc_files{i});
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


	% if we're successful reading, then move the file to a processed directory

	[path,name,ext]=fileparts(proc_files{i});

	disp(['Processing ' proc_files{i}]);

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

	mic_channel=find(amps==mic_trace);
	ephys_labels=setdiff(amps,mic_trace);

	ephys_channels=[];
	for j=1:length(ephys_labels)
		ephys_channels(j)=find(amps==ephys_labels(j));
	end

	if ~isempty(ttl_trace)
		ttl_data=aux(:,ttl_trace);
	end

	conditioned_data=data(:,mic_channel);
	if ~isempty(filtering)
		norm_data=filtfilt(b,a,conditioned_data);
	else
		norm_data=conditioned_data;
	end

	norm_data=norm_data./max(abs(norm_data));

	% short circuit song processing if current file is in the sleep interval

	if ~isempty(file_datenum) & length(sleep_window)==2

		% convert the sleep window times to datenum

		[year,month,day,hour]=datevec(file_datenum);

		% compare hour, are we in the window?

		if hour>=sleep_window(1) | hour<=sleep_window(2)

			if hour<=sleep_window(2)
				new_datenum=addtodate(file_datenum,-1,'day');
			else
				new_datenum=file_datenum;
			end

			sleep_foldername=fullfile(root_dir,birdid,recid,...
				datestr(new_datenum,folder_format));

			sleep_dir=fullfile(sleep_foldername,sleep_pre);

			% time elapsed since we processed the last sleep file (in seconds)?
			% match bird and recording ID

			birdid_idxs_tmp=strfind(sleep_birdid,birdid);
			recid_idxs_tmp=strfind(sleep_recid,recid)

			birdid_idxs=~cellfun(@isempty,birdid_idxs_tmp)
			recid_idxs=~cellfun(@isempty,recid_idxs_tmp)

			match=find(birdid_idxs&recid_idxs)

			% if there is a match, check the elapsed time,
			% otherwise if we haven't processed a file with this bird and recid
			% go ahead and store

			if ~isempty(match)
				time_elapsed=etime(datevec(file_datenum),datevec(sleep_timestamp(match)))
			else
				match=length(sleep_timestamp)+1;
				sleep_birdid{match}=birdid
				sleep_recid{match}=recid
				time_elapsed=(sleep_fileinterval*60)+1;
			end         

			% is it greater than the proposed file interval?

			if time_elapsed>=sleep_fileinterval*60

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

				audio_extraction=conditioned_data(1:stopsample);
				ephys_extraction=data(1:stopsample,ephys_channels);

				parsave(fullfile(sleep_dir,['sleepdata1_' name '.mat']),...
					ephys_extraction,audio_extraction,intan_fs,ephys_labels,file_datenum);

				% if we process the file store the new timestamp

				disp(['New timestamp:  ' datestr(file_datenum)]);


				sleep_timestamp(match)=file_datenum

			end

			% finish the loop, do not move on to song detection

			disp('Continuing to song detection...');

		end
	end

	image_dir=fullfile(foldername,image_pre);
	wav_dir=fullfile(foldername,wav_pre);
	data_dir=fullfile(foldername,data_pre);

	image_dir_ttl=fullfile(foldername,[image_pre '_ttl']);
	wav_dir_ttl=fullfile(foldername,[wav_pre '_ttl']);
	data_dir_ttl=fullfile(foldername[data_pre '_ttl']);

	if ~exist(image_dir,'dir')
		mkdir(image_dir);
	end

	if ~exist(wav_dir,'dir');
		mkdir(wav_dir);
	end

	if ~exist(data_dir,'dir');
		mkdir(data_dir);
	end

	if ~exist(image_dir_ttl,'dir')
		mkdir(image_dir_ttl);
	end

	if ~exist(wav_dir_ttl,'dir');
		mkdir(wav_dir_ttl);
	end

	if ~exist(data_dir_ttl,'dir');
		mkdir(data_dir_ttl);
	end

	[sonogram_im sonogram_f sonogram_t]=pretty_sonogram(norm_data,intan_fs,'n',500,'overlap',200,'low',1.5);

	startidx=max([find(sonogram_f<=disp_minfs)]);

	if isempty(startidx)
		startidx=1;
	end

	stopidx=min([find(sonogram_f>=disp_maxfs)]);

	if isempty(stopidx)
		stopidx=length(sonogram_f);
	end

	sonogram_im=sonogram_im(startidx:stopidx,:);
	sonogram_im=flipdim(sonogram_im,1);
	[f,t]=size(sonogram_im);
	im_son_to_vec=(length(norm_data)-350)/t;

	% if we have a TTL trace, extract using the TTL

	if ~isempty(ttl_trace)

		%idx=1:length(ttl_data)-1;

		%ttl_pts=find(ttl_data(idx)<.2&ttl_data(idx+1)>.5);
		ttl_pts=find(ttl_data>.5);

		if ~isempty(ttl_pts)

			ttl_idx=[0 find(diff(ttl_pts)>audio_pad*2*intan_fs) length(ttl_pts)];

			disp(['TTL detected in file:  ' proc_files{i}]);

			for j=1:length(ttl_idx)-1

				startpoint=floor(ttl_pts(ttl_idx(j)+1)-audio_pad*intan_fs);
				endpoint=ceil(ttl_pts(ttl_idx(j+1))+audio_pad*intan_fs);

				if startpoint<1, startpoint=1; end
				if endpoint>length(norm_data), endpoint=length(norm_data); end

				sonogram_filename=fullfile(image_dir,[ name '_ttl.gif' ]);

				norm_extraction=norm_data(startpoint:endpoint);
				audio_extraction=conditioned_data(startpoint:endpoint);
				ephys_extraction=data(startpoint:endpoint,ephys_channels);
				ttl_extraction=ttl_data(startpoint:endpoint);

				save_name=[ name '_chunk_' num2str(j) '_ttl' ];

				%if length(audio_extraction)<500
				%	continue;
				%end

				sonogram_im(1:10,ceil(startpoint/im_son_to_vec):ceil(endpoint/im_son_to_vec))=63;

				[chunk_sonogram_im chunk_sonogram_f chunk_sonogram_t]=pretty_sonogram(norm_extraction,intan_fs,'n',500,'overlap',350,'low',1.5);

				startidx=max([find(chunk_sonogram_f<=disp_minfs)]);
				stopidx=min([find(chunk_sonogram_f>=disp_maxfs)]);

				chunk_sonogram_im=chunk_sonogram_im(startidx:stopidx,:);
				chunk_sonogram_im=flipdim(chunk_sonogram_im,1);

				imwrite(uint8(chunk_sonogram_im),colors,fullfile(image_dir_ttl,[ save_name '.gif']),'gif');

				parsave(fullfile(data_dir_ttl,['songdet1_' save_name '.mat']),...
					ephys_extraction,audio_extraction,ttl_extraction,intan_fs,ephys_labels,file_datenum);

				% normalize audio to write out to wav file

				min_audio=min(norm_extraction(:));
				max_audio=max(norm_extraction(:));

				if min_audio + max_audio < 0
					norm_extraction=norm_extraction./(-min_audio);
				else
					norm_extraction=norm_extraction./(max_audio*(1+1e-3));
				end

				wavwrite(norm_extraction,intan_fs,fullfile(wav_dir_ttl,[ save_name '.wav']));

			end

			reformatted_im=im_reformat(sonogram_im,(ceil((length(audio_extraction)/intan_fs)/5)));
			imwrite(uint8(reformatted_im),colors,sonogram_filename,'gif');

		end

		disp('Continuing to song detection...');


	end

	% did we detect song?

	try
		[song_bin]=song_det(norm_data,intan_fs,minfs,maxfs,window,...
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

	son_to_vec=(length(norm_data)-noverlap)/(length(song_bin));

	% use diff to find non_continguous song bouts separated by the audio pad + 1 second

	song_idx=[0 find(diff(song_pts*son_to_vec)>intan_fs+audio_pad*2*intan_fs) length(song_pts)];
	sonogram_filename=fullfile(image_dir,[ name '.gif' ]);

	for j=1:length(song_idx)-1

		startpoint=floor((song_pts(song_idx(j)+1))*son_to_vec-audio_pad*intan_fs);
		endpoint=ceil((song_pts(song_idx(j+1)))*son_to_vec+audio_pad*intan_fs);

		if startpoint<1, startpoint=1; end
		if endpoint>length(norm_data), endpoint=length(norm_data); end

		norm_extraction=norm_data(startpoint:endpoint);
		audio_extraction=conditioned_data(startpoint:endpoint);
		ephys_extraction=data(startpoint:endpoint,ephys_channels);

		if ~isempty(ttl_trace)
			ttl_extraction=ttl_data(startpoint:endpoint);
		else
			ttl_extraction=[];
		end

		save_name=[ name '_chunk_' num2str(j) ];

		if length(audio_extraction)<500
			continue;
		end

		sonogram_im(1:10,ceil(startpoint/im_son_to_vec):ceil(endpoint/im_son_to_vec))=63;

		[chunk_sonogram_im chunk_sonogram_f chunk_sonogram_t]=pretty_sonogram(norm_extraction,intan_fs,'n',500,'overlap',350,'low',3.5);

		startidx=max([find(chunk_sonogram_f<=disp_minfs)]);
		stopidx=min([find(chunk_sonogram_f>=disp_maxfs)]);

		chunk_sonogram_im=chunk_sonogram_im(startidx:stopidx,:);
		chunk_sonogram_im=flipdim(chunk_sonogram_im,1);

		imwrite(uint8(chunk_sonogram_im),colors,fullfile(image_dir,[ save_name '.gif']),'gif');

		parsave(fullfile(data_dir,['songdet1_' save_name '.mat']),...
			ephys_extraction,audio_extraction,ttl_extraction,intan_fs,ephys_labels,file_datenum);

		% normalize audio to write out to wav file

		min_audio=min(norm_extraction(:));
		max_audio=max(norm_extraction(:));

		if min_audio + max_audio < 0
			norm_extraction=norm_extraction./(-min_audio);
		else
			norm_extraction=norm_extraction./(max_audio*(1+1e-3));
		end

		wavwrite(norm_extraction,intan_fs,fullfile(wav_dir,[ save_name '.wav']));

	end

	reformatted_im=im_reformat(sonogram_im,(ceil((length(audio_extraction)/intan_fs)/5)));
	imwrite(uint8(reformatted_im),colors,sonogram_filename,'gif');

end
end

function parsave(file,ephys_data,mic_data,ttl_data,fs,channels,start_datenum)

% still saving in this function in case we go back to parfor 

if nargin<6
    start_datenum=[];
end

save(file,'ephys_data','mic_data','ttl_data','fs','channels','start_datenum');

end
