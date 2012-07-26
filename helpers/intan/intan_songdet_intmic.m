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
%

% wrap with a bash script to run whenever files appear

% detect song, place into a directory sorted by BIRD ID then DATE (two directories!), get from filename
% if we can't parse a bird id then place into bird # (iterate until non-overlapping), if song is not
% detected then junk the file, otherwise place the detected song into BIRDID/DATE/MAT

% this way we can collect all raw data into a single directory, sort into new directories and delete the rest
% should save lots of space

% add song detection parameters as well

% while running the daemon this can be changed 

mic_pre=12;
minfs=2e3;
maxfs=6e3;
ratio_thresh=4;
window=250;
noverlap=0;
song_thresh=.27; % between .2 and .3 seems to work best (higher is more exlusive)
songduration=.8;
low=5;
high=10;
colors=hot;
disp_minfs=1;
disp_maxfs=10e3;
filtering=100; % changed to 100 from 700 as a more sensible default
intan_fs=25e3;
audio_pad=.2;
config_file='';
folder_format='yyyy-mm-dd';
image_pre='gif';
wav_pre='wav';
data_pre='mat';
delimiter='\_';
nosort=0;
subdir='pretty_bird';

%logname='song_det.log';
% where to place the parsed files

root_dir=fullfile(pwd,'..','..','data','intan_data'); % where will the detected files go
proc_dir=fullfile(pwd,'..','processed'); % where do we put the files after processing
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
end


% embed in infinite loop

intlisting=dir(fullfile(DIR,'*.int'));

proc_files={};
for i=1:length(intlisting)
	proc_files{i}=fullfile(DIR,intlisting(i).name);
end

parfor i=1:length(proc_files)

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

		% fourth is date

		file_datenum=datenum([tokens{4} tokens{5}(1:end-4)],'yymmddHHMMSS');

		% now create the folder it doesn't exist already

		foldername=fullfile(root_dir,birdid,recid,datestr(file_datenum,folder_format));	

	end

	try
		[t,amps,data,aux]=read_intan_data_cli(proc_files{i});
	catch err
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

	% standard song detection

	mic_channel=find(amps==mic_trace);
	ephys_labels=setdiff(amps,mic_trace);

	ephys_channels=[];
	for j=1:length(ephys_labels)
		ephys_channels(j)=find(amps==ephys_labels(j));
	end

	if ~isempty(filtering)
		conditioned_data=filtfilt(b,a,data(:,mic_channel));
	else
		conditioned_data=data(:,mic_channel);
	end

	conditioned_data=conditioned_data./max(abs(conditioned_data));

	% did we detect song?

	try
		[song_bin]=song_det(conditioned_data,intan_fs,minfs,maxfs,window,...
			noverlap,songduration,ratio_thresh,song_thresh);
	catch err
		disp([err]);
		disp('Song detection failed, continuing...');
		fclose('all');
		continue;
	end

	[sonogram_im sonogram_f sonogram_t]=pretty_sonogram(conditioned_data,intan_fs,'n',500,'overlap',350,'low',3.5);

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

	song_pts=find(song_bin>0);

	if isempty(song_pts)
		continue;
	else
		disp(['Song detected in file:  ' proc_files{i}]);
	end

	% if we're here, we've detected song

	if ~exist(foldername,'dir')
		mkdir(foldername);
	end

	% factor to move from sonogram coordinates to raw audio data coordinates

	image_dir=fullfile(foldername,image_pre);
	wav_dir=fullfile(foldername,wav_pre);
	data_dir=fullfile(foldername,data_pre);

	if ~exist(image_dir,'dir')
		mkdir(image_dir);
	end

	if ~exist(wav_dir,'dir');
		mkdir(wav_dir);
	end

	if ~exist(data_dir,'dir');
		mkdir(data_dir);
	end

	[f,t]=size(sonogram_im);
	son_to_vec=(length(conditioned_data)-noverlap)/(length(song_bin));
	im_son_to_vec=(length(conditioned_data)-350)/t;

	% use diff to find non_continguous song bouts

	song_idx=[0 find(diff(song_pts*son_to_vec)>intan_fs) length(song_pts)];
	sonogram_filename=fullfile(image_dir,[ name '.gif' ]);

	for j=1:length(song_idx)-1

		startpoint=floor((song_pts(song_idx(j)+1))*son_to_vec-audio_pad*intan_fs);
		endpoint=ceil((song_pts(song_idx(j+1)))*son_to_vec+audio_pad*intan_fs);

		if startpoint<1, startpoint=1; end
		if endpoint>length(conditioned_data), endpoint=length(conditioned_data); end

		audio_extraction=conditioned_data(startpoint:endpoint);
		ephys_extraction=data(startpoint:endpoint,ephys_channels);

		save_name=[ name '_chunk_' num2str(j) ];

		if length(audio_extraction)<500
			continue;
		end

		sonogram_im(1:10,ceil(startpoint/im_son_to_vec):ceil(endpoint/im_son_to_vec))=63;

		[chunk_sonogram_im chunk_sonogram_f chunk_sonogram_t]=pretty_sonogram(audio_extraction,intan_fs,'low',3.5);

		startidx=max([find(chunk_sonogram_f<=disp_minfs)]);
		stopidx=min([find(chunk_sonogram_f>=disp_maxfs)]);

		chunk_sonogram_im=chunk_sonogram_im(startidx:stopidx,:);
		chunk_sonogram_im=flipdim(chunk_sonogram_im,1);

		imwrite(uint8(chunk_sonogram_im),colors,fullfile(image_dir,[ save_name '.gif']),'gif');

		min_audio=min(audio_extraction(:));
		max_audio=max(audio_extraction(:));

		if min_audio + max_audio < 0
			audio_extraction=audio_extraction./(-min_audio);
		else
			audio_extraction=audio_extraction./(max_audio*(1+1e-3));
		end

		parsave(fullfile(data_dir,['songdet1_' save_name '.mat']),ephys_extraction,audio_extraction,intan_fs,ephys_labels);
		wavwrite(audio_extraction,intan_fs,fullfile(wav_dir,[ save_name '.wav']));

	end

	reformatted_im=im_reformat(sonogram_im,10);
	imwrite(uint8(reformatted_im),colors,sonogram_filename,'gif');

end

end



function parsave(file,ephys_data,mic_data,fs,channels)

	save(file,'ephys_data','mic_data','fs','channels');

end
