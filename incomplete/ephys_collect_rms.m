function [IFR,RMS,THRESHOLD,PROC_DATA]=ephys_murate(DIR,FILELIST,varargin)
%generates song-aligned single-unit rasters
%
%	ephys_visual_sua(EPHYS.data,HISTOGRAM,EPHYS.labels,varargin)
%
%	EPHYS
%	structure with the following fields
%
%	EPHYS.data
%	sound-aligned voltage traces from extracted_data.mat (should be the variable ephys_data)
%
%	EPHYS.labels
%	channel labels (i.e. the channel that corresponds to a given element in the cell array
%	ephys_data) from extracted_data.mat
%
%	HISTOGRAM
%	contour histogram returned by ephys_visual_histogram.m (or loaded from histogram.mat)
%
%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETER COLLECTION %%%%%%%%%%%%%%%%%

% TODO:  remove multi-channel support

nparams=length(varargin);

if mod(nparams,2)>0
	error('ephysPipeline:argChk','Parameters must be specified as parameter/value pairs!');
end

noise='car'; % none, nn for nearest neighbor, or car for common average
car_exclude=[];
savedir='rms';

min_f=1; % min frequency to show for song histogram
max_f=10e3; % max frequency
hist_colors='jet'; % colormap for histogram

figtitle='';


% the processing parameter MUST match what we're using for 
% spike detection

freq_range=[400 4e3]; % bandpassing <10e3 distorted results, reasoning that >800 Hz is fine for spikes < 1ms long
filt_type='bandpass'; % high,low or bandpass
filt_order=6;
filt_name='e';
car_trim=0;

songfs=[ 2e3 8e3 ]; % vector, lower to upper fs cutoffs for song
smoothing=.2; % smoothing of song energy ratio
rthresh=1; % ratio threshold (song/nonsong)
songlen=.08; % length of extraction that contains significant song energy
winsize=.005; % fft window size (in seconds) for song detection

%proc_fs=10e3;

% remove eps generation, too slow here...

for i=1:2:nparams
	switch lower(varargin{i})
		case 'noise'
			noise=varargin{i+1};
		case 'savedir'
			savedir=varargin{i+1};
		case 'car_exclude'
			car_exclude=varargin{i+1};
		case 'filt_type'
			filt_type=varargin{i+1};
		case 'filt_name'
			filt_name=varargin{i+1};
		case 'freq_range'
			freq_range=varargin{i+1};
		case 'method'
			method=varargin{i+1};
		case 'filt_order'
			filt_order=varargin{i+1};
		case 'car_trim'
			car_trim=varargin{i+1};
	end
end

% collect all filenames

if nargin<1 | isempty(DIR)

	DIR=pwd;
end

if nargin<2
	listing=dir(fullfile(DIR,'*.mat'));
	FILELIST={listing(:).name};
end

trim_filelisting={};

% original file listing

listing=dir(fullfile(DIR,'..','..','*.mat'));
rawlist={listing(:).name};
trim_rawlist={};

for i=1:length(rawlist)
	[~,filename,~]=fileparts(rawlist{i});
	trim_rawlist{i}=filename;
end

trim_filelist={};
for i=1:length(FILELIST)
	[~,filename,~]=fileparts(FILELIST{i});
	trim_filelist{i}=filename;
end

filemap=zeros(1,length(trim_filelist));
for i=1:length(trim_rawlist)

	% find original file this was extracted from, get RMS, append to extraction file (doesn't cost
	% much in memory to do so)

	% where's the match

	tmp=strfind(trim_filelist,trim_rawlist{i});
	
	% hits
	
	hits=find(~cellfun(@isempty,tmp));
	filemap(hits)=i;

end

if any(filemap==0)
	error(['Source file not identified for all extraction files...']);
end

% get the file mapping

[uniq_files]=unique(filemap);

% get list of original files

counter=1;


for i=1:length(uniq_files)

	matches=find(filemap==uniq_files(i));

	% first load the original file

	source_file=fullfile(DIR,'..','..',rawlist{uniq_files(i)});
	load(source_file,'ephys','audio');

	% process as we would for multi-unit spike detection

	PROC_DATA=ephys_denoise_signal(ephys.data,ephys.labels,ephys.labels,'method',noise,'car_exclude',car_exclude,'car_trim',car_trim);
	PROC_DATA=ephys_condition_signal(PROC_DATA,'s','freq_range',...
		freq_range,'filt_type',filt_type,'filt_order',filt_order,'filt_name',filt_name,'fs',ephys.fs);

	[b,a]=ellip(5,.2,40,[500]/(audio.fs/2),'high');

	% first get rms for all data, then data spanning silence/song only

	micdata=filtfilt(b,a,audio.data);
	[s,f,t]=spectrogram(micdata,round(winsize*audio.fs),0,[],audio.fs);

	spect_fs=1/(t(2)-t(1));

	energy=abs(s);

	minfs=min(find(f>=songfs(1)));
	maxfs=max(find(f<=songfs(2)));

	% get the ratio of song/nonsong energy

	songenergy=mean(energy(minfs:maxfs,:)); % song energy
	nonsongenergy=mean(energy([ 1:minfs maxfs:size(energy,1)],:)); % nonsong energy
	ratio=smooth(songenergy./nonsongenergy,smoothing*spect_fs); % simple moving average

	% write out points where we detected song offline, use this to calculate thresholds
	% for how many points did we detect song?

	songpts=double(ratio>rthresh);

	timebase=[1:length(micdata)]/audio.fs;
	songdet_offline=interp1(t,songpts,timebase,'nearest');

	% interpolate from song detection to original vector

	songpts=songdet_offline>0;
	silentpts=songdet_offline==0;

	rms.channels=ephys.labels;

	rms.allpts.standard=sqrt(mean(PROC_DATA.^2));
	rms.allpts.robust=median(abs(PROC_DATA)/.6745);

	rms.songpts.standard=sqrt(mean(PROC_DATA(songpts,:).^2));
	rms.songpts.robust=median(abs(PROC_DATA(songpts,:))/.6745);
	
	rms.silentpts.standard=sqrt(mean(PROC_DATA(silentpts,:).^2));
	rms.silentpts.robust=median(abs(PROC_DATA(silentpts,:))/.6745);

	for j=1:length(matches)

		disp(['File ' num2str(counter) ' of ' num2str(length(FILELIST))]);

		extraction_file=FILELIST{matches(j)};
		save(fullfile(DIR,extraction_file),'rms','source_file','-append');

		counter=counter+1;

	end
end
