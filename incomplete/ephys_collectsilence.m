function [EPHYS_DATA,MIC_DATA,START_DATENUM,CHANNELS]=ephys_collectsilence(DIR,varargin)
%
%
%
%

% DIR specifies directory to process, GAP specifies how long 

if nargin<1
	error('ephysPipeline:tfhistogram:notenoughparams','Need 1 argument to continue, see documentation');
end

nparams=length(varargin);

if mod(nparams,2)>0
	error('ephysPipeline:argChk','Parameters must be specified as parameter/value pairs!');
end

fs=25e3;
savedir='silence'; % where to save the new data
min_dist=4; % minimum distance from any vocalization (in seconds,typically 4 seconds)
seg_length=1; % how long should the segments be (typically 1 second)
max_trials=300; % maximum number of trials to collect

minfs=2e3; % the song 'band'
maxfs=6e3; % the song 'band'
ratio_thresh=4; % power ratio between song and non-song band
window=250; % window to calculate ratio in (samples)
noverlap=0; % just do no overlap, faster
song_thresh=.1; % between .2 and .3 seems to work best (higher is more exlusive)
songduration=.7; % moving average of ratio

micfilter=300; % mic high-pass

for i=1:2:nparams
	switch lower(varargin{i})
		case 'min_dist'
			min_dist=varargin{i+1};
		case 'seg_length'
			seg_length=varargin{i+1};
		case 'savedir'
			savedir=varargin{i+1};
		case 'max_trials'
			max_trials=varargin{i+1};
		case 'song_thresh'
			song_thresh=varargin{i+1};
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1) cycle through each file in the directory
% 2) for each file run song detection
% 3) collect all seg_length long samples that are min_dist from any vocalizations
% 4) aggregate and save for analysis

% use typical song detection, power in the relevant frequency bands
% first get the file list

segsamples=round(seg_length*fs);
distsamples=round(min_dist*fs);

filelist=dir(fullfile(DIR,'*.mat'));

channel_labels=[];

for i=1:length(filelist)

	load(fullfile(DIR,filelist(i).name),'channels');

	% check for any inconsistency in channel labels

	for j=1:length(channels)

		% loop and if any channels are not included in the channel_label vector, include!

		if ~any(channels(j)==channel_labels)
			channel_labels=[channel_labels channels(j)];
		end
	end
end

channel_labels=sort(channel_labels);

disp(['Found channels ' num2str(channel_labels)]);

% initialize variables

CHANNELS=channel_labels;
EPHYS_DATA=zeros(segsamples,max_trials,length(channel_labels),'single');
MIC_DATA=zeros(segsamples,max_trials);
START_DATENUM=zeros(1,max_trials);

% as in song extract, first cycle through all files to detect which channels are present

[b,a]=butter(2,[micfilter]/(fs/2),'high');

counter=1;
contflag=1;

for i=1:length(filelist)

	load(fullfile(DIR,filelist(i).name),'ephys_data','mic_data','channels','fs','start_datenum');
	
	norm_data=filtfilt(b,a,double(mic_data));
	norm_data=norm_data./abs(max(norm_data));

	% filter the data then run through song_det with extremely low threshold	
	
	[song_bin]=song_det(norm_data,fs,minfs,maxfs,window,...
		noverlap,songduration,ratio_thresh,song_thresh);

	% where was song detected (song_bin>0)

	son_to_vec=(length(norm_data)-noverlap)/(length(song_bin));
	song_idx=find(song_bin>0);
	song_idx=unique([1 song_idx length(song_bin)]);
	song_idx=round(song_idx*son_to_vec);

	% now get points that are more than min_dist+seg_length apart

	idxdiff=diff(song_idx);

	% this will return the intervals that start gaps long enough for extractions

	candidates=find(idxdiff>distsamples*2+segsamples);

	% return song detection indices in terms of original samples
	% start the extraction point at the beginning of the interval,
        % walk in seg_length steps until the end if possible



	for j=1:length(candidates)

		% trailing edge of the song epoch

		startpoint=song_idx(candidates(j));

		% end of the song epoch

		stoppoint=song_idx(candidates(j)+1);

		k=startpoint+distsamples+1;


		while k+segsamples-1<stoppoint
			
			silentstart=k;
			silentstop=k+segsamples-1;
		
			% collect and organize into a samples x trials x channels matrix
			
			EPHYS_DATA(:,counter,:)=ephys_data(silentstart:silentstop,:);

			% collect mic data

			MIC_DATA(:,counter)=mic_data(silentstart:silentstop);

			% timestamps for save file

			START_DATENUM(counter)=start_datenum;

			k=k+segsamples;

			counter=counter+1

			if counter>max_trials
				contflag=0;
				break;
			end
		end

		if ~contflag
			break;
		end

	end

	if ~contflag
		break;
	end
end

% delete unused 

EPHYS_DATA(:,counter:max_trials,:)=[];
MIC_DATA(:,counter:max_trials)=[];
START_DATENUM(counter:max_trials)=[];

%

if ~isempty(savedir)
	
	if ~exist(fullfile(pwd,savedir),'dir')
		mkdir(fullfile(pwd,savedir));
	end

	save(fullfile(pwd,savedir,'silencedata.mat'),'EPHYS_DATA','MIC_DATA','START_DATENUM','CHANNELS');
end
