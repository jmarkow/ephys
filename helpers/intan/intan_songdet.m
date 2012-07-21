function intan_songdet(DIR,INTAN_DIR,varargin)
%intan_songdet.m processes .mat files for singing by looking at the ratio
%between energy in a particular "song" band and energy at all other frequencies.
%If song is detected, the script searches for Intan files with related ephys data
%according to their time stamp and performs a fine alignment by syncing TTL sequences
%aligned to the microphone recording and ephys recording (via the Intan aux channel).
%The output is written to "extracted_song" in the current folder.  
%
%Examples:
%
%intan_songdet(pwd)
%
%will use .mat files in the current directory and call up a GUI to choose the directory
%with the Intan files for alignment, the output will be saved in "extracted_song" in the 
%current directory.
%
%
%	intan_songdet(DIR,INTAN_DIR,varargin)
%
%	DIR
%	directory that contains .mat files with the microphone data (default: pwd)
%
%	INTAN_DIR
%	directory with Intan files (leaving empty will call up a GUI to choose a dir)
%
%	the following may be specified as parameter/value pairs:
%
%		mic_trace
%		column in data.voltage that contains the microphone trace
%
%		ttl_trace
%		column in data.voltage that contains the TTL trace
%		
%		minfs
%		lowermost frequency in the song band
%
%		maxfs
%		uppermost frequency in the song band (must be >minfs)
%
%		ratio_tresh
%		threshold for song energy ratio (default:  3.5)
%
%		window
%		window to compute song energy ratio (default: 250 samples)
%
%		noverlap
%		overlap for computing song energy ratio (default: 0 samples)
%
%		song_thresh
%		threshold for song detection (default:  .3)
%
%		song_duration
%		smoothing duration for song energy ratio (default: 800 ms)
%
%		colors
%		colormap to use for generating sonograms (default:  hot)
%
%		intan_interval
%		save interval for Intan files (default: 60s, Intan default)
%
%		intan_fs
%		fs of Intan recording (default: 25e3, Intan default)
%
%		aux_trace
%		aux channel used for TTL sync signal (default: 6)
%
%		ephys_preview
%		define to write image files with aligned voltages (default: empty,
%		define to enable)
%
%		ephys_exclude
%		exclude electrodes for ephys_preview (default: empty)
%
%		debug
%		define to display windows with all relevant alignment information (default: empty)
%
%
%see also intan_align.m,read_intan_data_cli.m
%

% where are the Intan files?

if nargin<2 | isempty(INTAN_DIR)
	INTAN_DIR=uigetdir(pwd,'Select the directory with matching Intan files...');
end

if nargin<1 | isempty(DIR), DIR=pwd; end

nparams=length(varargin);

if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs!');
end

mic_trace=1;
ttl_trace=2;
logfile=fullfile(DIR,'log.txt');

minfs=1e3;
maxfs=4e3;
ratio_thresh=2.5;
window=250;
noverlap=0;
song_thresh=.5; % between .2 and .3 seems to work best (higher is more exlusive)
songduration=.5;
low=5;
high=10;
colors=hot;
intan_lag=1; %consistent lag on the Intan TTL 
intan_pad=2; %padding in seconds for the Intan extraction to account for sloppy clocks
	     %1-2 seems to be more than enough

intan_interval=60; %intan save interval in seconds
intan_fs=25e3; %sampling rate of the Intan chip
audio_fs=40e3; %audio sampling rate
audio_pad=.2; % how much padding on each side for the extractions
aux_trace=6; %which Intan TTL trace to process
ephys_preview=[]; %write out files with squared & smoothed voltage traces
ephys_exclude=[];
debug=[];

for i=1:2:nparams
	switch lower(varargin{i})
		case 'mic_trace'
			mic_trace=varargin{i+1};
		case 'ttl_trace'
			ttl_trace=varargin{i+1};
		case 'logfile'
			logfile=varargin{i+1};
		case 'window'
			window=varargin{i+1};
		case 'noverlap'
			noverlap=varargin{i+1};
		case 'minfs'
			minfs=varargin{i+1};
		case 'maxfs'
			maxfs=varargin{i+1};
		case 'ratio_thresh'
			ratio_thresh=varargin{i+1};
		case 'song_thresh'
			song_thresh=varargin{i+1};
		case 'colors'
			colors=varargin{i+1};
		case 'songduration'
			songduration=varargin{i+1};
		case 'intan_interval'
			intan_interval=varargin{i+1};
		case 'intan_fs'
			intan_fs=varargin{i+1};
		case 'intan_interval'
			intan_interval=varargin{i+1};
		case 'intan_pad'
			intan_pad=varargin{i+1};
		case 'ephys_preview'
			ephys_preview=varargin{i+1};
		case 'ephys_exclude'
			ephys_exclude=varargin{i+1};
		case 'debug'
			debug=varargin{i+1};
		otherwise
	end
end



dir_name=[];

while isempty(dir_name)

	dir_name=input('What would you like to name the directory to store the results?  ','s');

	if exist(fullfile(DIR,dir_name),'dir')
		disp('Directory exists!');

		response=[];

		while isempty(response)

			if exist(fullfile(DIR,dir_name,'parameters.mat'),'file')
				old=load(fullfile(DIR,dir_name,'parameters.mat'),...
					'window','noverlap','ratio_thresh','song_thresh',...
					'minfs','maxfs','songduration');
				flag=(old.window==window)&&(old.noverlap==noverlap)&&(old.ratio_thresh==ratio_thresh)...
					&&(old.song_thresh==song_thresh)&&(old.minfs==minfs)...
					&&(old.maxfs==maxfs)&&(old.songduration==songduration);
				if flag
					disp('Parameters match');
				else
					disp('Parameters do not match');
				end


			end

			response=input('Okay to proceed, files may be overwritten (y or n)?  ','s');

			switch lower(response(1))

				case 'y'
					break;
				case 'n'
					dir_name=[];
				otherwise
					response=[];
			end
		end

	end

end

mkdir(fullfile(DIR,dir_name));

save(fullfile(DIR,dir_name,'parameters.mat'),...
	'window','noverlap','ratio_thresh','song_thresh','minfs','maxfs','songduration');

save_dir=fullfile(DIR,dir_name,'mat');
image_dir=fullfile(DIR,dir_name,'gif');
wav_dir=fullfile(DIR,dir_name,'wav');
mkdir(save_dir);
mkdir(image_dir);
mkdir(wav_dir);

% basic idea:
%
% 1) process mic trace for singing
% 2) extract slaved TTL trace
% 3) search the directory with relevant Intan data and parse the timestamp 
%    in the filenames
% 4) use the timestamp from the mic trace to identiy the point in the Intan
%    data to grab
% 5) extract multi-channel Intan data and use the aux TTL for precise alignment

% standard anonymous function to be used in parfor for saving
% files 

% should add an option to create separate directories, standard query...

par_save = @(FILE,data) save([FILE],'data');

%%

% it is assumed that the TTL trace is the column after the mic trace, if not, specify 
% as a parameter/value pair

if isempty(ttl_trace)
	ttl_trace=mic_trace+1;
end

% the mic traces should be stored in the mat files

if isempty(DIR)
	[filename,pathname]=uigetfile('*.mat','MAT-files (*.mat)',...
		'Pick .mat files with mic and TTL traces','MultiSelect','on');
	for i=1:length(filename)
		files_to_proc{i}=fullfile(pathname,filename{i});
	end
else
	files_to_proc_pre=dir(fullfile(DIR,'*.mat'));
	files_to_proc={files_to_proc_pre(:).name};
end
	
% need to filter out the low_fs noise in the mic recording

match_score=[];

[b,a]=butter(5,[600/(audio_fs/2)],'high');

parfor i=1:length(files_to_proc)

	ax=[];
	data=[];

	disp(['Processing ' files_to_proc{i}]);

	% now we need to pull out time and extract the appropriate data from Intan

	try
		data=getfield(load(files_to_proc{i},'data'),'data');
	catch	
		continue;
	end

	fs=data.sampling_rate;

	%[b,a]=butter(3,(600/(fs/2)),'high');

	audio_data=data.voltage(:,mic_trace);

	if length(audio_data)<window
		continue;
	end

	% run song detection, check for power in the songband relative to power outside

	song_bin=song_det(audio_data,fs,minfs,maxfs,window,noverlap,songduration,ratio_thresh,song_thresh);

	song_pts=find(song_bin>0);

	if isempty(song_pts)
		continue;
	else
		disp(['Song detected in file:  ' files_to_proc{i}]);
	end

	% create a large sonogram with detected song marked with a black bar

	audio_data=filtfilt(b,a,audio_data);

	[sonogram_im,F,T]=spectrogram(audio_data,window,noverlap,[]);
	sonogram_im=1e4*abs(flipdim(sonogram_im,1));
	sonogram_im=log(sonogram_im+2);
	sonogram_im(sonogram_im>high)=high;
	sonogram_im(sonogram_im<low)=low;
	sonogram_im=sonogram_im-low;
	sonogram_im=sonogram_im/(high-low);
	sonogram_im=63*(1-sonogram_im);

	% factor to move from sonogram coordinates to raw audio data coordinates

	son_to_vec=(length(audio_data)-noverlap)/(length(song_bin));

	% use diff to find non_continguous song bouts

	song_idx=[0 find(diff(song_pts*son_to_vec)>fs) length(song_pts)];

	[path,name,ext]=fileparts(files_to_proc{i});

	sonogram_filename=fullfile(image_dir,[ name '.gif' ]);

	for j=1:length(song_idx)-1

		startpoint=floor((song_pts(song_idx(j)+1))*son_to_vec-audio_pad*fs);
		endpoint=ceil((song_pts(song_idx(j+1)))*son_to_vec+audio_pad*fs);

		if startpoint<1, startpoint=1; end
		if endpoint>length(audio_data), endpoint=length(audio_data); end

		audio_extraction=audio_data(startpoint:endpoint);
		save_name=[ name '_chunk_' num2str(j) ];
		
		if length(audio_extraction)<500
			continue;
		end

		if exist(fullfile(save_dir,[ save_name '_aggregated.mat']),'file')
			disp([ save_name '_aggregated.mat' ' already processed, skipping...']);
			continue;
		end

		ttl_extraction=double(data.voltage(startpoint:endpoint,ttl_trace)>1);

		sonogram_im(1:10,ceil(startpoint/son_to_vec):ceil(endpoint/son_to_vec))=0;


		[chunk_sonogram_im,F,T]=spectrogram(audio_extraction,500,350,[]);
		chunk_sonogram_im=1e4*abs(flipdim(chunk_sonogram_im,1));
		chunk_sonogram_im=log(chunk_sonogram_im+2);
		chunk_sonogram_im(chunk_sonogram_im>high)=high;
		chunk_sonogram_im(chunk_sonogram_im<low)=low;
		chunk_sonogram_im=chunk_sonogram_im-low;
		chunk_sonogram_im=chunk_sonogram_im/(high-low);
		chunk_sonogram_im=63*(1-chunk_sonogram_im);

		imwrite(uint8(chunk_sonogram_im),colors,fullfile(image_dir,[ save_name '.gif']),'gif');
				
		% the actual time our extraction starts

		starting_time=addtodate(datenum(data.start_time),round(data.time(startpoint))+intan_lag,'second');

		% resample data acquired by the nidaq to match Intan data

		resample_factor=intan_fs/fs;
		[p,q]=rat(resample_factor);

		res_mic_data=resample(audio_extraction,p,q);
		res_mic_ttl=double(resample(ttl_extraction,p,q)>.5);

		sound_len=length(res_mic_data)-1; % length of the sound in samples

		[intan_extraction,original_ttl,intan_ttl,channels,match_score_tmp]=intan_extract(INTAN_DIR,intan_interval,intan_fs,aux_trace,...
			res_mic_ttl,starting_time,sound_len,intan_pad);
		
		% write the high fs audio data out as a .wav

		min_audio=min(audio_extraction(:));
		max_audio=max(audio_extraction(:));

		if min_audio + max_audio < 0
			audio_extraction=audio_extraction./(-min_audio);
		else
			audio_extraction=audio_extraction./(max_audio*(1+1e-3));
		end
	
		%disp([num2str(min(audio_extraction))]);
		%disp([num2str(max(audio_extraction))]);

		wavwrite(audio_extraction,fs,fullfile(wav_dir,[ save_name '.wav']));
		
		% write out the resampled audio data and ttl aligned to Intan data

		if ~isempty(intan_extraction)

			parsave(fullfile(save_dir,[ save_name '_aggregated.mat']),...
				intan_extraction,original_ttl,intan_ttl,res_mic_ttl,res_mic_data,intan_fs,channels);

			% write out squared and smoothed traces to a PDF for previewing

			if ~isempty(ephys_preview)
				write_ephys_preview(fullfile(save_dir,[ save_name '_ephys']),...
					intan_extraction,res_mic_data,intan_fs,ephys_exclude);
			end

			match_score=[match_score match_score_tmp];

			if ~isempty(debug)
				figure();

				ax(1)=subplot(3,1,1);plot(res_mic_ttl);
				ylabel('NiDAQ TTL');
				ax(2)=subplot(3,1,2);plot(original_ttl);
				ylabel('Intan guess');
				ax(3)=subplot(3,1,3);plot(intan_ttl);
				ylabel('Adjusted with xcorr');

				linkaxes([ax(1) ax(3)],'x');
			end



		end

	end

	reformatted_im=im_reformat(sonogram_im,10);
	imwrite(uint8(reformatted_im),colors,sonogram_filename,'gif');

end

save(fullfile(DIR,dir_name,'match_scores.mat'),'match_score');

end

%%%%%%%%%%%%%%%

function song_idx=song_det(audio,fs,minfs,maxfs,window,noverlap,songduration,ratio_thresh,song_thresh)

[s,f,t]=spectrogram(audio,window,noverlap,[],fs);

% take the power and find our fs band

power=abs(s);
min_idx=max(find(f<=minfs));
max_idx=min(find(f>=maxfs));

% take the song/nonsong power ratio

song=mean(power(min_idx:max_idx,:),1);
nonsong=mean(power([1:min_idx-1 max_idx+1:end],:),1)+eps;

song_ratio=song./nonsong;
%song_det=smooth(double(song_ratio>ratio_thresh),window);

% convolve with a moving average filter

filt_size=round((fs*songduration)/(window-noverlap));
mov_filt=ones(1,filt_size)*1/filt_size;
song_detvec=conv(double(song_ratio>ratio_thresh),mov_filt,'same');

% where is the threshold exceeded?

song_idx=song_detvec>song_thresh;

end

%%%%%%%%%%%%%%

function [intan_extraction,original_ttl,intan_ttl,channels,match_score]=intan_extract(intan_dir,interval,intan_fs,aux_trace,...
		res_mic_ttl,starting_time,sound_len,padding)

% time is passed as a datenum, we need to find the Intan
% file that recorded data at the same time

% parse the intan filename, which is a stamp of the TIME 
% the file is opened, writing appears to have a slight
% delay (thus making the Intan data 1-2 s delayed relative
% to nidaq according to the file's timestamp)

% and calculate approximately which data we want to extract

% thus the time has to be within the interval TIME_STAMP-save_interval

pre_file_list=dir(fullfile(intan_dir,'*.int')); %intan files have the extension int

file_list={};

for i=1:length(pre_file_list)
	file_list{i}=fullfile(intan_dir,pre_file_list(i).name);
end

intan_extraction=[];
original_ttl=[];
intan_ttl=[];
channels=[];
match_score=[];

% should convert the Intan grab to samples not seconds (rounding errors are causing the 
% the Intan data to have an extra sample occassionally, just specify in samples!)

for i=1:length(file_list)

	[path,name,ext]=fileparts(file_list{i});
	tokens=regexp(name,'\_','split');

	% the second token is the date, third the time

	intan_start=datenum([tokens{end-1} tokens{end}],'yymmddHHMMSS');

	% can we grab the entire sound with the padding?

	% intan_start<micdata_start-padding & micdata_start+sound_len+padding<intan_start+save interval (time
	% collected in each Intan file)


	flag_1=addtodate(starting_time,-padding,'second')>intan_start;
	flag_2=addtodate(starting_time,ceil(sound_len/intan_fs)+padding,'second')<addtodate(intan_start,interval,'second');

	if flag_1 && flag_2

		% estimated offset is the time elapsed between the starting time of the imc
		% trace and the start of the Intan file

		offset=etime(datevec(starting_time),datevec(intan_start));
		disp(['Checking file: '  file_list{i}]);
		disp(['Assuming offset:  ' num2str(offset)]); 
		[t,channels,intan_data,aux] = read_intan_data_cli(file_list{i});

		% extract at the estimated point-fudge factor in samples

		startpoint=round((offset-padding)*intan_fs);

		% the initial extraction (endpoint+fudge factor in secs, conv to samples)

		stoppoint=startpoint+sound_len+(padding*intan_fs);

		if stoppoint>length(t)
			continue;
		end

		intan_ttl=double(aux(startpoint:stoppoint,6));

		original_ttl=intan_ttl;

		disp('Adjusting alignment...');

		%intan_ttl(intan_ttl==0)=-1;
		[align_score,lags]=xcorr(intan_ttl,res_mic_ttl);
		[dummy,delay]=max(align_score);
		
		align_point=lags(delay);
		startpoint=align_point+startpoint;

		disp([num2str([length(intan_ttl) length(res_mic_ttl) align_point])]);

		[samples,traces]=size(intan_data);

		if startpoint<1 || (startpoint+sound_len)>samples
			disp('Lag too short or too large, skipping alignment...');
			return;
		end

		intan_ttl=double(aux(startpoint:startpoint+sound_len,aux_trace));
		intan_extraction=intan_data(startpoint:startpoint+sound_len,:);

		match_intan=intan_ttl;
		match_intan(match_intan==0)=-1;
		match_nidaq=res_mic_ttl;
		match_nidaq(match_nidaq==0)=-1;

		match_score=sum(match_intan.*match_nidaq)/length(match_intan);
		disp(['Match score ' num2str(match_score)]);

		if match_score<.75
			disp('Warning!  Match score is less than expected, skipping alignment...');
			intan_ttl=[];
			intan_extraction=[];
		end

		return;

	end

end


end

%%%%%%%%%%%%%%%

function parsave(file,ephys_data,original_intan_ttl,adjusted_intan_ttl,resampled_mic_ttl,mic_data,fs,channels)

save(file,'ephys_data','original_intan_ttl','adjusted_intan_ttl','resampled_mic_ttl','mic_data','fs','channels');

end


function write_ephys_preview(save_file,intan_data,mic_data,fs,exclude)

% spectrogram parameters

low=6;
high=13;
smoothing=.01; % smoothing window in seconds
subplot_space=.05;

if nargin<5, exclude=[]; end

% for display just square and smooth after car noise rejection

good_electrodes=setdiff([1:16],exclude);
nplots=length(good_electrodes)+3;

[samples,traces]=size(intan_data);

% figure dimensions

width=samples/10.5e3
height=nplots*.9


% make a large subplot with nelectrodes abutting subplots

[mic_sonogram,F,T]=spectrogram(mic_data,500,350,[],fs,'yaxis');

mic_sonogram=1e4*abs(mic_sonogram);
mic_sonogram=log(mic_sonogram+2);
mic_sonogram(mic_sonogram>high)=high;
mic_sonogram(mic_sonogram<low)=low;
mic_sonogram=mic_sonogram-low;
mic_sonogram=mic_sonogram/(high-low);
mic_sonogram=64*(1-mic_sonogram);

all_traces_fig=figure('Visible','off');
ax(1)=subplot(nplots,1,1:3);
image(T,F,mic_sonogram);colormap(hot);
set(ax(1),'ydir','norm');
ylabel('Hz');
box off;
set(gca,'xtick',[],'tickdir','out');
hold on;

for i=1:length(good_electrodes)
	
	ax(i+1)=subplot(nplots,1,i+3);

	% car rejection

	noise_estimate=mean(intan_data(:,good_electrodes),2);

	% square and smooth

	plot_data=smooth((intan_data(:,good_electrodes(i))-noise_estimate).^2,smoothing*fs);
	plot([1:samples]./fs,plot_data);
	ylabel(['Elec. ' num2str(good_electrodes(i))]);

	%if i>1
	%	prev_position=get(subplot(nplots,1,i+3-1),'position');
	%	curr_position=get(ax(i+1),'position');
	%	new_bot=curr_position(3)+subplot_space;
	%	curr_position(3)=new_bot;
	%	set(ax(i+1),'position',curr_position);
	%end
		
	%set(gca,'xtick',[],'ytick',[min(plot_data) max(plot_data)],'tickdir','out');
	box off
	
	

end

axis tight
linkaxes(ax,'x');

for i=1:length(good_electrodes)-1
	subplot(nplots,1,i+3);
	ylimits=ylim();
	set(gca,'xtick',[],'tickdir','out');
end

subplot(nplots,1,nplots);
xlabel('Time (in seconds');

set(all_traces_fig,'PaperSize',[width height])
set(all_traces_fig,'PaperPosition',[ .1 .1 width-.2 height-.2 ]);
%orient
%orient(all_traces_fig,'portrait');

print(all_traces_fig,'-depsc2',fullfile([save_file '.eps']));
close([all_traces_fig]);

end
