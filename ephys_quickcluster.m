function ephys_cluster(DIR,varargin)
%extracts and aligns renditions of a template along with slaved ephys
%
%example:
%
%ephys_cluster(pwd)
%
%First the script prompts the user to create a directory or continue a
%previous run, then
%the user selects a .mat file that contains the template vocalization and is prompted to 
%draw a bounding box around the template.  All of the sound files in the same directory
%are checked for spectral similarity to the template, and the user then manually cuts clusters
%to choose the cluster of sounds similar to the template (the cluster with a mean highest score
%in the feature dimensions is the likeliest candidate).  Finally, the cluster is saved to 
%a directory specified by the user in extracted_data.mat.  The results can be visualized with
%ephys_visual_mua.m (for multi-unit data).
%
%
%	ephys_cluster(DIR,varargin)
%	
%	DIR
%	directory that contains the extracted files (default: pwd)
%
%	the following may be specified as parameter/value pairs:
%
%		fs
%		sampling rate for aligned data (25e3, default Intan)
%
%		min_f
%		lowermost frequency for template spectrogram (default: 1)
%
%		max_f
%		uppermost frequency for template spectrogram (default: 10e3)
%		
%		colors
%		colormap for template spectrogram (default: hot)
%
%		padding
%		only relevant if you are using ephys_cluster to generate a template
%		for the pipeline, this will force the standalone sound clustering
%		daemon to add a pad before and after an extraction (two element vector
%		for seconds before and after extractions, in seconds)
%
%
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

spect_thresh=.1; % deprecated, this parameter is no longer used
colors='hot';
min_f=1;
max_f=9e3;
subset='';
padding=[]; % padding that will be saved with the template, in seconds (relevant for the pipeline only)
   	    % two elements vector, both specify in seconds how much time before and after to extract
	    % e.g. [.2 .2] will extract 200 msec before and after the extraction point when clustering
	    % sounds through the pipeline
hit_thresh=.3;
downfact=4; % speed up matching
% smscore parameters, THESE MUST MATCH THE PIPELINE PARAMETERS IN EPHYS_PIPELINE.CFG, OTHERWISE
% THE FEATURE COMPUTATION BETWEEN THE TEMPLATE AND CANDIDATE SOUNDS WILL NOT
% BE APPROPRIATELY MATCHED

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SMSCORE PIPELINE PARAMETERS %%%%%%%%

n=1024;
overlap=1000;
filter_scale=10;
downsampling=5;
train_classifier=1;
source='audio';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETER COLLECTION  %%%%%%%%%%%%%%

nparams=length(varargin);

if mod(nparams,2)>0
	error('ephysPipeline:argChk','Parameters must be specified as parameter/value pairs!');
end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'spect_thresh'
			spect_thresh=varargin{i+1};
		case 'colors'
			colors=varargin{i+1};
		case 'masks'
			masks=varargin{i+1};
		case 'time_range'
			time_range=varargin{i+1};
		case 'subset'
			subset=varargin{i+1};
		case 'padding'
			padding=varargin{i+1};
		case 'lowfs'
			lowfs=varargin{i+1};
		case 'highfs'
			highfs=varargin{i+1};
		case 'source'
			source=varargin{i+1};
		case 'hit_thresh'
			hit_thresh=varargin{i+1};
		case 'downfact'
			downfact=varargin{i+1};
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TODO features to read from config file, need to make this play nice with changes to smscore...


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DIRECTORY CHECK %%%%%%%%%%%%%%%%%%%%


if nargin<1 | isempty(DIR)
	DIR=pwd;
end

prev_run_listing={};

listing=dir(fullfile(DIR));

% all embedded directories could be previous runs

for i=1:length(listing)
	if listing(i).isdir & listing(i).name(1)~='.'
		prev_run_listing{end+1}=listing(i).name;
	end
end

proc_dir=[];

% check for previous runs

if ~isempty(prev_run_listing)
	response=[];
	while isempty(response)
		response=input('Would you like to go to a (p)revious run or (c)reate a new one?  ','s');
		
		switch lower(response(1))

			case 'p'
				dir_num=menu('Which directory would you like to use?',prev_run_listing);

				if isempty(dir_num), continue; end

				dir_name=prev_run_listing{dir_num};
				proc_dir=fullfile(DIR,dir_name);

			case 'c'

			otherwise
				response=[];
		end

	end
end

% prompt the user for a directory name if necessary

if isempty(proc_dir)

	dir_name=[];

	while isempty(dir_name)

		dir_name=input('What would you like to name the new directory?  ','s');

		if exist(fullfile(DIR,dir_name),'dir')
			warning('ephysPipeline:ephysCluster:direxist','Directory exists!');
			dir_name=[];
		end

	end

	proc_dir=fullfile(DIR,[ dir_name ]);
	mkdir(proc_dir);

end

%tokens=regexp(proc_dir,'\_','split');
%template_name=tokens{1};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TEMPLATE CHECK %%%%%%%%%%%%%%%%%%%%%



% check for previously extracted templates

if ~exist(fullfile(proc_dir,'template_data.mat'),'file')

	template=select_template(fullfile(DIR),source);

	% compute the features of the template
	%
	% % generate a nice sonogram of the selected template

	[b,a]=ellip(5,.2,40,[300]/(template.fs/2),'high');

	template_fig=figure('Visible','off');
	[template_image,f,t]=pretty_sonogram(filtfilt(b,a,template.data),template.fs,'N',1024,'overlap',1000,'low',1);

	startidx=max([find(f<=min_f);1]);

	if isempty(startidx)
		startidx=1;
	end

	stopidx=min([find(f>=max_f);length(f)]);

	if isempty(stopidx)
		stopidx=length(f);
	end

	imagesc(t,f(startidx:stopidx),template_image(startidx:stopidx,:));
	set(gca,'ydir','normal');

	xlabel('Time (in s)');
	ylabel('Fs');
	colormap(colors);
	multi_fig_save(template_fig,proc_dir,'template','png');

	close([template_fig]);

	save(fullfile(proc_dir,'template_data.mat'),'template');


else
	disp('Loading stored template...');
	load(fullfile(proc_dir,'template_data.mat'),'template');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GET DIFFERENCE SCORES %%%%%%%%%%%%%%

% have we computed the difference between the template and the sound data?

skip=0;
response=[];
if exist(fullfile(proc_dir,'cluster_data.mat'),'file')
	disp('Looks like you have computed the scores before...');

	while isempty(response)
		response=input('Would you like to (r)ecompute or (s)kip to clustering?  ','s');	
		switch (lower(response))
			case 'r'
				skip=0;
			case 's'
				skip=1;
			otherwise
				response=[];
		end
	end
end

% if we haven't computed the scores, do it!

if ~skip

	% collect all of the relevant .mat files

	pre_files_to_proc=dir(fullfile(DIR,'*.mat'));

	for i=1:length(pre_files_to_proc)

		files_to_proc{i}=fullfile(DIR,pre_files_to_proc(i).name);
	end

	disp('Comparing sound files to the template (this may take a minute)...');

	% take a subset if the user has passed the option
	%
	
	if length(subset)==1
		disp(['Will use ' num2str(subset*100) '% of the available files']);
		%selection=randsample(1:length(files_to_proc),floor(length(files_to_proc)*subset));
		%selection=sort(selection);
		selection=round(linspace(1,length(files_to_proc),...
			floor(length(files_to_proc)*subset)));
		files_to_proc=files_to_proc(selection);
	elseif length(subset>1)
		disp('Will use the user provided subset');

		subset(subset>length(files_to_proc))=[];
		files_to_proc=files_to_proc(subset);
	end

	template_match(template,files_to_proc,fullfile(proc_dir,'cluster_data.mat'),source,hit_thresh,downfact);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CLUSTERING GUI %%%%%%%%%%%%%%%%%%%%%

% do we need to cluster again?

load(fullfile(proc_dir,'cluster_data.mat'),'filenames','hit_locs');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% HIT EXTRACTION %%%%%%%%%%%%%%%%%%%%%


skip=0;
response=[];
if exist(fullfile(proc_dir,'extracted_data.mat'),'file')
	disp('Looks like you have extracted the data before..');

	while isempty(response)
		response=input('Would you like to (r)eextract or (s)kip?  ','s');	
		switch (lower(response))
			case 'r'
				skip=0;
			case 's'
				skip=1;
			otherwise
				response=[];
		end
	end
end


if ~skip
	disp(['Saving data to ' proc_dir]);
	extract_hits(hit_locs,filenames,length(template.data),padding,proc_dir,dir_name);
end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TEMPLATE SELECT %%%%%%%%%%%%%%%%%%%%

function [TEMPLATE]=select_template(DIR,SOURCE)

pause(.001); % inserting 1 msec pause since uigetfile does not always open without it, not sure why...

response=[];

while isempty(response)
	[filename,pathname]=uigetfile('*.mat','Pick a sound file to extract the template from',fullfile(DIR));
	
	is_legacy=check_legacy(fullfile(pathname,filename));

	if is_legacy
		load(fullfile(pathname,filename),'mic_data','fs');
		audio.data=mic_data;
		audio.fs=fs;

		clearvars mic_data fs;
	else
		load(fullfile(pathname,filename),SOURCE);
	end

	switch lower(SOURCE(1))
		case 'a'
			tmp=audio;
		case 'p'
			tmp=playback;
	end

	TEMPLATE.data=spectro_navigate(tmp.data);
	TEMPLATE.fs=tmp.fs;

	response2=[];
	while isempty(response2)
		
		response2=input('(C)ontinue with selected template or (s)elect another sound file?  ','s');

		switch lower(response2(1))
			case 'c'
				response=1;
			case 's'
				response=[];
			otherwise
				response2=[];
		end

	end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% COMPARE FEATURES WITH TEMPLATE %%%%%


function template_match(TEMPLATE,TARGET_FILES,SAVEFILE,SOURCE,THRESH,DOWNFACT)

% do the template matching here...

% simple matched filter, don't worry about covariance term

template_filt=downsample(detrend(TEMPLATE.data(end:-1:1)),DOWNFACT);

% normalize convolution by autocorrelation (limit to max 1 for perfect match)

template_normfac=max(xcorr(template_filt,100));
hit_locs={};

[nblanks formatstring]=progressbar(100);
fprintf(1,['Progress:  ' blanks(nblanks)]);

for i=1:length(TARGET_FILES)

	% load the features of the sound data

	fprintf(1,formatstring,round((i/length(TARGET_FILES))*100));
	
	target=[];

	[path,name,ext]=fileparts(TARGET_FILES{i});

	load(TARGET_FILES{i},SOURCE);

	switch lower(SOURCE(1))
		case 'a'
			target_data=audio;
		case 'p'
			target_data=playback;
		otherwise
	end

	score=conv(downsample(detrend(target_data.data),DOWNFACT),template_filt,'same')./template_normfac;
	figure(1);plot(score);
	warning('off','signal:findpeaks:largeMinPeakHeight');
	[~,hits]=findpeaks(score,'minpeakheight',THRESH,'minpeakdistance',100);
	warning('on','signal:findpeaks:largeMinPeakHeight');
	hit_locs{i}=hits*DOWNFACT;

end

fprintf(1,'\n');

filenames=TARGET_FILES;
save(SAVEFILE,'hit_locs','filenames');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%% the grand finale, extract the data!

% need to adapt to grab aligned sound data in a sample x trials matrix
% and a cell array of matrices for the Intan data (aligned for each electrodes)
%
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DATA EXTRACTION %%%%%%%%%%%%%%%%%%%%


function extract_hits(SELECTED_PEAKS,FILENAMES,TEMPLATESIZE,PADDING,OUT_DIR,EXTNAME)

USED_FILENAMES={};

disp(['Extracting cluster (' num2str(sum(cellfun(@length,SELECTED_PEAKS))) ' peaks):   ']);

[nblanks formatstring]=progressbar(100);
counter=0;

% check for all possible channels across the whole day, a matrix will be filled with zeros if the channel
% gets knocked out somehow...

tempsize=round(TEMPLATESIZE/2);

image_dir=fullfile(OUT_DIR,'gif');
wav_dir=fullfile(OUT_DIR,'wav');
data_dir=fullfile(OUT_DIR,'mat');

mkdir(image_dir);
mkdir(wav_dir);
mkdir(data_dir);

disp_minfs=0;
disp_maxfs=9e3;
colors='hot';

dirstruct=struct('image',image_dir,'wav',wav_dir,'data',data_dir);

for i=1:length(SELECTED_PEAKS)

	datastruct=load(FILENAMES{i});
	[~,filename,~]=fileparts(FILENAMES{i});

	if length(SELECTED_PEAKS{i})<1
		continue;
	end

	ext_pts=zeros(length(SELECTED_PEAKS{i}),2);

	for j=1:length(SELECTED_PEAKS{i})
		
		% padding?
		
		ext_pts(j,1)=(SELECTED_PEAKS{i}(j)-tempsize)/datastruct.audio.fs;
		ext_pts(j,2)=(SELECTED_PEAKS{i}(j)+tempsize)/datastruct.audio.fs;

		% get all the data types save, write out gif, mat, wav etc.

	end

	frontend_dataextract(filename,datastruct,dirstruct,ext_pts,...
		disp_minfs,disp_maxfs,colors,'playback',1,'',[ '_' EXTNAME ],1);
	frontend_dataextract(filename,datastruct,dirstruct,ext_pts,...
		disp_minfs,disp_maxfs,colors,'audio',0,'',[ '_' EXTNAME ],1);

end

end

