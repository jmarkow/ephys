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
%		n
%		spectral feature score spectrogram window, if you are using the pipeline
%		this MUST match the pipeline parameters in ephys_pipeline.cfg (default: 1024)
%
%		overlap
%		spectral feature score spectrogram overlap, must match ephys_pipeline.cfg (default: 1000)
%
%		filter_scale
%		spectral feature score smoothing window size, must match ephys_pipeline.cfg (default: 10)
%
%		downsampling
%		spectral feature downsampling factor, must match ephys_pipeline.cfg (default: 5)
%
%
%
%see also songdet.m,ephys_visual_mua.m,ephys_visual_sua.m,ephys_pipeline_smscore.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

spect_thresh=.1; % deprecated, this parameter is no longer used
fs=25e3; % sampling
colors='hot';
min_f=1;
max_f=10e3;
time_range=[0 inf];
subset='';
padding=[]; % padding that will be saved with the template, in seconds (relevant for the pipeline only)
   	    % two elements vector, both specify in seconds how much time before and after to extract
	    % e.g. [.2 .2] will extract 200 msec before and after the extraction point when clustering
	    % sounds through the pipeline

% smscore parameters, THESE MUST MATCH THE PIPELINE PARAMETERS IN EPHYS_PIPELINE.CFG, OTHERWISE
% THE FEATURE COMPUTATION BETWEEN THE TEMPLATE AND CANDIDATE SOUNDS WILL NOT
% BE APPROPRIATELY MATCHED

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SMSCORE PIPELINE PARAMETERS %%%%%%%%

n=1024;
overlap=1000;
filter_scale=10;
downsampling=5;

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
		case 'fs'
			fs=varargin{i+1};
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
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TODO features to read from config file, need to make this play nice with changes to smscore...


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DIRECTORY CHECK %%%%%%%%%%%%%%%%%%%%


if nargin<1 | isempty(DIR)
	DIR=uigetdir(pwd,'Select the directory with .mat files to process...');
end

prev_run_listing={};

listing=dir(fullfile(DIR));

% all embedded directories could be previous runs

for i=1:length(listing)
	if listing(i).isdir
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

	proc_dir=fullfile(DIR,[ dir_name '_MANUALCLUST']);
	mkdir(proc_dir);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TEMPLATE CHECK %%%%%%%%%%%%%%%%%%%%%



% check for previously extracted templates

if ~exist(fullfile(proc_dir,'template_data.mat'),'file')

	TEMPLATE=select_template(fullfile(DIR));

	% compute the features of the template

	disp('Computing the spectral features of the template');
	template_features=ephys_pipeline_smscore(TEMPLATE,fs,...
		'n',n,'overlap',overlap,'filter_scale',filter_scale,'downsampling',downsampling);
	save(fullfile(proc_dir,'template_data.mat'),'TEMPLATE','template_features','padding');

else
	disp('Loading stored template...');
	load(fullfile(proc_dir,'template_data.mat'),'TEMPLATE');
	
	disp('Computing the spectral features of the template');
	template_features=ephys_pipeline_smscore(TEMPLATE,fs,...
		'n',n,'overlap',overlap,'filter_scale',filter_scale,'downsampling',downsampling);
	save(fullfile(proc_dir,'template_data.mat'),'TEMPLATE','template_features','padding');

end

% generate a nice sonogram of the selected template

template_fig=figure('Visible','off');
[template_image,f,t]=pretty_sonogram(TEMPLATE,fs,'N',1024,'overlap',1000);

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

% get the template size so we can extract hits of the same size

[junk,templength]=size(template_features{1});
templength=templength-1;

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

		if ~isempty(findstr('match_scores.mat',pre_files_to_proc(i).name))
			continue;
		end

		files_to_proc{i}=fullfile(DIR,pre_files_to_proc(i).name);
	end

	disp('Computing features for all sounds...');

	sound_file_features(DIR,files_to_proc,n,overlap,filter_scale,downsampling);

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

	template_match(template_features,files_to_proc,fullfile(proc_dir,'cluster_data.mat'),templength);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CLUSTERING GUI %%%%%%%%%%%%%%%%%%%%%

% TODO update cluster gui so that it has parity with the spike sorting GUI (more advanced)


property_names={'cos','derivx', 'derivy', 'amp','product','curvature'};
save(fullfile(proc_dir,'cluster_data.mat'),'property_names','proc_dir','-append');

% do we need to cluster again?

skip=0;
response=[];
if exist(fullfile(proc_dir,'cluster_results.mat'),'file')
	disp('Looks like you have clustered the data before..');

	while isempty(response)
		response=input('Would you like to (r)ecluster or (s)kip?  ','s');	
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
	uiwait(new_data_plotter(fullfile(proc_dir,'cluster_data.mat'),fullfile(proc_dir,'cluster_results.mat')));
end

load(fullfile(proc_dir,'cluster_results.mat'),'sorted_syllable');
load(fullfile(proc_dir,'cluster_data.mat'),'filenames');

act_templatesize=length(TEMPLATE);

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
	[mic_data ephys_data channels used_filenames]=extract_hits(sorted_syllable,filenames,...
		act_templatesize,spect_thresh,time_range,fs,n,overlap,downsampling);

	disp(['Saving data to ' fullfile(proc_dir,'extracted_data.mat')]);

	save(fullfile(proc_dir,'extracted_data.mat'),'used_filenames','mic_data','ephys_data','time_range','channels','-v7.3');
end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TEMPLATE SELECT %%%%%%%%%%%%%%%%%%%%

function TEMPLATE=select_template(DIR)

pause(.001); % inserting 1 msec pause since uigetfile does not always open without it, not sure why...

response=[];

while isempty(response)
	[filename,pathname]=uigetfile('*.mat','Pick a sound file to extract the template from',fullfile(DIR));
	load(fullfile(pathname,filename),'mic_data');
	TEMPLATE=spectro_navigate(mic_data);

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

% function to compute the spectral features for all the pertinent wav files

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% COMPUTE FEATURES %%%%%%%%%%%%%%%%%%%


function sound_file_features(DIR,SOUND_FILES,N,OVERLAP,FILTER_SCALE,DOWNSAMPLING)

par_save = @(FILE,features) save([FILE],'features');

if ~exist(fullfile(DIR,'syllable_data'),'dir')
	mkdir(fullfile(DIR,'syllable_data'));
end

parfor i=1:length(SOUND_FILES)

	input_file=SOUND_FILES{i};
	disp([input_file])

	[path,name,ext]=fileparts(input_file);
	output_file=fullfile(DIR,'syllable_data',[ name '_score.mat']);

	if exist(output_file,'file'), continue; end

	disp(['Computing features for ' input_file]);

	% simply read in the file and score it

	% getfield hack to get around parfor errors

	data=load(input_file,'mic_data','fs');
	if ~isfield(data,'mic_data')
		warning('ephysPipeline:ephysCluster:errorsoundfile','Problem encountered with %s',input_file);
		continue;
	end
	
	sound_data=data.mic_data;
	fs=data.fs;

	if length(sound_data)<N
		warning('ephysPipeline:ephysCluster:shortsound','Sound extraction too short in %s, skipping...',input_file);
		continue;
	end

	sound_features=ephys_pipeline_smscore(sound_data,fs,...
		'n',N,'overlap',OVERLAP,'filter_scale',FILTER_SCALE,'downsampling',DOWNSAMPLING);

	% save for posterity's sake

	par_save(output_file,sound_features);

end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% COMPARE FEATURES WITH TEMPLATE %%%%%


function template_match(TEMPLATE,TARGET_FILES,SAVEFILE,TEMPLATESIZE)

% do the template matching here...

%disp('Comparing the target sounds to the template...');


parfor i=1:length(TARGET_FILES)

	% load the features of the sound data

	target=[];

	[path,name,ext]=fileparts(TARGET_FILES{i});

	input_file=fullfile(path,'syllable_data',[ name '_score.mat']);	
	
	try
		target=getfield(load(input_file,'features'),'features');
	catch
		warning('ephysPipeline:ephysCluster:errortemplatefile','Troubling reading %s',input_file);
		continue;
	end

	[junk,targetlength]=size(target{1});

	score_temp={};
	temp_mat=[];

	for j=1:length(target)
		score_temp{j}=[];

		for k=1:targetlength-TEMPLATESIZE
			score_temp{j}=[score_temp{j} sum(sum(abs(target{j}(:,k:k+TEMPLATESIZE)-TEMPLATE{j})))];
		end

		score_temp{j}=score_temp{j}-mean(score_temp{j});
		score_temp{j}=score_temp{j}/std(score_temp{j});
		score_temp{j}(score_temp{j}>0)=0;
		score_temp{j}=abs(score_temp{j});

	end

	attributes=length(score_temp);
	product_score=score_temp{1};
	
	for j=2:attributes, product_score=product_score.*score_temp{j}; end

	if length(product_score)<3
		variableCellArray{i}=temp_mat;
		peakLocation{i}=[];
		continue;
	end
	
	warning('off','signal:findpeaks:largeMinPeakHeight');
	[pks,locs]=findpeaks(product_score,'MINPEAKHEIGHT',.005);
	warning('on','signal:findpeaks:largeMinPeakHeight');
	
	if isempty(locs)
		variableCellArray{i}=temp_mat;
		peakLocation{i}=[];
		continue; 
	end

	curvature=gradient(gradient(product_score));

	for j=1:attributes, temp_mat(:,j)=log(score_temp{j}(locs)); end

	temp_mat(:,attributes+1)=log(product_score(locs));
	temp_mat(:,attributes+2)=log(abs(curvature(locs)));

	peakLocation{i}=locs;
	variableCellArray{i}=temp_mat;

end

warning('on','signal:findpeaks:largeMinPeakHeight');
filenames=TARGET_FILES;

empty_coords=find(cellfun(@isempty,variableCellArray));
variableCellArray(empty_coords)=[];
peakLocation(empty_coords)=[];
filenames(empty_coords)=[];

disp([length(variableCellArray)]);

save(SAVEFILE,'variableCellArray','peakLocation','filenames');

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


function [MIC_DATA EPHYS_DATA CHANNELS USED_FILENAMES]=extract_hits(SELECTED_PEAKS,FILENAMES,TEMPLATESIZE,SPECT_THRESH,TIME_RANGE,fs,N,OVERLAP,DOWNSAMPLING)

TEMPLATESIZE=TEMPLATESIZE+N;
USED_FILENAMES={};
MIC_DATA=[];

load(FILENAMES{1});

if ~exist('channels','var');
	[samples,nchannels]=size(ephys_data);
	CHANNELS=1:nchannels;
else
	CHANNELS=channels;
end

disp(['Extracting cluster (' num2str(length(SELECTED_PEAKS)) ' peaks):   ']);
disp('Preallocating matrices (this may take a minute)...');

[nblanks formatstring]=progressbar(100);
counter=0;

% check for all possible channels across the whole day, a matrix will be filled with zeros if the channel
% gets knocked out somehow...

channel_labels=[];

for i=1:length(SELECTED_PEAKS)
	data=load(FILENAMES{i});

	if ~isfield(data,'ephys_data')
		ephys_data=data.intan_data;
	else
		ephys_data=data.ephys_data;
    	end
    
    	channels=data.channels;

	for j=1:length(channels)

		% loop and if any channels are not included in the channel_label vector, include!

		if ~any(channels(j)==channel_labels)
			channel_labels=[channel_labels channels(j)];
		end
	end

end

channel_labels=sort(channel_labels);
disp(['Found channels ' num2str(channel_labels)]);

parfor i=1:length(SELECTED_PEAKS)

	if length(SELECTED_PEAKS{i})<1
		continue;
	end

	data=load(FILENAMES{i});
	sound_data=single(data.mic_data);
	
	if ~isfield(data,'ephys_data')
		ephys_data=data.intan_data;
	else
		ephys_data=data.ephys_data;
    	end
        	
	fs=data.fs;
	for j=1:length(SELECTED_PEAKS{i})
		
		peakLoc=SELECTED_PEAKS{i}(j);

		% the startpoint needs to be adjusted using the following formula
		% peaklocation*(WINDOWLENGTH-OVERLAP)*SUBSAMPLING-WINDOWLENGTH

		startpoint=(peakLoc*(N-OVERLAP)*DOWNSAMPLING)-N;
		endpoint=startpoint+TEMPLATESIZE;

		if startpoint/fs>=TIME_RANGE(1) & endpoint/fs<=TIME_RANGE(2)

			if length(sound_data)>endpoint && startpoint>0
				
				counter=counter+1;
			
			end

		end
	end
end

disp(['Found ' num2str(counter) ' trials ']);

%%%%

EPHYS_DATA=zeros(TEMPLATESIZE+1,counter,length(channel_labels),'single');
MIC_DATA=zeros(TEMPLATESIZE+1,counter,'single');

disp('Extracting data');
fprintf(1,['Progress:  ' blanks(nblanks)]);

trial=1;
for i=1:length(SELECTED_PEAKS)

	fprintf(1,formatstring,round((i/length(SELECTED_PEAKS))*100));

	if length(SELECTED_PEAKS{i})<1
		continue;
	end

	data=load(FILENAMES{i});
	sound_data=single(data.mic_data);
	
	if ~isfield(data,'ephys_data')
		ephys_data=data.intan_data;
	else
		ephys_data=data.ephys_data;
    	end
	 
	fs=data.fs;
	channels=data.channels;

	for j=1:length(SELECTED_PEAKS{i})
		
		peakLoc=SELECTED_PEAKS{i}(j);

		% the startpoint needs to be adjusted using the following formula
		% peaklocation*(WINDOWLENGTH-OVERLAP)*SUBSAMPLING-WINDOWLENGTH

		startpoint=(peakLoc*(N-OVERLAP)*DOWNSAMPLING)-N;
		endpoint=startpoint+TEMPLATESIZE;

		if startpoint/fs>=TIME_RANGE(1) & endpoint/fs<=TIME_RANGE(2)

			if length(sound_data)>endpoint && startpoint>0

				USED_FILENAMES{end+1}=FILENAMES{i};
                		MIC_DATA(:,trial)=single(sound_data(startpoint:endpoint));               
                
				% if we have differences in channel number, how to resolve?

				for k=1:length(channels)
                    
					ch_idx=find(channels(k)==channel_labels);
					EPHYS_DATA(:,trial,ch_idx)=single(ephys_data(startpoint:endpoint,k));

				end

				trial=trial+1;

			end

		end
	end
end

CHANNELS=channel_labels;

fprintf('\n');

end

