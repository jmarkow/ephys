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
colors='hot';
min_f=1;
max_f=10e3;
time_range=[0 inf];
subset='';
padding=[]; % padding that will be saved with the template, in seconds (relevant for the pipeline only)
   	    % two elements vector, both specify in seconds how much time before and after to extract
	    % e.g. [.2 .2] will extract 200 msec before and after the extraction point when clustering
	    % sounds through the pipeline
lowfs=[];
highfs=[];

% smscore parameters, THESE MUST MATCH THE PIPELINE PARAMETERS IN EPHYS_PIPELINE.CFG, OTHERWISE
% THE FEATURE COMPUTATION BETWEEN THE TEMPLATE AND CANDIDATE SOUNDS WILL NOT
% BE APPROPRIATELY MATCHED

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SMSCORE PIPELINE PARAMETERS %%%%%%%%

n=1024;
overlap=1000;
filter_scale=10;
downsampling=5;
train_classifier=1;

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
		case 'train_classifier'
			train_classifier=varargin{i+1};
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

	template=select_template(fullfile(DIR));

	% compute the features of the template

	disp('Computing the spectral features of the template');
	[template.features template.feature_parameters]=ephys_pipeline_smscore(template.data,template.fs,...
		'n',n,'overlap',overlap,'filter_scale',filter_scale,'downsampling',downsampling,'lowfs',lowfs,'highfs',highfs);
	save(fullfile(proc_dir,'template_data.mat'),'template','padding');

else
	disp('Loading stored template...');
	load(fullfile(proc_dir,'template_data.mat'),'template');
	
	disp('Computing the spectral features of the template');
	[template.features template.feature_parameters]=ephys_pipeline_smscore(template.data,template.fs,...
		'n',n,'overlap',overlap,'filter_scale',filter_scale,'downsampling',downsampling,'lowfs',lowfs,'highfs',highfs);
	save(fullfile(proc_dir,'template_data.mat'),'template','padding');

end

% generate a nice sonogram of the selected template

template_fig=figure('Visible','off');
[template_image,f,t]=pretty_sonogram(template.data,template.fs,'N',1024,'overlap',1000,'low',1);

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

[junk,templength]=size(template.features{1});
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

	sound_file_features(DIR,files_to_proc,n,overlap,filter_scale,downsampling,lowfs,highfs);

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

	template_match(template.features,files_to_proc,fullfile(proc_dir,'cluster_data.mat'),templength);

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

load(fullfile(proc_dir,'cluster_results.mat'),'sorted_syllable','syllable_data','cluster');
load(fullfile(proc_dir,'cluster_data.mat'),'filenames');

if train_classifier

	disp('Training classifier on your selection...');

	% fix for MATLAB 2010a complaining about too many iterations...enforce that method=smo
	% switched to quadratic kernel function 5/28/13, linear was found to be insufficient in edge-cases

	cluster_choice=cluster.choice;

	classobject=svmtrain(syllable_data(:,[1:6]),cluster.labels,'method','smo','kernel_function','quadratic');
	save(fullfile(proc_dir,'classify_data.mat'),'classobject','cluster_choice');

end

act_templatesize=length(template.data);

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
	[agg_audio agg_ephys agg_ttl used_filenames]=extract_hits(sorted_syllable,filenames,...
		act_templatesize,spect_thresh,time_range,n,overlap,downsampling,padding);

	disp(['Saving data to ' fullfile(proc_dir,'extracted_data.mat')]);

	save(fullfile(proc_dir,'extracted_data.mat'),'used_filenames','agg_audio','agg_ephys','agg_ttl','used_filenames','-v7.3');
end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TEMPLATE SELECT %%%%%%%%%%%%%%%%%%%%

function [TEMPLATE]=select_template(DIR)

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
		load(fullfile(pathname,filename),'audio');
	end

	TEMPLATE.data=spectro_navigate(audio.data);
	TEMPLATE.fs=audio.fs;

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


function sound_file_features(DIR,SOUND_FILES,N,OVERLAP,FILTER_SCALE,DOWNSAMPLING,LOWFS,HIGHFS)

par_save = @(FILE,features,parameters,TTL) save([FILE],'features','parameters','TTL');

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

	is_legacy=check_legacy(input_file);
	
	if is_legacy
		proc=load(input_file,'mic_data','fs');
		proc.audio.data=mic_data;
		proc.audio.fs=fs;
	else
		proc=load(input_file,'audio');
	end
	
	if ~isfield(proc,'audio')
		warning('ephysPipeline:ephysCluster:errorsoundfile','Problem encountered with %s',input_file);
		continue;
	end
	
	if length(proc.audio.data)<N
		warning('ephysPipeline:ephysCluster:shortsound','Sound extraction too short in %s, skipping...',input_file);
		continue;
	end

	[sound_features,parameters]=ephys_pipeline_smscore(proc.audio.data,proc.audio.fs,...
		'n',N,'overlap',OVERLAP,'filter_scale',FILTER_SCALE,'downsampling',DOWNSAMPLING,'lowfs',LOWFS,'highfs',HIGHFS);

	% save for posterity's sake

	if ~isempty(LOWFS) & ~isempty(HIGHFS)
		TTL=1;
	else
		TTL=0;
	end

	par_save(output_file,sound_features,parameters,TTL);

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

	disp([TARGET_FILES{i}])

	for j=1:length(target)

		%template=TEMPLATE{j}-median(TEMPLATE{j});
		%template=template./mad(template);

		%targ=target{j}-median(target{j});
		%targ=target{j}./mad(target{j});

		template=TEMPLATE{j};
		targ=target{j};
		score_temp{j}=zeros(1,targetlength-TEMPLATESIZE);

		for k=1:targetlength-TEMPLATESIZE
			score_temp{j}(k)=[sum(sum(abs(targ(:,k:k+TEMPLATESIZE)-template)))];
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


function [MIC EPHYS TTL USED_FILENAMES]=extract_hits(SELECTED_PEAKS,FILENAMES,TEMPLATESIZE,SPECT_THRESH,TIME_RANGE,...
		N,OVERLAP,DOWNSAMPLING,PADDING)

TEMPLATESIZE=TEMPLATESIZE+N;
USED_FILENAMES={};

MIC=[];
TTL=[];
EPHYS=[];

disp(['Extracting cluster (' num2str(length(SELECTED_PEAKS)) ' peaks):   ']);
disp('Preallocating matrices (this may take a minute)...');

[nblanks formatstring]=progressbar(100);
counter=0;

% check for all possible channels across the whole day, a matrix will be filled with zeros if the channel
% gets knocked out somehow...

all_labels=[];
all_ports=[];

% in case we have no ephys or ttl data
%
 
ephys.labels=[];
ephys.ports=[];
ephys.fs=[];
ephys.data=[];
ttl.data=[];
ttl.fs=[];

for i=1:length(SELECTED_PEAKS)

	is_legacy=check_legacy(FILENAMES{i});

	if is_legacy		
		load(FILENAMES{i},'channels');
		labels=channels;
		ports=repmat('A',[1 length(channels)]);
	else
		load(FILENAMES{i},'ephys');

		labels=ephys.labels;	
		ports=ephys.ports;
	end	

	for j=1:length(labels)

		label_chk=labels(j)==all_labels;
		port_chk=ports(j)==all_ports;

		if isempty(label_chk), label_chk=0; end
		if isempty(port_chk), port_chk=0; end

		% loop and if any channels are not included in the channel_label vector, include!

		if ~any(label_chk&port_chk)
			all_labels=[all_labels labels(j)];
			all_ports=[all_ports ports(j)];
		end

	end

end

disp(['Found channels:  ']);

for i=1:length(all_labels)
	fprintf(1,'%i%s ',all_labels(i),all_ports(i));
end

fprintf(1,'\n');

parfor i=1:length(SELECTED_PEAKS)

	if length(SELECTED_PEAKS{i})<1
		continue;
	end

	is_legacy=check_legacy(FILENAMES{i});

	if is_legacy

		data=load(FILENAMES{i},'mic_data','fs');

		data.audio.data=mic_data;
		data.audio.fs=fs;

		data=rmfield(data,'mic_data');

	else
		data=load(FILENAMES{i},'audio');
	
	end	

	for j=1:length(SELECTED_PEAKS{i})
		
		peakLoc=SELECTED_PEAKS{i}(j);

		% the startpoint needs to be adjusted using the following formula
		% peaklocation*(WINDOWLENGTH-OVERLAP)*SUBSAMPLING-WINDOWLENGTH

		startpoint=(peakLoc*(N-OVERLAP)*DOWNSAMPLING)-N;
		endpoint=startpoint+TEMPLATESIZE;

		if startpoint/data.audio.fs>=TIME_RANGE(1) & endpoint/data.audio.fs<=TIME_RANGE(2)

			if length(data.audio.data)>endpoint && startpoint>0	
				counter=counter+1;
			end

		end
	end
end

disp(['Found ' num2str(counter) ' trials ']);

%%%%

if is_legacy
	load(FILENAMES{1},'fs','ttl_data');
	audio.fs=fs;
	ephys.fs=fs;

	if ~exist('ttl','var')
		ttl_data=[];
	end

else
	load(FILENAMES{1},'audio','ephys','ttl');
end	

EPHYS.data=zeros(TEMPLATESIZE+1,counter,length(all_labels),'single');
EPHYS.fs=ephys.fs;
MIC.data=zeros(TEMPLATESIZE+1,counter,'single');
MIC.fs=audio.fs;

if ~isempty(ttl.data)
	TTL.data=zeros(TEMPLATESIZE+1,counter,'single');
else
	TTL.data=[];
end

TTL.fs=ephys.fs;

EPHYS.labels=all_labels;
EPHYS.ports=all_ports;

disp('Extracting data');
fprintf(1,['Progress:  ' blanks(nblanks)]);

trial=1;
eflag=1;
tflag=1;

for i=1:length(SELECTED_PEAKS)

	fprintf(1,formatstring,round((i/length(SELECTED_PEAKS))*100));

	if length(SELECTED_PEAKS{i})<1
		continue;
	end

	is_legacy=check_legacy(FILENAMES{i});

	if is_legacy

		load(FILENAMES{i},'mic_data','fs','ephys_data','channels','ttl_data');
		
		audio.data=single(mic_data);
		audio.fs=fs;

		ephys.data=ephys_data;
		ephys.labels=channels;
		ephys.fs=fs;

		if ~exist('ttl_data','var')
			ttl_data=zeros(size(ephys_data));
		else
			ttl.data=ttl_data;
			ttl.fs=fs;
		end

		clearvars mic_data fs ephys_data channels;

	else
		load(FILENAMES{i},'audio','ephys','ttl');
	end	

	if ~exist('ephys','var')
		eflag=0;
	end

	if ~exist('ttl','var')
		tflag=0;
	end

	if audio.fs~=ephys.fs & eflag
		error('Audio (%g) and ephys (%g) sampling rates are not equal for file %s',...
			audio.fs,ephys.fs,FILENAMES{i});
	end
	
	for j=1:length(SELECTED_PEAKS{i})
		
		peakLoc=SELECTED_PEAKS{i}(j);

		% the startpoint needs to be adjusted using the following formula
		% peaklocation*(WINDOWLENGTH-OVERLAP)*SUBSAMPLING-WINDOWLENGTH

		startpoint=(peakLoc*(N-OVERLAP)*DOWNSAMPLING)-N;
		endpoint=startpoint+TEMPLATESIZE;

		if startpoint/audio.fs>=TIME_RANGE(1) & endpoint/audio.fs<=TIME_RANGE(2)

			if length(audio.data)>endpoint && startpoint>0

				USED_FILENAMES{end+1}=FILENAMES{i};
                		MIC.data(:,trial)=single(audio.data(startpoint:endpoint));               
              
				if ~isempty(ttl.data) & tflag
					TTL.data(:,trial)=single(ttl.data(startpoint:endpoint));
				end

				% if we have differences in channel number, how to resolve?

				if ~eflag
					continue;
				end
				
				for k=1:length(ephys.labels)
                    
					label_chk=ephys.labels(k)==all_labels;
					port_chk=ephys.ports(k)==all_ports;

					ch_idx=find(label_chk&port_chk);

					EPHYS.data(:,trial,ch_idx)=single(ephys.data(startpoint:endpoint,k));

				end

				trial=trial+1;

			end

		end
	end
end

fprintf('\n');

end

