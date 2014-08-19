function ephys_pipeline_cluster_standalone(TEMPLATEFILE,DATAFILE,CLASSIFYFILE,CONFIG)
%
%
%
%


parameters=ephys_pipeline_readconfig(CONFIG);

% TODO update to allow padding

colors=hot(63); % uint8 colormap

disp_minfs=1;
disp_maxfs=10e3;
parameters.smscore_n=1024;
padding=[];

% takes the classifier object, template features and data features and extracts any instances
% of the template as determined by the classifier 

is_legacy=check_legacy(TEMPLATEFILE);

if is_legacy
	load(TEMPLATEFILE,'template_features','TEMPLATE','padding'); % we'll get warnings if padding DNE
	template.features=template_features;
	template.data=TEMPLATE;
	clearvars template_features TEMPLATE;
else
	load(TEMPLATEFILE,'template','padding');
end

% if for some reason we cannot load datafile (e.g. if the pipeline crashes, then delete it so it's 
% recomputed on the next run)

try
	load(DATAFILE,'features','features_parameters');
catch
	warning('ephysPipeline:clusterstandalone:loadingattempt1','Could not load %s, pausing 30 seconds',DATAFILE);
	pause(30);
	try
		load(DATAFILE,'features','features_parameters');
	catch
		warning('ephysPipeline:clusterstandalone:loadingattempt2','Could not load %s, deleting',DATAFILE);
		delete(DATAFILE);
		return;
	end
end

if exist('features_parameters','var') & isfield(template.feature_parameters)
	chk1=template.feature_parameters.low_cutoff==features_parameters.low_cutoff;
	chk2=template.feature_parameters.high_cutoff==features_parameters.high_cutoff;
	
	if (~chk1)|(~chk2)
		error('Mismatch parameters between template (lo %g hi %g) and data (lo %g hi %g)',...
			template.feature_parameters.low_cutoff,template.feature_parameters.high_cutoff,...
			features_parameters.low_cutoff,features_parameters.high_cutoff);
	end
else
	chk=size(features{1},1)==size(template.features{1},1);
	if ~chk
		error('Feature size mismatch between template (%g) and data (%g)',...
			size(features{1},1),size(template.features{1},1));
	end
end

load(CLASSIFYFILE,'cluster_choice','classobject');

if length(padding)==2
	disp(['Will pad extractions using the following pads:  ' num2str(padding)]);
end

[junk,templength]=size(template.features{1});
templength=templength-1;
fulltemplength=length(template.data);

% before extracting make sure data file exists

[junk,featureslength]=size(features{1});

score_temp={};
temp_mat=[];

% we'll write the done signal (.dotfile) to the datafile path

[templatepath,junk1,junk2]=fileparts(TEMPLATEFILE);

tokens=regexp(templatepath,filesep,'split');

% the template name is the directory the template data is stored in, this will
% be a prefix for the done signal

templatename=tokens{end};

[datapath,file,ext]=fileparts(DATAFILE);

rawfile=fullfile(datapath,'..',[ file(1:end-6) '.mat']);


if ~exist(rawfile,'file')
	return;
end

savedir=fullfile(datapath,'..',templatename);

imagedir=fullfile(savedir,'gif');
wavdir=fullfile(savedir,'wav');
matdir=fullfile(savedir,'mat');

if ~exist(imagedir,'dir')
	mkdir(imagedir);
end

if ~exist(wavdir,'dir');
	mkdir(wavdir);
end

if ~exist(matdir,'dir');
	mkdir(matdir);
end

fid=fopen(fullfile(savedir,'.songextraction'),'w');
fclose(fid);

% compute in a sliding 1 sec window perhaps...

disp(['Checking ' ' for matches to ' TEMPLATEFILE]);

error_count=0;

for i=1:length(features)
	
	score_temp{i}=[];
	energy=[];

	for j=1:featureslength-templength
		score_temp{i}=[score_temp{i} sum(sum(abs(features{i}(:,j:j+templength)-template.features{i})))];
		energy=[energy mean(mean(features{4}(:,j:j+templength)))];
	end

	% taking only points where we have significant energy, otherwise long bouts of silence
	% will screw up our normalization

	% assume that silence is the min point
	
	silent=min(energy);
	
	% where is energy greater than the min + 1 STD (very conservative)

	idxs=energy>(silent+.5*std(energy));

	%find(idxs)
	
	score_temp{i}=score_temp{i}-mean(score_temp{i}(idxs));
	
	%score_temp{i}=score_temp{i}-median(score_temp{i});

	% divide IQR by 1<x<2 to scale how conservative the clustering is, increasing
	% x makes the clustering LESS conservative, 1.349-->SD

	%scale=iqr(score_temp{i}(idxs));
	scale=std(score_temp{i}(idxs));
	
	%scale=iqr(score_temp{i}):
	score_temp{i}=score_temp{i}/scale;
	score_temp{i}(score_temp{i}>0)=0;
	score_temp{i}=abs(score_temp{i});


end

attributes=length(score_temp);
product_score=score_temp{1};

for i=2:attributes, product_score=product_score.*score_temp{i}; end

%figure();plot(product_score)
% need to tag the file as processed at the end of this fiasco

if length(product_score)<3
	disp('too short');
	write_done_signal(savedir,[file ext]);
	return;
end

[pks,locs]=findpeaks(product_score,'MINPEAKHEIGHT',.005);

if isempty(locs)
	disp('no locs');
	write_done_signal(savedir,[file ext]);
	return;
end

curvature=gradient(gradient(product_score));

% for each peak create matrix with rows for observations and columns features
% this will get passed to the classifier object, thumbs up/down, then extract if necessary

for j=1:attributes, feature_mat(:,j)=log(score_temp{j}(locs)); end

feature_mat(:,attributes+1)=log(product_score(locs));
feature_mat(:,attributes+2)=log(abs(curvature(locs)));

peak_locs=locs;

% classify

labels=svmclassify(classobject,feature_mat);

% instances where we have a hit

hits=find(labels==cluster_choice);

if isempty(hits)
	disp('no hits');
	write_done_signal(savedir,[file ext]);
	return;
end

% else, extract

templength=templength+parameters.smscore_n;


is_legacy=check_legacy(rawfile);

if is_legacy
	load(rawfile,'audio_data','ephys_data','channels','fs','ttl_data','t','start_datenum');
	
	ephys.data=ephys_data;
	ephys.fs=fs;
	ephys.t=t;
	ephys.ports='';

	audio.audio.data;
	audio.fs=fs;
	audio.t=t;
	
	ttl.data=ttl_data;
	file_datenum=start_datenum;

else
	load(rawfile,'audio','ephys','ttl','file_datenum');
end

[b,a]=ellip(5,.2,80,[700]/(audio.fs/2),'high');


if audio.fs~=ephys.fs
	error('Unequal audio (%g) and ephys (%g) sampling rates not yet supported',audio.fs,ephys.fs);
end

% write images to 'gif', wav to 'wav' and 'mat' to mat again

[sonogram_im sonogram_f sonogram_t]=pretty_sonogram(filtfilt(b,a,double(audio.data)),audio.fs,'n',500,'overlap',450,'low',2.5);
startidx=max([find(sonogram_f<=disp_minfs)]);
stopidx=min([find(sonogram_f>=disp_maxfs)]);
sonogram_im=sonogram_im(startidx:stopidx,:);
sonogram_im=flipdim(sonogram_im,1);

[f,t]=size(sonogram_im);
im_son_to_vec=(length(audio.data)-450)/t;

% make sure we don't pull out overlapping data

startpoints=hits.*(parameters.smscore_n-parameters.smscore_overlap)*parameters.smscore_downsampling;

% prevent overlap by only using hits sufficiently spaced from each other

prev_endpoint=0;

store_audio=audio;
store_ephys=ephys;
store_ttl=ttl;

disp(['Extracting ' num2str(length(hits)) ' hits']);

for i=1:length(hits)

	hitloc=peak_locs(hits(i));

	startpoint=(hitloc*(parameters.smscore_n-parameters.smscore_overlap)*parameters.smscore_downsampling);
	endpoint=startpoint+fulltemplength;

	if length(padding)==2
		startpoint=startpoint-floor(padding(1)*audio.fs);
		endpoint=endpoint+ceil(padding(2)*audio.fs);
	end

	if startpoint-prev_endpoint<=0
		disp('Hits are overlapping, skipping...');
		continue;
	end

	if length(store_audio.data)>endpoint && startpoint>0	

		audio.data=store_audio.data(startpoint:endpoint);
		ephys.data=store_ephys.data(startpoint:endpoint,:);
		audio.t=store_audio.t(startpoint:endpoint);
		ephys.t=store_ephys.t(startpoint:endpoint);

		%audio.fs=fs;
		%ephys.fs=fs;

		if isfield(store_ttl,'data') & ~isempty(store_ttl.data)
			ttl.data=store_ttl.data(startpoint:endpoint);
			ttl.t=store_ttl.t(startpoint:endpoint);
		else
			ttl.data=[];
			ttl.t=[];
		end

		savename=[ file(1:end-6) '_' templatename '_' num2str(i)];

		save(fullfile(matdir,[savename '.mat']),'audio','ephys','ttl','file_datenum');

		% write out the extraction

		sonogram_im(1:10,ceil(startpoint/im_son_to_vec):ceil(endpoint/im_son_to_vec))=63;

		[chunk_sonogram_im chunk_sonogram_f chunk_sonogram_t]=pretty_sonogram(filtfilt(b,a,double(audio.data)),audio.fs,'low',2.5);

		startidx=max([find(chunk_sonogram_f<=disp_minfs)]);
		stopidx=min([find(chunk_sonogram_f>=disp_maxfs)]);
		chunk_sonogram_im=chunk_sonogram_im(startidx:stopidx,:);
		chunk_sonogram_im=flipdim(chunk_sonogram_im,1);
		imwrite(uint8(chunk_sonogram_im),colors,fullfile(imagedir,[savename '.gif']),'gif');
		wavwrite(audio.data,audio.fs,fullfile(wavdir,[savename '.wav']));
		prev_endpoint=endpoint;
	
	end

	

end

imwrite(uint8(sonogram_im),colors,fullfile(imagedir,[ file(1:end-6) '.gif']),'gif');

write_done_signal(savedir,[file ext]);

end

% need this so the bash script does not process the same file for the same template again

function write_done_signal(savedir,dataname)

fid=fopen(fullfile(savedir,[ '.'  dataname ]),'w');
fclose(fid);

end

