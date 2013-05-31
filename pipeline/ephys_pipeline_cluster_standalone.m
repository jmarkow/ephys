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
parameters.fs=25e3;
padding=[];
% takes the classifier object, template features and data features and extracts any instances
% of the template as determined by the classifier 

load(TEMPLATEFILE,'template_features','TEMPLATE','padding'); % we'll get warnings if padding DNE

% if for some reason we cannot load datafile (e.g. if the pipeline crashes, then delete it so it's 
% recomputed on the next run)

try
	load(DATAFILE,'features');
catch
	warning('ephysPipeline:clusterstandalone:loadingattempt1','Could not load %s, pausing 30 seconds',DATAFILE);
	pause(30);
	try
		load(DATAFILE,'features');
	catch
		warning('ephysPipeline:clusterstandalone:loadingattempt2','Could not load %s, deleting',DATAFILE);
		delete(DATAFILE);
		return;
	end
end

load(CLASSIFYFILE,'cluster_choice','classobject');

if length(padding)==2
	disp(['Will pad extractions using the following pads:  ' num2str(padding)]);
end

[junk,templength]=size(template_features{1});
templength=templength-1;
fulltemplength=length(TEMPLATE);

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

rawfile

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

for i=1:length(features)
	
	score_temp{i}=[];
	energy=[];

	for j=1:featureslength-templength
		score_temp{i}=[score_temp{i} sum(sum(abs(features{i}(:,j:j+templength)-template_features{i})))];
		energy=[energy mean(mean(features{4}(:,j:j+templength)))];
	end

	% taking only points where we have significant energy, otherwise long bouts of silence
	% will screw up our normalization

	% assume that silence is the min point
	
	silent=min(energy)
	
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
data=load(rawfile,'mic_data','ephys_data','channels','fs');

% write images to 'gif', wav to 'wav' and 'mat' to mat again

[sonogram_im sonogram_f sonogram_t]=pretty_sonogram(data.mic_data,data.fs,'n',500,'overlap',450,'low',2.5);
startidx=max([find(sonogram_f<=disp_minfs)]);
stopidx=min([find(sonogram_f>=disp_maxfs)]);
sonogram_im=sonogram_im(startidx:stopidx,:);
sonogram_im=flipdim(sonogram_im,1);

[f,t]=size(sonogram_im);
im_son_to_vec=(length(data.mic_data)-450)/t;

% make sure we don't pull out overlapping data

startpoints=hits.*(parameters.smscore_n-parameters.smscore_overlap)*parameters.smscore_downsampling;

% prevent overlap by only using hits sufficiently spaced from each other

%keeppoints=find(diff(startpoints)>fulltemplength+sum(padding)*parameters.fs);
%hits=hits(keeppoints);

%

prevendpoint=0;

for i=1:length(hits)

	hitloc=peak_locs(hits(i));

	startpoint=(hitloc*(parameters.smscore_n-parameters.smscore_overlap)*parameters.smscore_downsampling);
	endpoint=startpoint+fulltemplength;

	if length(padding)==2
		startpoint=startpoint-floor(padding(1)*parameters.fs);
		endpoint=endpoint+ceil(padding(2)*parameters.fs);
	end

	if startpoint-prevendpoint<=0
		disp('Hits are overlapping, skipping...');
		continue;
	end

	if length(data.mic_data)>endpoint && startpoint>0	

		mic_data=data.mic_data(startpoint:endpoint);
		ephys_data=data.ephys_data(startpoint:endpoint,:);
		channels=data.channels;
		fs=data.fs;

		savename=[ file(1:end-6) '_' templatename '_' num2str(i)];

		save(fullfile(matdir,[savename '.mat']),'mic_data','ephys_data','channels','fs');

		% write out the extraction

		sonogram_im(1:10,ceil(startpoint/im_son_to_vec):ceil(endpoint/im_son_to_vec))=63;

		[chunk_sonogram_im chunk_sonogram_f chunk_sonogram_t]=pretty_sonogram(mic_data,fs,'low',2.5);

		startidx=max([find(chunk_sonogram_f<=disp_minfs)]);
		stopidx=min([find(chunk_sonogram_f>=disp_maxfs)]);
		chunk_sonogram_im=chunk_sonogram_im(startidx:stopidx,:);
		chunk_sonogram_im=flipdim(chunk_sonogram_im,1);
		imwrite(uint8(chunk_sonogram_im),colors,fullfile(imagedir,[savename '.gif']),'gif');
		wavwrite(mic_data,fs,fullfile(wavdir,[savename '.wav']));
		prevendpoint=endpoint;
	
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

