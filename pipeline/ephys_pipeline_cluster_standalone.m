function cluster_standalone(TEMPLATEFILE,DATAFILE,CLASSIFYFILE)
%
%
%
%

colors=hot(63); % uint8 colormap

disp_minfs=1e3;
disp_maxfs=7.5e3;
NW=1024;
SR=25e3;

% takes the classifier object, template features and data features and extracts any instances
% of the template as determined by the classifier 

load(TEMPLATEFILE,'template_features','TEMPLATE');
load(DATAFILE,'features');
load(CLASSIFYFILE,'cluster_choice','classobject');

% ANNOYING we have to train the classifier each time in the compiled version since objects are not supported...

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

for i=1:length(features)
	score_temp{i}=[];

	for j=1:featureslength-templength
		score_temp{i}=[score_temp{i} sum(sum(abs(features{i}(:,j:j+templength)-template_features{i})))];
	end

	score_temp{i}=score_temp{i}-mean(score_temp{i});
	score_temp{i}=score_temp{i}/std(score_temp{i});
	score_temp{i}(score_temp{i}>0)=0;
	score_temp{i}=abs(score_temp{i});

end

attributes=length(score_temp);
product_score=score_temp{1};

for i=2:attributes, product_score=product_score.*score_temp{i}; end

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

templength=templength+NW;
data=load(rawfile,'mic_data','ephys_data','channels','fs');

% write images to 'gif', wav to 'wav' and 'mat' to mat again

[sonogram_im sonogram_f sonogram_t]=pretty_sonogram(data.mic_data,data.fs,'n',500,'overlap',450,'low',3.5);
startidx=max([find(sonogram_f<=disp_minfs)]);
stopidx=min([find(sonogram_f>=disp_maxfs)]);
sonogram_im=sonogram_im(startidx:stopidx,:);
sonogram_im=flipdim(sonogram_im,1);

[f,t]=size(sonogram_im);
im_son_to_vec=(length(data.mic_data)-450)/t;

for i=1:length(hits)

	hitloc=peak_locs(hits(i));

	startpoint=(hitloc*(NW-1e3)*5);
	endpoint=startpoint+fulltemplength;

	if length(data.mic_data)>endpoint && startpoint>0	

		mic_data=data.mic_data(startpoint:endpoint);
		ephys_data=data.ephys_data(startpoint:endpoint,:);
		channels=data.channels;
		fs=data.fs;

		savename=[ file(1:end-6) '_' templatename '_' num2str(i)];

		save(fullfile(matdir,[savename '.mat']),'mic_data','ephys_data','channels','fs');

		% write out the extraction

		sonogram_im(1:10,ceil(startpoint/im_son_to_vec):ceil(endpoint/im_son_to_vec))=63;

		[chunk_sonogram_im chunk_sonogram_f chunk_sonogram_t]=pretty_sonogram(mic_data,fs,'low',3.5);

		startidx=max([find(chunk_sonogram_f<=disp_minfs)]);
		stopidx=min([find(chunk_sonogram_f>=disp_maxfs)]);
		chunk_sonogram_im=chunk_sonogram_im(startidx:stopidx,:);
		chunk_sonogram_im=flipdim(chunk_sonogram_im,1);
		imwrite(uint8(chunk_sonogram_im),colors,fullfile(imagedir,[savename '.gif']),'gif');
		wavwrite(mic_data,fs,fullfile(wavdir,[savename '.wav']));
	
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

