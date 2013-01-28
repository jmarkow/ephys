function ephys_pipeline_sua_track_standalone(CANDIDATEFILE,CANDIDATECLUST,CONFIG)
%
%
%
%
%
%

% read in a template FILE from BIRDID/spike_data/CELLID
% need a TEMPLATECONFIG with the unit cluster (other options for the future)
% need a CANDIDATEFILE with the waveform and ISI to compare
% CLUSTER for the CANDIDATE must be specified (parse with sed)
% need global config for location of boundary points

if isdeployed
	tmp=CANDIDATECLUST;
	clear CANDIDATECLUST;
	CANDIDATECLUST=str2num(tmp);
	clear tmp;
end

global_parameters=ephys_pipeline_readconfig(CONFIG);

[path,file,ext]=fileparts(CANDIDATEFILE);
tokens=regexp(file,'[0-9]+','match');
candidate_channel=str2num(tokens{1});

tokens=regexp(path,filesep,'split');

% check for other templates this could match, if no template exists
% create a new directory

% unit file simply specifies what cluster to look in, perhaps later add option to hook
% into postproc for LFP and coherence processing

% TODO:  deal with non-default directory structures (or, simply enforce the standard ones more rigorously)

suatype=tokens{end};
extractionstring=tokens{end-4};
datestring=tokens{end-6};
recstring=tokens{end-7};
birdidstring=tokens{end-8};

% cellid is the directory that contains the templatefile

bookkeeping_dir=fullfile(global_parameters.bookkeeping_dir);
savedir=fullfile(bookkeeping_dir,birdidstring,recstring,['ch' num2str(candidate_channel) ]);

% check for any possible templates in the channel
%%%%%

savedir_list=dir(savedir);
template_list={};

for i=1:length(savedir_list)
	if savedir_list(i).isdir && savedir_list(i).name(1)~='.'
		template_list{end+1}=fullfile(savedir,savedir_list(i).name);
	end
end

template_status=zeros(1,length(template_list),'int8');

candidate=load(CANDIDATEFILE,'clusterwindows','clusterisi');

candidate_wave=zscore(mean(candidate.clusterwindows{CANDIDATECLUST},2));
candidate_isi=candidate.clusterisi{CANDIDATECLUST};

% load the training data, specify the cluster boundary

trainfile=fullfile(global_parameters.bookkeeping_dir,'train_data','training_data.mat');

if exist(trainfile,'file')
	load(trainfile,'train_matrix','class');
else
	warning('ephysPipeline:suatracking:notrainingdata','Could not find training data %s',trainfile);
	return;
end

% thumbs up or down?

traindata(:,1)=atanh(train_matrix(:,1));
traindata(:,2)=train_matrix(:,2);

class=class+1;

[C,err,P,logp,coeff]=classify(traindata(class==2,:),traindata,class,'quadratic');

threshold=quantile(P(:,1),.73); % with the current training set this produces 5% error rate
n=(2*threshold)/(1+2*threshold);
p=1/(1+2*threshold);
priors=[ n p ];
%priors=[ .5 .5 ];

filedir=path;

for i=1:length(template_list)

	% in each directory use the earliest extraction as the template

	extract_list={};
	cellid_list={};
	template_extraction=dir(template_list{i});

	for j=1:length(template_extraction)
		if template_extraction(j).isdir && template_extraction(j).name(1)~='.'
			extract_list{end+1}=fullfile(template_list{i},template_extraction(j).name);
			cellid_list{end+1}=template_extraction(j).name;
		end
	end

	% check for any potential matches

	for j=1:length(extract_list)

		fid=fopen(fullfile(extract_list{j},'cellinfo.txt'),'r');

		readdata=textscan(fid,'%s%[^\n]','commentstyle','#',...
			'delimiter','\t','MultipleDelimsAsOne',1);	

		fclose(fid);

		% read in cluster and channel number from cellinfo.txt

		cluster=str2num(readdata{2}{find(strcmpi(readdata{1},'cluster:'))});
		channel=str2num(readdata{2}{find(strcmpi(readdata{1},'channel:'))});
		
		% sua file should match, otherwise skip

		clustfile=fullfile(extract_list{j},'raster',['sua_channels ' num2str(channel) '.mat']);

		if exist(clustfile,'file')	
			template=load([clustfile],'clusterwindows','clusterisi');
		else
			warning('ephysPipeline:singleunittrack:nosuafile','Could not find sua file at %s',clustfile);
			continue;
		end

		template_wave=zscore(mean(template.clusterwindows{cluster},2));
		template_isi=template.clusterisi{cluster};

		if length(template_wave)~=length(candidate_wave)

			% if the spike waveforms do not have equal length bail

			disp('Waveforms not of equal length, cannot process...');
			print_done_signal(filedir,CANDIDATEFILE);
			return;
		end

		[r,lags]=xcov(template_wave,candidate_wave,'coeff');
		wavescore=atanh(max(r));
		isiscore=sqrt(kld(template_isi,candidate_isi,'jsd',1));

		[C]=classify([wavescore isiscore],traindata,class,'quadratic',priors);

		% break after we read the template, presumably the first file

		break; 

	end

	template_status(i)=C;

end

namebase=global_parameters.unit_name;
template_status=template_status>1;

if all(template_status==0) || isempty(template_status)
	dirnum=length(template_list)+1;
	cellid=[ namebase '_' num2str(dirnum) ];
	fprintf('New cell found, cluster number %i\n\n',dirnum);
elseif sum(template_status)==1
	dirnum=find(template_status);
	[path,file,ext]=fileparts(template_list{dirnum});
	cellid=file;
	fprintf('Match found in %s\n\n',template_list{i});
elseif sum(template_status)>1
	disp('More than one match, using first match!');
	dirnum=find(template_status)
	dirnum=dirnum(1);
	[path,file,ext]=fileparts(template_list{dirnum});
	cellid=file;
	fprintf('Match found in %s\n\n',template_list{i});
end

savedir=fullfile(savedir,cellid,[ datestring ' (' extractionstring ',' lower(suatype) ')']);

% copy any single unit rasters for safe keeping...

% 2 is accept, 1 reject

if ~exist(fullfile(savedir,'raster'),'dir')
	mkdir(fullfile(savedir,'raster'));
end

fprintf('Extraction:\t%s\nDate:\t%s\nRec:\t%s\nBird:\t%s\n\n',...
	extractionstring,datestring,recstring,birdidstring);

fid=fopen(fullfile(savedir,'cellinfo.txt'),'w');

fprintf(fid,'Path:\t\t%s\nRec ID:\t\t%s\nBird ID:\t%s\nDate:\t\t%s\nCell ID:\t%s\nExtraction:\t%s\nChannel:\t%g\nCluster:\t%g',...
	filedir,recstring,birdidstring,datestring,cellid,extractionstring,candidate_channel,CANDIDATECLUST);

fclose(fid);

try
	copyfile(fullfile(filedir,['ephys_sua_freqrange*electrode_' num2str(candidate_channel) ...
		'_raster_cluster_' num2str(CANDIDATECLUST) '.*' ]),fullfile(savedir,'raster'));
	copyfile(fullfile(filedir,['ephys_sua_freqrange*electrode_' num2str(candidate_channel) ...
		'_stats_cluster_*' ]),fullfile(savedir,'raster'));
	copyfile(fullfile(filedir,['sua_channels ' num2str(candidate_channel) '.mat']),fullfile(savedir,'raster'));
catch
	warning('ephysPipeline:singleunitpostproc:nocopy','Could not copy single unit raster file');
end

% print out filename with . in front to signal that we've processed it

print_done_signal(filedir,CANDIDATEFILE);

end

function print_done_signal(FILEDIR,CANDIDATEFILE)

[path,file,ext]=fileparts(CANDIDATEFILE);

fid=fopen(fullfile(FILEDIR,[ '.' file ext]),'w');
fclose(fid);

end

% multiple numbers in config file are read in as strings
% if we have designated to skip coherence, that's it
