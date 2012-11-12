function ephys_pipeline_sua_track_standalone(TEMPLATEFILE,TEMPLATECONFIG,CANDIDATEFILE,CANDIDATECLUST,CONFIG)
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

disp([num2str(nargin)]);

if isdeployed
	tmp=CANDIDATECLUST;
	clear CANDIDATECLUST;
	CANDIDATECLUST=str2num(tmp);
	clear tmp;
end

disp([TEMPLATEFILE])
disp([TEMPLATECONFIG])
disp([CANDIDATEFILE])
disp([num2str(CANDIDATECLUST)])
disp([CONFIG])

global_parameters=ephys_pipeline_readconfig(CONFIG);
unit_parameters=ephys_pipeline_readconfig(TEMPLATECONFIG);

[path,file,ext]=fileparts(TEMPLATEFILE);

tokens=regexp(file,'[0-9]+','match');
channel=str2num(tokens{1});

tokens=regexp(path,filesep,'split');

% path containing the unit file is the CELLID

cellid=tokens{end};

[path,file,ext]=fileparts(CANDIDATEFILE);
tokens=regexp(path,filesep,'split');

% unit file simply specifies what cluster to look in, perhaps later add option to hook
% into postproc for LFP and coherence processing

% TODO:  deal with non-default directory structures (or, simply enforce the standard ones more rigorously)

suatype=tokens{end}
extractionstring=tokens{end-4}
datestring=tokens{end-6}
recstring=tokens{end-7}
birdidstring=tokens{end-8}

% cellid is the directory that contains the templatefile

bookkeeping_dir=fullfile(global_parameters.bookkeeping_dir);
savedir=fullfile(bookkeeping_dir,birdidstring,cellid,...
	[ datestring ' (' extractionstring ',' lower(suatype) ')']);

disp(['Post-processing ' cellid]);
disp(['Will save in directory ' savedir]);

filedir=path;

% copy any single unit rasters for safe keeping...

template=load(TEMPLATEFILE,'clusterwindows','clusterisi');
candidate=load(CANDIDATEFILE,'clusterwindows','clusterisi');

template_wave=zscore(mean(template.clusterwindows{unit_parameters.cluster},2));
template_isi=template.clusterisi{unit_parameters.cluster};

% candidate cluster is read in by filename tagging (meets all simple criteria)

candidate_wave=zscore(mean(candidate.clusterwindows{CANDIDATECLUST},2));
candidate_isi=candidate.clusterisi{CANDIDATECLUST};

% check peak correlation of scaled waveforms

[r,lags]=xcov(template_wave,candidate_wave,'coeff');
wavescore=atanh(max(r));
isiscore=sqrt(kld(template_isi,candidate_isi,'jsd',1));

load(fullfile(global_parameters.bookkeeping_dir,'train_data','training_data.mat'),...
	'train_matrix','class');

% thumbs up or down?

traindata(:,1)=atanh(train_matrix(:,1));
traindata(:,2)=train_matrix(:,2);

class=class+1;

%[C]=classify([wavescore isiscore],traindata,class+1)

[C,err,P,logp,coeff]=classify(traindata(class==2,:),traindata,class,'quadratic');

threshold=quantile(P(:,1),.9);
n=(2*threshold)/(1+2*threshold);
p=1/(1+2*threshold);
priors=[ n p ];
%priors=[ .5 .5 ];

[C]=classify([wavescore isiscore],traindata,class,'quadratic',priors);

% create a dotfile at the end of this
% 2 is accept, 1 reject

if C==2

	disp('Single-unit match');
	
	if ~exist(fullfile(savedir,'raster'),'dir')
		mkdir(fullfile(savedir,'raster'));
	end

	
	fprintf('Extraction:\t%s\nDate:\t%s\nRec:\t%s\nBird:\t%s\n\n',...
		extractionstring,datestring,recstring,birdidstring);

	fid=fopen(fullfile(savedir,'cellinfo.txt'),'w');

	fprintf(fid,'Bird ID:\t%s\nDate:\t\t%s\nCell ID:\t%s\nExtraction:\t%s\nChannel:\t%g\nCluster:\t%g',...
		birdidstring,datestring,cellid,extractionstring,channel,CANDIDATECLUST);

	fclose(fid);

	try
		copyfile(fullfile(filedir,['ephys_sua_freqrange*electrode_' num2str(channel) ...
			'_raster_cluster_' num2str(CANDIDATECLUST) '.*' ]),fullfile(savedir,'raster'));
		copyfile(fullfile(filedir,['ephys_sua_freqrange*electrode_' num2str(channel) ...
			'_stats_cluster_*' ]),fullfile(savedir,'raster'));
		copyfile(fullfile(filedir,['sua_channels ' num2str(channel) '.mat']),fullfile(savedir,'raster'));
	catch
		warning('ephysPipeline:singleunitpostproc:nocopy','Could not copy single unit raster file');
	end
end

% print out filename with . in front to signal that we've processed it

fid=fopen(fullfile(filedir,[ '.' cellid '_' file ext]),'w');
fclose(fid);

% multiple numbers in config file are read in as strings
% if we have designated to skip coherence, that's it
