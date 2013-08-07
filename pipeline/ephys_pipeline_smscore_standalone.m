function smscore_standalone(FILENAME,CONFIG)

[path,file,ext]=fileparts(FILENAME);

if ext=='.mat'
	load(FILENAME,'mic_data','fs');
elseif ext=='.wav'
	[mic_data,fs]=wavread(FILENAME);
else 
	error('Cannot recognize file type!');
end

TTL=0;
lowfs=[];
highfs=[];

% read CONFIG file for N OVERLAP AND SIGMA

parameters=ephys_pipeline_readconfig(CONFIG);

% are processing feedback data?

if length(file)>2 & strcmp(file(end-2:end),'ttl')
	TTL=1;
	lowfs=parameters.smscore_ttl_lowfs;
	highfs=parameters.smscore_ttl_highfs;
end

fprintf('Parameters\n\n%-10s%-10s%-10s%-10s\n\n','N','Overlap','FiltScale','Downsampling');
fprintf('%-10d%-10d%-10d%-10d\n\n\n',...
	parameters.smscore_n,parameters.smscore_overlap,parameters.smscore_filter_scale,parameters.smscore_downsampling);

features=ephys_pipeline_smscore(mic_data,fs,'n',parameters.smscore_n,'overlap',parameters.smscore_overlap,...
	'filter_scale',parameters.smscore_filter_scale,'downsampling',parameters.smscore_downsampling,...
	'norm_amp',parameters.smscore_norm_amp,'lowfs',lowfs,'highfs',highfs);


[path,file,ext]=fileparts(FILENAME);

if ~exist(fullfile(path,'syllable_data'),'dir')
	mkdir(fullfile(path,'syllable_data'));
end

save(fullfile(path,'syllable_data',[ file '_score.mat']),'features','TTL','lowfs','highfs');
