function smscore_standalone(FILENAME,CONFIG)

[path,file,ext]=fileparts(FILENAME);

if ext=='.mat'
	
	% check for filetype

	is_legacy=check_legacy(FILENAME);
	
	if is_legacy
		load(FILENAME,'mic_data','fs');
		audio.data=mic_data;
		audio.fs=fs;
		clearvars fs mic_data;
	else
		load(FILENAME,'audio');
	end

elseif ext=='.wav'
	[audio.data,audio.fs]=wavread(FILENAME);
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

	% if we match the bird name, use this bird's specific frequency parameters
	%

	TTL=1;
	tmp=fieldnames(parameters);

	[~,file,~]=fileparts(FILENAME);
	tokens=regexp(file,parameters.delimiter,'split');
	birdid=tokens{2};

	matches=find(~cellfun(@isempty,strfind(tmp,birdid)));
	
	if length(matches)>2
		error('Too many matches in configuration file %g',length(matches));
	end

	if length(matches)<2
		disp('Did not detect bird specific parameters');
		lowfs=parameters.smscore_ttl_lowfs;
		highfs=parameters.smscore_ttl_highfs;
	else
		for i=1:length(matches)
			fieldname=tmp{matches(i)};

			if strcmp(fieldname,[ birdid '_smscore_ttl_lowfs' ])
				lowfs=parameters.(fieldname);
			elseif strcmp(fieldname,[ birdid '_smscore_ttl_highfs' ])
				highfs=parameters.(fieldname);
			else
				error('Did not recognize fieldname %s',fieldname);
			end
		end
	end
	fprintf('Set bird specific parameters for %s:\nLow\t%g\nHigh\t%g\n',birdid,lowfs,highfs);

end

fprintf('Parameters\n\n%-10s%-10s%-10s%-10s\n\n','N','Overlap','FiltScale','Downsampling');
fprintf('%-10d%-10d%-10d%-10d\n\n\n',...
	parameters.smscore_n,parameters.smscore_overlap,parameters.smscore_filter_scale,parameters.smscore_downsampling);

[features,features_parameters]=ephys_pipeline_smscore(audio.data,audio.fs,'n',parameters.smscore_n,'overlap',parameters.smscore_overlap,...
	'filter_scale',parameters.smscore_filter_scale,'downsampling',parameters.smscore_downsampling,...
	'norm_amp',parameters.smscore_norm_amp,'lowfs',lowfs,'highfs',highfs);


[path,file,ext]=fileparts(FILENAME);

if ~exist(fullfile(path,'syllable_data'),'dir')
	mkdir(fullfile(path,'syllable_data'));
end

save(fullfile(path,'syllable_data',[ file '_score.mat']),'features','features_parameters','TTL','parameters');
