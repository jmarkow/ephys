function ephys_songextract_aggregate(FILEDIR,SAVEDIR,CONFIG,SLEEPSTATUS)
%ephys_songextract_aggregate.m is part of the ephys pipeline
%and acts as a wrapper script for computing fields, multi-unit
%and (soon) single-unit rasters
%
%	ephys_songextract_aggregate(FILEDIR,SAVEDIR)
%
%	FILEDIR
%	directory with .mat files to process
%
%	SAVEDIR
%	base directory to store results
%

% do we want to run ephys scripts after aggregating?

% CONFIG FILE LIVES HERE!

if nargin<4
	SLEEPSTATUS=0;
end

if isdeployed
	SLEEPSTATUS=str2double(SLEEPSTATUS);
end

if nargin<3
	error('Need config file to continue!');
end

parameters=ephys_pipeline_readconfig(CONFIG);

% read config file one line at a time

listing=dir(fullfile(FILEDIR,'*.mat'));

% where do we want to store the output?

if length(listing)<1
	disp(['No files found in ' FILEDIR ' bailing...' ]);
	return;
end

% pre-allocate the matrix, store as a single

is_legacy=check_legacy(fullfile(FILEDIR,listing(1).name));

if is_legacy
	load(fullfile(FILEDIR,listing(1).name),'ephys_data','channels','fs');
	
	ephys.labels=channels;
	ephys.ports='';
	ephys.data=ephys_data;
	ephys.fs=fs;

	audio.fs=ephys.fs;

	clearvars ephys_data channels fs;
else
	load(fullfile(FILEDIR,listing(1).name),'ephys','audio','playback');
end

EPHYS.fs=ephys.fs;
AUDIO.fs=audio.fs;

EPHYS.t=ephys.t;
AUDIO.t=audio.t;

if exist('playback','var') & ~isempty(playback.data)	
	PLAYBACK.t=playback.t;
	PLAYBACK.fs=playback.fs;
end

% number of samples and channels

[samples,channels]=size(ephys.data);
ntrials=length(listing);

disp(['Aggregating trials in: ' FILEDIR]); 
disp(['Number of trials: ' num2str(ntrials)]);

if SLEEPSTATUS

	disp(['Identified sleep data...']);
	parameters.trial_win=parameters.trial_win_sleep;
end

dircount=1;

% extract in a sliding window, generate histograms and trigger mua, sua and lfp processing

savefun=@(filename,agg_ephys,agg_audio,agg_ttl,agg_playback,agg_file_datenum,agg_rms) ...
	save(filename,'agg_ephys','agg_audio','agg_ttl','agg_playback','agg_file_datenum','agg_rms','-v7.3');

for i=0:parameters.trial_win:ntrials
	
	savedir=fullfile(FILEDIR,[SAVEDIR '_' num2str(dircount)]);

	if ~exist(savedir,'dir')
		mkdir(savedir);
	else
		if exist(fullfile(savedir,'.done_extraction'),'file')
			disp(['Directory ' num2str(dircount) ' finished, continuing']);
			dircount=dircount+1;
			continue;
		end
	end

	currtrials=i+1:i+parameters.trial_win;
	currtrials(currtrials>ntrials)=[];

	all_labels=[];
	all_ports='';

	disp(['Identifying channels...']);

	rms_listing={};

	for j=1:length(currtrials)

		tmpfile=fullfile(FILEDIR,listing(currtrials(j)).name);

		if is_legacy
			load(tmpfile,'channels');
			ephys.labels=channels;
			ephys.ports=repmat('A',[1 length(channels)]);
		else
			load(tmpfile,'ephys');
		end

		if isempty(whos('rms','-file',tmpfile)) & ~SLEEPSTATUS
			rms_listing{end+1}=listing(currtrials(j)).name;
		end

		% check for any inconsistency in channel labels

		for k=1:length(ephys.labels)

			label_chk=ephys.labels(k)==all_labels;
			port_chk=ephys.ports(k)==all_ports;

			% loop and if any channels are not included in the channel_label vector, include!

			if ~any(label_chk&port_chk)	
				all_labels=[all_labels ephys.labels(k)];
				all_ports=[all_ports ephys.ports(k)];
			end
		end
	end
	
	if ~SLEEPSTATUS
		disp(['Computing RMS for ' num2str(length(rms_listing)) ' files...']);
		ephys_collect_rms(FILEDIR,rms_listing);
	end

	%[all_labels=sort(all_labels);

	disp(['Found channels:  ']);

	for j=1:length(all_labels)
		fprintf(1,'%i%s ',all_labels(j),all_ports(j));
	end

	fprintf(1,'\n');

	EPHYS.data=zeros(samples,length(currtrials),length(all_labels),'single');
	AUDIO.data=zeros(samples,length(currtrials));
	PLAYBACK.data=zeros(samples,length(currtrials));
	FILE_DATENUM=zeros(1,length(currtrials));
	TTL.data=zeros(samples,length(currtrials));

	if ~SLEEPSTATUS
		
		RMS.standard=zeros(length(all_labels),length(currtrials));
		RMS.robust=zeros(length(all_labels),length(currtrials));
		RMS.labels=all_labels(:);
		RMS.ports=all_ports(:);

	else
		RMS=[];
	
	end
	EPHYS.labels=all_labels;
	EPHYS.ports=all_ports;

	for j=1:length(currtrials)

		disp(['Processing trial ' num2str(currtrials(j))]);

		if is_legacy
		
			load(fullfile(FILEDIR,listing(currtrials(j)).name),'ephys_data','mic_data','channels','fs','start_datenum','ttl_data','rms');
			
			ephys.data=ephys_data;
			ephys.fs=fs;
			ephys.labels=channels;

			audio.data=mic_data;
			audio.fs=fs;

			ttl.data=ttl_data;
			ttl.fs=fs;

			file_datenum=start_datenum;

			clearvars ephys_data mic_data channels fs start_datenum ttl_data;

		else
			load(fullfile(FILEDIR,listing(currtrials(j)).name),'ephys','audio','file_datenum','ttl','playback','rms');
		end

		AUDIO.data(:,j)=audio.data;
		
		if ~exist('ttl','var') | ~isfield(ttl,'data') | isempty(ttl.data) 
			ttl.data=ones(size(audio.data)).*NaN;
		end

		if ~exist('playback','var') |  ~isfield(playback,'data') | isempty(playback.data) 
			playback.data=ones(size(audio.data)).*NaN;
		end

		TTL.data(:,j)=ttl.data;
		PLAYBACK.data(:,j)=playback.data;

		for k=1:length(ephys.labels)

			ch_idx=find(ephys.labels(k)==all_labels);

			if isempty(ch_idx)
				continue;
			end

			EPHYS.data(:,j,ch_idx)=single(ephys.data(:,k));

		end

		if ~SLEEPSTATUS
			for k=1:length(rms.channels)
				ch_idx=find(rms.channels(k)==all_labels);
				RMS.standard(ch_idx,j)=rms.standard(k);
				RMS.robust(ch_idx,j)=rms.robust(k);
			end
		end

		if exist('file_datenum','var')
			FILE_DATENUM(j)=file_datenum;
		else

			% attempt to find the datenum

			tokens=regexpi(listing(currtrials(j)).name,parameters.delimiter,'split');
			FILE_DATENUM(j)=datenum([tokens{end-3} tokens{end-2}],'yymmddHHMMSS');

		end

	end

	% clear any unused data

	if sum(all(isnan(TTL.data)))==length(currtrials)
		TTL.data=[];
	end

	if sum(all(isnan(PLAYBACK.data)))==length(currtrials)
		PLAYBACK.data=[];
	end

	% save the collected data

	savefun(fullfile(savedir,'aggregated_data.mat'),EPHYS,AUDIO,TTL,PLAYBACK,FILE_DATENUM,RMS);
	
	% collect the firing rate for each directory

	if SLEEPSTATUS==1
		fr=ephys_cat_ratetracker(EPHYS,AUDIO,FILE_DATENUM,'savedir',savedir,...
			'savefile','aggregated_stats.mat','channels',EPHYS.labels,'bound',[]);
	else
		fr=ephys_cat_ratetracker(EPHYS,AUDIO,FILE_DATENUM,'savedir',savedir,...
			'savefile','aggregated_stats.mat','channels',EPHYS.labels,'rms',RMS);
	end


	if SLEEPSTATUS

		% write out a file to indicate that this is sleep data

		fid=fopen(fullfile(savedir,'SLEEP_DATA'),'w');
		fclose(fid);

	else
		histogram=ephys_visual_histogram(AUDIO.data,'fs',AUDIO.fs);
		save(fullfile(savedir,'histogram.mat'),'histogram','-v7.3');
	end

	% collect the firing rate for each empty directory, concatenate this later on for showing 
	% movement of RMS, etc. etc., across the day

	% touch signals to multi-unit, single-unit, and fields processing

	if ntrials>parameters.trial_min

		% multi-unit signal

		fid=fopen(fullfile(savedir,'.mua_signal'),'w');
		fclose(fid);

		% signal-unit signal

		fid=fopen(fullfile(savedir,'.sua_signal'),'w');
		fclose(fid);

		% lfp signal

		fid=fopen(fullfile(savedir,'.lfp_signal'),'w');
		fclose(fid);

	end

	if i+parameters.trial_win<=ntrials

		% don't process again if we have a full window

		fid=fopen(fullfile(savedir,'.done_extraction'),'w');
		fclose(fid);

	end

	dircount=dircount+1;


end

clearvars EPHYS AUDIO;

% if there aren't too many directories, go ahead and concatenate the data

dircount=dircount-1; % remove the extra dir

if dircount<=parameters.cat_limit

	disp('Concatenating directory data...');

	ephys_cat_collect_data(FILEDIR);	

	disp('Collecting firing rate across all trials on all channels...');

	% use entire trial for spike threshold detection if sleep data

	disp('Directory 1');

	load(fullfile(FILEDIR,[SAVEDIR '_' num2str(1)],'aggregated_stats.mat'),'fr')
	cat_fr=fr;
	clear fr;

	for i=2:dircount

		disp(['Directory ' num2str(i)]);

		load(fullfile(FILEDIR,[SAVEDIR '_' num2str(i)],'aggregated_stats.mat'),'fr');

		if length(fr)~=length(cat_fr)
			disp('Channel mismatch for firing rates, bailing!');
			return;
		end

		for j=1:length(cat_fr)

			cat_fr(j).time_series=[cat_fr(j).time_series fr(j).time_series];
			cat_fr(j).threshold=[cat_fr(j).threshold fr(j).threshold];
			cat_fr(j).rms.edge=[cat_fr(j).rms.edge fr(j).rms.edge];
			cat_fr(j).rms.song=[cat_fr(j).rms.song fr(j).rms.song];
			cat_fr(j).datenums=[cat_fr(j).datenums fr(j).datenums];

		end

		clear fr;

	end

	% generate plots first, i/o is a little slow here

	disp('Generating plots...');

	ephys_visual_cat_mutrack(cat_fr,'savedir',fullfile(FILEDIR,'ephys_cat'));

	mkdir(fullfile(FILEDIR,'ephys_cat'))
	save(fullfile(FILEDIR,'ephys_cat','cat_fr.mat'),'cat_fr','-v7.3');

	clear cat_fr;


end
