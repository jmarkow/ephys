function [agg_ephys,histogram,agg_file_datenum]=ephys_pipeline_dataload(PROCDIR,FILENAME)
%
%
%

if ~exist(fullfile(PROCDIR,'histogram.mat'),'file')
	histogram=[];
end

try	
	tmp=whos('-file',fullfile(PROCDIR,FILENAME));
	is_legacy=ismember('ephys.data',{tmp(:).name});

	if is_legacy
		load(fullfile(PROCDIR,FILENAME),'EPHYS_DATA','CHANNELS','START_agg_file_datenum');

		if exist(fullfile(PROCDIR,'histogram.mat'),'file')
			load(fullfile(PROCDIR,'histogram.mat'),'HISTOGRAM');
		end

		agg_ephys.data=EPHYS_DATA;
		agg_ephys.labels=CHANNELS;
		agg_file_datenum=START_agg_file_datenum;
		clearvars EPHYS_DATA CHANNELS HISTOGRAM START_agg_file_datenum;

	else
		load(fullfile(PROCDIR,FILENAME),'agg_ephys','agg_file_datenum');

		if exist(fullfile(PROCDIR,'histogram.mat'),'file')
			load(fullfile(PROCDIR,'histogram.mat'),'histogram');
		end

	end
catch
	try	
		tmp=whos('-file',fullfile(PROCDIR,FILENAME));
		is_legacy=ismember('ephys.data',{tmp(:).name});

		disp('Pausing for 60 seconds and will retry');
		pause(60);

		if is_legacy
			
			load(fullfile(PROCDIR,FILENAME),'EPHYS_DATA','CHANNELS','START_agg_file_datenum');

			if exist(fullfile(PROCDIR,'histogram.mat'),'file')
				load(fullfile(PROCDIR,'histogram.mat'),'histogram');
			end

			agg_ephys.data=ephys.data;
			agg_ephys.labels=EPHYS_CHANNELS;
			agg_file_datenum=start_datenum;

			clearvars EPHYS_DATA CHANNELS HISTOGRAM START_agg_file_datenum;


		else
			load(fullfile(PROCDIR,FILENAME),'agg_ephys','agg_file_datenum');

			if exist(fullfile(PROCDIR,'histogram.mat'),'file')
				load(fullfile(PROCDIR,'histogram.mat'),'histogram');
			end
		end
	catch

		disp('Could not properly load files, bailing!');
		return;
	end

end
