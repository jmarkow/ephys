function [STORE_EPHYS STORE_AUDIO STORE_PLAYBACK STORE_DATENUM]=ephys_collect_data(DIR,varargin)
% run in the same directory as ephys_* directories created by the songextraction scripts, this
% will create another directory with data aggregated across all directories, useful to check 
% firing rates, singing rates, etc. across directories
%
% see also ephys_pipeline_songextract_aggregate.m
%

nparams=length(varargin);

if mod(nparams,2)>0
	error('ephysPipeline:argChk','Parameters must be specified as parameter/value pairs!');
end

ephys_fs=10e3; % downsampling rate
downfilt=5e3; % corner frequency of anti-aliasing filter
savedir='ephys_cat'; % save directory

% remove eps generation, too slow here...

for i=1:2:nparams
	switch lower(varargin{i})
		case 'ephys_fs'
			ephys_fs=varargin{i+1};
		case 'downfilt'
			downfilt=varargin{i+1};
		case 'savedir'
			savedir=varargin{i+1};
	end
end

if nargin<1 | isempty(DIR)
	DIR=pwd;
end

savefun=@(filename,cat_ephys,cat_audio,cat_playback,cat_file_datenum,cat_rms) ...
	save(filename,'cat_ephys','cat_audio','cat_playback','cat_file_datenum','cat_rms','-v7.3');

% collect and downsample all ephys, audio, and playback data

[status,result]=unix(['find ' DIR ' -type f -name "aggregated_data.mat"']);
ephys_file=regexp(result,'\n','split');
ephys_file(end)=[];

disp('Loading data file 1...');

load(ephys_file{1},'agg_ephys','agg_file_datenum','agg_playback','agg_audio','agg_rms');

[b,a]=butter(4,downfilt/(agg_ephys.fs/2),'low');

disp('Downsampling...');

agg_ephys.data=single(downsample(filtfilt(b,a,double(agg_ephys.data)),agg_ephys.fs/ephys_fs));

STORE_EPHYS=agg_ephys;
STORE_DATENUM=agg_file_datenum;
STORE_PLAYBACK=agg_playback;
STORE_AUDIO=agg_audio;
STORE_EPHYS.fs=ephys_fs;
STORE_EPHYS.t=downsample(STORE_EPHYS.t,agg_ephys.fs/ephys_fs);
STORE_RMS=agg_rms;

clearvars agg_ephys agg_file_datenum agg_playback agg_audio;

for i=2:length(ephys_file)

	disp(['Loading data file ' num2str(i) '...']);

	load(ephys_file{i},'agg_ephys','agg_file_datenum','agg_playback','agg_audio','agg_rms');

	% get song range

	disp('Downsampling and storing data...');

	agg_ephys.data=single(downsample(filtfilt(b,a,double(agg_ephys.data)),agg_ephys.fs/ephys_fs));
	STORE_EPHYS.data=[STORE_EPHYS.data agg_ephys.data];

	clear agg_ephys;

	STORE_DATENUM=[STORE_DATENUM agg_file_datenum];
	STORE_PLAYBACK.data=[STORE_PLAYBACK.data agg_playback.data];
	STORE_AUDIO.data=[STORE_AUDIO.data agg_audio.data];

	if ~isempty(STORE_RMS)
		STORE_RMS.allpts.standard=[STORE_RMS.allpts.standard agg_rms.allpts.standard];
		STORE_RMS.allpts.robust=[STORE_RMS.allpts.robust agg_rms.allpts.robust];
		STORE_RMS.songpts.standard=[STORE_RMS.songpts.standard agg_rms.songpts.standard];
		STORE_RMS.songpts.robust=[STORE_RMS.songpts.robust agg_rms.songpts.robust];
		STORE_RMS.silentpts.standard=[STORE_RMS.silentpts.standard agg_rms.silentpts.standard];
		STORE_RMS.silentpts.robust=[STORE_RMS.silentpts.robust agg_rms.silentpts.robust];
	end
	clearvars agg_file_datenum agg_playback agg_audio agg_rms;

end

if ~isempty(savedir)
	mkdir(fullfile(DIR,savedir))
	savefun(fullfile(DIR,savedir,'cat_data.mat'),STORE_EPHYS,STORE_AUDIO,STORE_PLAYBACK,STORE_DATENUM,STORE_RMS);
end
