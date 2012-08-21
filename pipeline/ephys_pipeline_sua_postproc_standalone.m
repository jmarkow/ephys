function ephys_pipeline_sua_postproc_standalone(UNITFILE,CONFIG)
%
%
%
%
%
%

global_parameters=ephys_pipeline_readconfig(CONFIG);
unit_parameters=ephys_pipeline_readconfig(UNITFILE);

% run su_lfp_coherence_spect
% run su_lfp_coherence_tf
% run su_lfp_coherence_ifr

% perhaps have options to plot other spike stats, autocorrelation, etc.

[path,file,ext]=fileparts(UNITFILE);
tokens=regexp(path,filesep,'split');

if length(tokens)<12
	error('ephysPipeline:singleunitpostproc:wrongdir','Unit file in wrong directory for post processing.');
end

% TODO:  deal with non-default directory structures (or, simply enforce the standard ones more rigorously)

extractionstring=tokens{end-2};
datestring=tokens{end-4};
recstring=tokens{end-5};
birdidstring=tokens{end-6};

bookkeeping_dir=fullfile(filesep,tokens{1:end-6},recstring,global_parameters.bookkeeping_dir);
savedir=fullfile(bookkeeping_dir,unit_parameters.cellid,[ datestring ' (' extractionstring ')']);

disp(['Post-processing ' unit_parameters.cellid]);
disp(['Will save in directory ' savedir]);

filedir=path;

if exist(fullfile(path,'histogram.mat'))
	load(fullfile(path,'histogram.mat'),'HISTOGRAM');
else
	error('ephysPipeline:singleunitpostproc:nohistogram','No histogram in same directory as unit file.');
end

% copy any single unit rasters for safe keeping...

if ~exist(fullfile(savedir,'raster'),'dir')
	mkdir(fullfile(savedir,'raster'));
end

try
	copyfile(fullfile(filedir,'sua',['ephys_sua_freqrange*electrode_' num2str(unit_parameters.channel) ...
		'_raster_cluster_' num2str(unit_parameters.cluster) '*' ]),fullfile(savedir,'raster'));
	copyfile(fullfile(filedir,'sua',['ephys_sua_freqrange*electrode_' num2str(unit_parameters.channel) ...
		'_stats_cluster_' num2str(unit_parameters.cluster) '*' ]),fullfile(savedir,'raster'));
	copyfile(fullfile(filedir,'sua',['sua_channels*' num2str(unit_parameters.channel) '.mat']),fullfile(savedir,'raster'));
catch
	warning('ephysPipeline:singleunitpostproc:nocopy','Could not copy single unit raster file');
end

% multiple numbers in config file are read in as strings

if length(unit_parameters.lfp_channels)>1
	lfp_channels=str2num(unit_parameters.lfp_channels);
else
	lfp_channels=unit_parameters.lfp_channels;
end

tmpvec=[0:25:100];
bands_to_check=[1 tmpvec(2:end)];


for i=lfp_channels

	fprintf('SU Channel:\t%i\n',unit_parameters.channel);
	fprintf('SU Cluster:\t%i\n',unit_parameters.cluster);
	fprintf('LFP Channel:\t%i\n',i);

	[abscoh,m_coherence_null,m_coherence_err,freqs]=ephys_su_lfp_coherence_spect(i,unit_parameters.channel,...
		unit_parameters.cluster,'filedir',filedir,'savedir',savedir,'alpha',global_parameters.coh_alpha);

	% perform checks

	for j=1:length(bands_to_check)-1

		fprintf('FS band:\t%g-%g\n',bands_to_check(j),bands_to_check(j+1));


		startidx=max(find(freqs<=bands_to_check(j)));

		if isempty(startidx)
			startidx=1;
		end

		stopidx=min(find(freqs>=bands_to_check(j+1)));

		if isempty(stopidx)
			stopidx=length(freqs);
		end

		minline=abscoh(startidx:stopidx)-m_coherence_err;
		
		if any(minline>m_coherence_null)

			% if any point crosses significance, generate IFR plots...

			disp(['Computing IFR-triggered LFPs for the following frequency band:  '...
			       	num2str(bands_to_check(j:j+1))]);

			ephys_su_lfp_coherence_ifr(i,unit_parameters.channel,unit_parameters.cluster,'filedir',filedir,...
				'savedir',savedir,'freq_range',[bands_to_check(j) bands_to_check(j+1)]);
			ephys_su_lfp_coherence_spikestats(i,unit_parameters.channel,unit_parameters.cluster,'filedir',filedir,...
				'savedir',savedir,'freq_range',[bands_to_check(j) bands_to_check(j+1)]);

		end


	end

	ephys_su_lfp_coherence_tf(i,unit_parameters.channel,unit_parameters.cluster,HISTOGRAM,'filedir',filedir,'savedir',savedir);

end





