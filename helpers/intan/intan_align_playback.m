function intan_align_playback(DIR,INTAN_DIR,varargin)
%
%
%
%
%

mic_trace=1;
ttl_trace=2;
aux_trace=6;
intan_fs=25e3;
fs=40e3;
debug=0;

nparams=length(varargin);

if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs');
end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'mic_trace'
			mic_trace=varargin{i+1};
		case 'ttl_trace'
			ttl_trace=varargin{i+1};
		case 'aux_trace'
			aux_trace=varargin{i+1};
		case 'spect_thresh'
			spect_thresh=varargin{i+1};
		case 'sr'
			SR=varargin{i+1};
		case 'colors'
			colors=varargin{i+1};
		case 'debug'
			debug=varargin{i+1};
	end
end


% load select data, extract mic and ttl data

[filename,pathname]=uigetfile('*.mat','Pick a sound file to extract the template from',fullfile(DIR));
load(fullfile(pathname,filename),'data');

mic_data=data.voltage(:,mic_trace);
ttl_data=data.voltage(:,ttl_trace);

% 1) select sample from NIDAQ (use spectro gui to select sound and corresponding TTL)
% 2) resample to INTAN rate (25e3)
% 3) match digital signal
% 4) read directory of Intan files, match and extract

[templatesound,extractidxs]=spectro_navigate(mic_data);
ttl_template=ttl_data(extractidxs(1):extractidxs(2));

% resample to Intan rate

resample_factor=intan_fs/fs;
[p,q]=rat(resample_factor);

res_mic_data=resample(audio_extraction,p,q);
res_mic_ttl=double(resample(ttl_extraction,p,q)>.5);

figure();

ax(1)=subplot(2,1,1);
plot(res_mic_data);
ax(2)=subplot(2,1,2);
plot(res_mic_ttl);

linkaxes(ax,'x');

% with TTL, template, blow through all Intan files and extract hits via xcorr

if isempty(INTAN_DIR)
	INTAN_DIR=uigetdir(pwd,'Select the directory with matching Intan files...');
end

% extract all files with the Intan extension, .int

pre_file_list=dir(fullfile(intan_dir,'*.int')); 
file_list={};

for i=1:length(pre_file_list)
	file_list{i}=fullfile(intan_dir,pre_file_list(i).name);
end

% length of the extraction

template_length=length(res_mic_ttl);

for i=1:length(file_list)

	% extract intan data

	[t,channels,intan_data,aux] = read_intan_data_cli(file_list{i});

	intan_ttl=aux(:,aux_trace);

	[align_score,lags]=xcorr(intan_ttl,res_mic_ttl);

	% peaks? align_score>x

	peaks=find(align_score>threshold);

	% time between peaks?

	peakstoextract=peaks((diff(peaks)<min_spacing)+1)=[];

	for j=1:length(peakstoextract)

		% indices for extraction from Intan

		alignpoint=lags(peakstoextract(j));

		stoppoint=alignpoint+templatelength;
		
		ephys_extraction=intan_data(startidx:stoppoint);
		ttl_extraction=intan_ttl(startidx:stoppoint);





	end

end

% 
%

