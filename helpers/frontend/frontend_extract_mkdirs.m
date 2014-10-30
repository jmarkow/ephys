function frontend_extract_mkdirs(FOLDER,IMAGE,WAVE,DATA,TTL,PLAYBACK)
%frontend_mkdirs makes directories for the frontend
%
%
%
%

image_dir=fullfile(FOLDER,IMAGE);
wav_dir=fullfile(FOLDER,WAVE);
data_dir=fullfile(FOLDER,DATA);

if TTL
	image_dir_ttl=fullfile(FOLDER,[IMAGE '_ttl']);
	wav_dir_ttl=fullfile(FOLDER,[WAVE '_ttl']);
	data_dir_ttl=fullfile(FOLDER,[DATA '_ttl']);
end


if PLAYBACK
	image_dir_pback=fullfile(FOLDER,[IMAGE '_pback']);
	wav_dir_pback=fullfile(FOLDER,[WAVE '_pback']);
	data_dir_pback=fullfile(FOLDER,[DATA '_pback']);
end

if ~exist(image_dir,'dir')
	mkdir(image_dir);
end

if ~exist(wav_dir,'dir');
	mkdir(wav_dir);
end

if ~exist(data_dir,'dir');
	mkdir(data_dir);
end

if TTL
	if ~exist(image_dir_ttl,'dir')
		mkdir(image_dir_ttl);
	end

	if ~exist(wav_dir_ttl,'dir');
		mkdir(wav_dir_ttl);
	end

	if ~exist(data_dir_ttl,'dir');
		mkdir(data_dir_ttl);
	end
end

if PLAYBACK
	if ~exist(image_dir_pback,'dir')
		mkdir(image_dir_pback);
	end

	if ~exist(wav_dir_pback,'dir');
		mkdir(wav_dir_pback);
	end

	if ~exist(data_dir_pback,'dir');
		mkdir(data_dir_pback);
	end
end

