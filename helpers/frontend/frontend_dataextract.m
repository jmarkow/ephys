function frontend_extract_data(FILENAME,DATA,DIRS,EXT_PTS,DISP_MINFS,DISP_MAXFS,ISTTL,COLORS)
%extracts data given a set of points (e.g. TTL or SONG EXTRACTION)
%
%
%
%
%
%

fs=DATA.audio.fs;

suffix='';
if ISTTL
	suffix='_ttl';
end

[sonogram_im sonogram_f sonogram_t]=pretty_sonogram(DATA.audio.norm_data,fs,'n',500,'overlap',100,'low',1.5);

startidx=max([find(sonogram_f<=DISP_MINFS)]);

if isempty(startidx)
	startidx=1;
end

stopidx=min([find(sonogram_f>=DISP_MAXFS)]);

if isempty(stopidx)
	stopidx=length(sonogram_f);
end

sonogram_im=sonogram_im(startidx:stopidx,:);
sonogram_im=flipdim(sonogram_im,1);
[f,t]=size(sonogram_im);
im_son_to_vec=(length(DATA.audio.norm_data)-350)/t;

data_types={'ephys','ttl','digout','digin','adc','aux','audio','playback'};

savefun=@(filename,datastruct) save(filename,'-struct','datastruct','-v7.3');

for i=1:size(EXT_PTS,1)

	startpoint=EXT_PTS(i,1);
	endpoint=EXT_PTS(i,2);

	if startpoint<1, startpoint=1; end
	if endpoint>length(DATA.audio.norm_data), endpoint=length(DATA.audio.norm_data); end

	sonogram_filename=fullfile(DIRS.image,[ FILENAME suffix '.gif' ]);

	% cut out the extraction

	EXTDATA=DATA;

	for j=1:length(data_types)
		if ~isempty(EXTDATA.(data_types{j}).data)
			EXTDATA.(data_types{j}).data=EXTDATA.(data_types{j}).data(startpoint:endpoint,:);
			EXTDATA.(data_types{j}).t=EXTDATA.(data_types{j}).t(startpoint:endpoint);
		end
	end

	EXTDATA.audio.norm_data=EXTDATA.audio.norm_data(startpoint:endpoint);

	save_name=[ FILENAME '_chunk_' num2str(i) suffix ];

	sonogram_im(1:10,ceil(startpoint/im_son_to_vec):ceil(endpoint/im_son_to_vec))=62;

	[chunk_sonogram_im chunk_sonogram_f chunk_sonogram_t]=pretty_sonogram(EXTDATA.audio.norm_data,fs,'n',500,'overlap',300,'low',1.5);

	startidx=max([find(chunk_sonogram_f<=DISP_MINFS)]);
	stopidx=min([find(chunk_sonogram_f>=DISP_MAXFS)]);

	chunk_sonogram_im=chunk_sonogram_im(startidx:stopidx,:);
	chunk_sonogram_im=flipdim(chunk_sonogram_im,1);

	imwrite(uint8(chunk_sonogram_im),colormap([ COLORS '(63)']),fullfile(DIRS.image,[ save_name '.gif']),'gif');

	
	% normalize audio to write out to wav file

	min_audio=min(EXTDATA.audio.norm_data(:));
	max_audio=max(EXTDATA.audio.norm_data(:));

	if min_audio + max_audio < 0
		EXTDATA.audio.norm_data=EXTDATA.audio.norm_data./(-min_audio);
	else
		EXTDATA.audio.norm_data=EXTDATA.audio.norm_data./(max_audio*(1+1e-3));
	end

	wavwrite(EXTDATA.audio.norm_data,fs,fullfile(DIRS.wav,[ save_name '.wav']));

	EXTDATA.audio=rmfield(EXTDATA.audio,'norm_data');
	savefun(fullfile(DIRS.data,['songdet1_' save_name '.mat']),EXTDATA);

	clear EXTDATA;

end

reformatted_im=im_reformat(sonogram_im,(ceil((length(DATA.audio.data)/fs)/10)));
imwrite(uint8(reformatted_im),colormap([ COLORS '(63)']),sonogram_filename,'gif');

