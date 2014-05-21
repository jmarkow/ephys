function frontend_extract_data(FILENAME,DATA,DIRS,IDX,DISP_MINFS,DISP_MAXFS,ISTTL)
%extracts data given a set of points (e.g. TTL or SONG EXTRACTION)
%
%
%
%
%
%

fs=DATA.fs;
ephys_labels=DATA.ephys_labels;
file_datenum=DATA.datenum;

suffix='';
if isttl
	suffix='_ttl';
end

[sonogram_im sonogram_f sonogram_t]=pretty_sonogram(DATA.norm_data,fs,'n',500,'overlap',200,'low',1.5);

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
im_son_to_vec=(length(DATA.norm_data)-350)/t;

for i=1:length(IDX)-1

	startpoint=floor(ttl_pts(ttl_idx(i)+1)-audio_pad*fs);
	endpoint=ceil(ttl_pts(ttl_idx(i+1))+audio_pad*fs);

	if startpoint<1, startpoint=1; end
	if endpoint>length(DATA.norm_data), endpoint=length(DATA.norm_data); end

	sonogram_filename=fullfile(DIRS.image,[ FILENAME suffix '.gif' ]);

	norm_extraction=DATA.norm_data(startpoint:endpoint);
	audio_extraction=DATA.conditioned_data(startpoint:endpoint);
	ephys_extraction=DATA.data(startpoint:endpoint,ephys_channels);
	
	if ~isempty(DATA.ttl_data)
		ttl_extraction=DATA.ttl_data(startpoint:endpoint);
	end

	save_name=[ FILENAME '_chunk_' num2str(i) suffix ];

	sonogram_im(1:10,ceil(startpoint/im_son_to_vec):ceil(endpoint/im_son_to_vec))=63;

	[chunk_sonogram_im chunk_sonogram_f chunk_sonogram_t]=pretty_sonogram(norm_extraction,fs,'n',500,'overlap',350,'low',1.5);

	startidx=max([find(chunk_sonogram_f<=DISP_MINFS)]);
	stopidx=min([find(chunk_sonogram_f>=DISP_MAXFS)]);

	chunk_sonogram_im=chunk_sonogram_im(startidx:stopidx,:);
	chunk_sonogram_im=flipdim(chunk_sonogram_im,1);

	imwrite(uint8(chunk_sonogram_im),colors,fullfile(DIRS.image,[ save_name '.gif']),'gif');
	save(fullfile(DIRS.data,['songdet1_' save_name '.mat']),...
		'ephys_extraction','audio_extraction','ttl_extraction','fs','ephys_labels','file_datenum');

	% normalize audio to write out to wav file

	min_audio=min(norm_extraction(:));
	max_audio=max(norm_extraction(:));

	if min_audio + max_audio < 0
		norm_extraction=norm_extraction./(-min_audio);
	else
		norm_extraction=norm_extraction./(max_audio*(1+1e-3));
	end

	wavwrite(norm_extraction,fs,fullfile(DIRS.wav,[ save_name '.wav']));

end

reformatted_im=im_reformat(sonogram_im,(ceil((length(audio_extraction)/fs)/5)));
imwrite(uint8(reformatted_im),colors,sonogram_filename,'gif');

