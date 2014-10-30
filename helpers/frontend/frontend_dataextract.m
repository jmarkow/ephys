function frontend_extract_data(FILENAME,DATA,DIRS,EXT_PTS,DISP_MINFS,DISP_MAXFS,COLORS,SOURCE,DATASAVE,PREFIX,SUFFIX,SKIP)
%extracts data given a set of points (e.g. TTL or SONG EXTRACTION)
%
%
%
%
%
%


fs=DATA.(SOURCE).fs;

if nargin<12
	SKIP=0;
end

if nargin<11 
	SUFFIX='';
end

if nargin<10
	PREFIX='songdet1_';
end

if nargin<9 | isempty(DATASAVE)
	DATASAVE=0;
end

if nargin<8 | isempty(SOURCE)
	SOURCE='audio';
end

[b,a]=ellip(5,.2,40,[300/(fs/2)],'high'); 
if ~isfield(DATA.(SOURCE),'norm_data')
	DATA.(SOURCE).norm_data=filtfilt(b,a,DATA.(SOURCE).data);
end

[sonogram_im sonogram_f sonogram_t]=pretty_sonogram(DATA.(SOURCE).norm_data,fs,'n',500,'overlap',100,'low',1.5);

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
im_son_to_vec=(length(DATA.(SOURCE).norm_data)-100)/t;

data_types={'ephys','ttl','digout','digin','adc','aux','audio','playback'};

savefun=@(filename,datastruct) save(filename,'-struct','datastruct','-v7.3');
sonogram_filename=fullfile(DIRS.image,[ PREFIX FILENAME SUFFIX '.gif' ]);

for i=1:size(EXT_PTS,1)

	startpoint=EXT_PTS(i,1);
	endpoint=EXT_PTS(i,2);

	if startpoint<1 & SKIP
		continue;
	elseif startpoint<1 & ~SKIP
	       	startpoint=1; 
	end
	
	if endpoint>length(DATA.(SOURCE).norm_data) & SKIP
	       continue;
       	elseif endpoint>length(DATA.(SOURCE).norm_data) & ~SKIP
	       endpoint=length(DATA.(SOURCE).norm_data); 
       	end
	
	% cut out the extraction

	EXTDATA=DATA;

	for j=1:length(data_types)
		if ~isempty(EXTDATA.(data_types{j}).data)
			EXTDATA.(data_types{j}).data=EXTDATA.(data_types{j}).data(startpoint:endpoint,:);
			EXTDATA.(data_types{j}).t=EXTDATA.(data_types{j}).t(startpoint:endpoint);
		end
	end

	EXTDATA.(SOURCE).norm_data=EXTDATA.(SOURCE).norm_data(startpoint:endpoint);

	save_name=[ PREFIX FILENAME '_chunk_' num2str(i) SUFFIX ];

	sonogram_im(1:10,ceil(startpoint/im_son_to_vec):ceil(endpoint/im_son_to_vec))=62;

	[chunk_sonogram_im chunk_sonogram_f chunk_sonogram_t]=pretty_sonogram(EXTDATA.(SOURCE).norm_data,fs,'n',500,'overlap',300,'low',1.5);

	startidx=max([find(chunk_sonogram_f<=DISP_MINFS)]);
	stopidx=min([find(chunk_sonogram_f>=DISP_MAXFS)]);

	chunk_sonogram_im=chunk_sonogram_im(startidx:stopidx,:);
	chunk_sonogram_im=flipdim(chunk_sonogram_im,1);
	[f,t]=size(chunk_sonogram_im);
	chunk_im_son_to_vec=(length(EXTDATA.(SOURCE).data)-300)/t;

	if ~isempty(EXTDATA.ttl.data)
		ttl_points=find(EXTDATA.ttl.data>.5);
		ttl_son=round(ttl_points/chunk_im_son_to_vec);
		ttl_son(ttl_son<1|ttl_son>size(chunk_sonogram_im,2))=[];
		chunk_sonogram_im(1:10,round(ttl_son))=62;
	end

	imwrite(uint8(chunk_sonogram_im),colormap([ COLORS '(63)']),fullfile(DIRS.image,[ save_name '.gif']),'gif');

	
	% normalize audio to write out to wav file

	min_audio=min(EXTDATA.(SOURCE).norm_data(:));
	max_audio=max(EXTDATA.(SOURCE).norm_data(:));

	if min_audio + max_audio < 0
		EXTDATA.(SOURCE).norm_data=EXTDATA.(SOURCE).norm_data./(-min_audio);
	else
		EXTDATA.(SOURCE).norm_data=EXTDATA.(SOURCE).norm_data./(max_audio*(1+1e-3));
	end

	wavwrite(EXTDATA.(SOURCE).norm_data,fs,fullfile(DIRS.wav,[ save_name '.wav']));
	
	% remove all data used for plotting (i.e. norm_data)

	for j=1:length(data_types)
		if isfield(EXTDATA.(data_types{j}),'norm_data')
			EXTDATA.(data_types{j})=rmfield(EXTDATA.(data_types{j}),'norm_data');
		end
	end

	%EXTDATA.(SOURCE)=rmfield(EXTDATA.(SOURCE),'norm_data');

	if DATASAVE
		savefun(fullfile(DIRS.data,[ save_name '.mat']),EXTDATA);
	end

	clear EXTDATA;
end

reformatted_im=im_reformat(sonogram_im,(ceil((length(DATA.(SOURCE).data)/fs)/10)));
imwrite(uint8(reformatted_im),colormap([ COLORS '(63)']),sonogram_filename,'gif');

