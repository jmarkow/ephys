function file2nex(varargin)
% frontend takes an Intan RHD or RHD file (*.rhd|*.rha) and converts it to 
% a continuously digitized binary file suitable for importing into Plexon 
% Offline sorter.  To run, navigate to the directory with the Intan file
% then run frontend; in the command window.  The script will prompt for the 
% Intan file(s) to convert, then it will ask for the channels you want to include in the
% new file.  The output consists of two files, one with a .bin extension (the binary data
% ) and another with an .ofi extension (the parameter settings for Plexon).
%
%
%

buffersize=1e6;
[filename,pathname]=uigetfile({'*.rhd';'*.mat'},'Pick an intan or MATLAB file to extract the template from',pwd,'MultiSelect','on');

% determine whether or nor we have multiple files

nfiles=length(filename);

% specify filename using uiputfile

[storefile,storepath]=uiputfile('*.nex','New filename');

filenames=filename;
clearvars filename;

% read in first data file

[~,~,ext]=fileparts(filenames{1});

switch lower(ext)
	case '.rhd'
		datastruct=frontend_readdata(fullfile(pathname,filenames{1}));
	case '.mat'
		datastruct=load(fullfile(pathname,filenames{1}));
		is_legacy=check_legacy(fullfile(pathname,filenames{1}));

	otherwise
end

if is_legacy
	datastruct.ephys.data=datastruct.ephys_data;
	datastruct.ephys.fs=datastruct.fs;
	datastruct.ephys.ports=repmat('A',[1 size(datastruct.ephys.data,2)]);
	datastruct.ephys.labels=datastruct.channels;
	datastruct=rmfield(datastruct,'ephys_data');
	datastruct=rmfield(datastruct,'fs');
	datastruct=rmfield(datastruct,'channels');
end

fs=datastruct.ephys.fs;
% build menu to select ports and channels

dialogtxt={};

for i=1:length(datastruct.ephys.labels)
	dialogtxt{i}=sprintf('Ch %i (Port %s)',datastruct.ephys.labels(i),datastruct.ephys.ports(i));
end

[selection,status]=listdlg('PromptString','Select channels to include',...
	'SelectionMode','multi',...
	'ListString',dialogtxt);
portselection=datastruct.ephys.ports(selection);
channelselection=datastruct.ephys.labels(selection);

% format the data for writing

disp('Reformatting data (this may take a minute)..');

store_datavec=zeros(buffersize,length(channelselection),'double'); % signed int16, initialize to reasonable number, cut off unused samples
counter=0;
trial_stamps=[];

for ii=1:nfiles

	disp( [ 'Grabbing data from: ' filenames{ii} ] );
	
	[~,~,ext]=fileparts(filenames{1});

	if ii>1
		switch lower(ext)
			case '.rhd'
				datastruct=frontend_readdata(fullfile(pathname,filenames{ii}));
			case '.mat'
				is_legacy=check_legacy(fullfile(pathname,filenames{ii}));
				datastruct=load(fullfile(pathname,filenames{ii}));
			otherwise
		end

		% convert legacy data if necessary


		if is_legacy
			datastruct.ephys.data=datastruct.ephys_data;
			datastruct.ephys.fs=datastruct.fs;
			datastruct.ephys.ports=repmat('A',[1 size(datastruct.ephys.data,2)]);
			datastruct.ephys.labels=datastruct.channels;
			datastruct=rmfield(datastruct,'ephys_data');
			datastruct=rmfield(datastruct,'fs');
			datastruct=rmfield(datastruct,'channels');
		end


	end

	datasamples=size(datastruct.ephys.data,1);
	nchannels=length(selection);
	datalen=datasamples*length(selection);
	datavec=zeros(1,datalen);

	% get indices

	chidx=[];
	for i=1:length(datastruct.ephys.labels)
		flag1=datastruct.ephys.labels(i)==channelselection;
		flag2=datastruct.ephys.ports(i)==portselection;	
		chidx=[chidx find(flag1&flag2)];
	end

	% stitch multiple files if necessary

	trial_stamps=[trial_stamps counter+1]; % every trial starts at COUNTER+1

	for j=1:length(chidx)
		store_datavec(counter+1:counter+datasamples,j)=datastruct.ephys.data(:,chidx(j));
	end


	counter=counter+datasamples;
	clearvars datastruct;

end

% delete unnecessary data

if counter<buffersize
	store_datavec(counter:buffersize,:)=[];
end

disp('Writing data...');

nexfile=nexCreateFileData(40e3);

for i=1:size(store_datavec,2)
	nexfile=nexAddContinuous(nexfile,1/fs,fs,store_datavec(:,i),...
		[ 'CH ' num2str(channelselection(i)) ' (Port ' portselection(i) ')' ]);
end

nexfile=nexAddEvent(nexfile,trial_stamps/fs,'trialtimes');
writeNexFile(nexfile,fullfile(storepath,storefile));

