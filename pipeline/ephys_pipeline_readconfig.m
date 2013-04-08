function [PARAMETERS]=ephys_pipeline_readconfig(FILE)
%
%
%
%
%


fid=fopen(FILE,'r');
readdata=textscan(fid,'%s%[^\n]','commentstyle','#');

idxs=find(strcmp('freq_range',readdata{1}));

if ~isempty(idxs)

	PARAMETERS.freq_range={};

	for i=1:length(idxs)
		PARAMETERS.freq_range{end+1}=str2num(readdata{2}{idxs(i)});
	end

end

readdata{1}(idxs)=[];
readdata{2}(idxs)=[];

idxs=find(strcmp('spike_weights',readdata{1}));

if ~isempty(idxs)
	PARAMETERS.spike_weights=[];
	PARAMETERS.spike_weights=str2num(readdata{2}{idxs(1)});
end

readdata{1}(idxs)=[];
readdata{2}(idxs)=[];

idxs=find(strcmp('spike_freq_range',readdata{1}));

if ~isempty(idxs)
	PARAMETERS.spike_freq_range=[];
	PARAMETERS.spike_freq_range=str2num(readdata{2}{idxs(1)});
end

readdata{1}(idxs)=[];
readdata{2}(idxs)=[];

idxs=find(strcmp('spike_window',readdata{1}));

if ~isempty(idxs)
	PARAMETERS.spike_window=[];
	PARAMETERS.spike_window=str2num(readdata{2}{idxs(1)});
end

readdata{1}(idxs)=[];
readdata{2}(idxs)=[];


idxs=find(strcmp('spike_cluststart',readdata{1}));

if ~isempty(idxs)
	PARAMETERS.spike_cluststart=[];
	PARAMETERS.spike_cluststart=str2num(readdata{2}{idxs(1)});
end

readdata{1}(idxs)=[];
readdata{2}(idxs)=[];

for i=1:length(readdata{1})
	PARAMETERS.(readdata{1}{i})=str2double(readdata{2}{i});

	% if this returns NaN, then the parameter value was specified as a string,
	% read as such

	if isnan(PARAMETERS.(readdata{1}{i}))

		if strcmp(readdata{2}{i},'[]')
			PARAMETERS.(readdata{1}{i})=str2num(readdata{2}{i});
		else
			PARAMETERS.(readdata{1}{i})=readdata{2}{i};
		end

	end


end





