function [PARAMETERS]=ephys_pipeline_readconfig(FILE)
%
%
%
%
%


fid=fopen(FILE,'r');
readdata=textscan(fid,'%s%[^\n]','commentstyle','#');

PARAMETERS.freq_range={};
idxs=find(strcmp('freq_range',readdata{1}));

for i=1:length(idxs)
	PARAMETERS.freq_range{end+1}=str2num(readdata{2}{idxs(i)});

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





