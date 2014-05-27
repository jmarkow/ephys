function [T,AMPS,DATA,AUX,PARAMETERS]=frontend_fileparse(FILENAME)
%frontend_readdata returns t amps data and aux
%
%	t
%	time vector
%
%	amps
%	channel ids
%
%	data
%	data matrix
%
%	aux
%	auxiliary channels (if they exist)
%
%	parameters
%	parameters (if they exist)
%
%

[path,filename,ext]=fileparts(FILENAME);

% switch based on ext

switch lower(ext)

	case '.int'
		[T,AMPS,DATA,AUX]=read_intan_data_cli(FILENAME);
		PARAMETERS=[];
	case '.rhd'
		% insert RHD extraction code when ready
	
	% add option for open ephys?
	otherwise
		error('Could not recognize file type %s ', ext);
end
