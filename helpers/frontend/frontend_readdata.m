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

% TODO: change RHA and RHD file handling to pause for .3 seconds to see if file is being written,
% should speed things up considerably

switch lower(ext)

	case '.int'
		
		[T,AMPS,DATA,AUX]=read_intan_data_cli(FILENAME);
		PARAMETERS=[];
	
	case '.rhd'

		


		% pack most of the data into a misc structure, the rest can be unpacked into
		% T AMPS DATA and AUX as before
	
		% insert RHD extraction code when ready

	case '.continuous'

		% open ephys extraction code here

	% add option for open ephys?
	otherwise
		error('Could not recognize file type %s ', ext);
end
