function [T,AMPS,DATA,AUX_AMPS,AUX,PARAMETERS,DIG,ADC]=frontend_fileparse(FILENAME)
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

T=[];
AMPS=[];
DATA=[];
AUX_AMPS=[];
AUX=[];
PARAMETERS=[];
DIG=[];
ADC=[];

switch lower(ext)

	case '.int'
		
		[T,AMPS,DATA,AUX]=read_intan_data_cli(FILENAME);	
		AUX_AMPS=1:size(AUX,2);
		
	case '.rhd'

		[amp,aux_input,PARAMETERS,notes,supply_voltage,ADC,dig_in,dig_out,temp_sensor]=...
			read_intan_data_cli_rhd2000(FILENAME);

		T=amp.t;
		DATA=amp.data;
		AMPS=cat(1,amp.channels(:).native_order)+1;
		PARAMETERS.notes=notes;
		PARAMETERS.amps=AMPS;
		DIG.IN=dig_in;
		DIG.OUT=dig_out;

		AUX=aux_input.data;
		AMPS_AUX=cat(1,aux_input.channels(:).native_order)+1;

		% pack most of the data into a misc structure, the rest can be unpacked into
		% T AMPS DATA and AUX as before
	
		% insert RHD extraction code when ready

	case '.continuous'

		% open ephys extraction code here

	% add option for open ephys?
	otherwise
		error('Could not recognize file type %s ', ext);
end
