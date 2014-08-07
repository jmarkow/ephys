function [DATASTRUCT]=frontend_readdata(FILENAME)
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

% pre-allocating everything the stupid way

DATASTRUCT.ephys.labels=[];
DATASTRUCT.ephys.ports=[];
DATASTRUCT.ephys.data=[];
DATASTRUCT.ephys.fs=[];
DATASTRUCT.ephys.t=[];

DATASTRUCT.aux.labels=[];
DATASTRUCT.aux.ports=[];
DATASTRUCT.aux.data=[];
DATASTRUCT.aux.fs=[];
DATASTRUCT.aux.t=[];

DATASTRUCT.digin.labels=[];
DATASTRUCT.digin.data=[];
DATASTRUCT.digin.fs=[];
DATASTRUCT.digin.t=[];

DATASTRUCT.digout.labels=[];
DATASTRUCT.digout.data=[];
DATASTRUCT.digout.fs=[];
DATASTRUCT.digout.t=[];

DATASTRUCT.adc.labels=[];
DATASTRUCT.adc.data=[];
DATASTRUCT.adc.fs=[];
DATASTRUCT.adc.t=[];

switch lower(ext)

	case '.int'
		
		[t,amps,data,aux]=read_intan_data_cli(FILENAME);

		DATASTRUCT.ephys.t=t;
		DATASTRUCT.ephys.labels=amps;
		DATASTRUCT.ephys.data=data;
		DATASTRUCT.ephys.fs=round(1./(t(2)-t(1)));

		DATASTRUCT.aux.t=DATASTRUCT.ephys.t;
		DATASTRUCT.aux.labels=1:size(AUX,2);
		DATASTRUCT.aux.data=aux;
		DATASTRUCT.aux.fs=DATASTRUCT.ephys.fs;

	case '.rhd'

		[amp,aux_input,parameters,notes,supply_voltage,adc,dig_in,dig_out,temp_sensor]=...
			read_intan_data_cli_rhd2000(FILENAME);

		DATASTRUCT.ephys.labels=cat(1,amp.channels(:).native_order)';
		DATASTRUCT.ephys.ports=cat(1,amp.channels(:).port_prefix)';	
		DATASTRUCT.ephys.data=amp.data';
		DATASTRUCT.ephys.fs=parameters.amplifier_sample_rate;
		DATASTRUCT.ephys.t=amp.t(:);

		DATASTRUCT.parameters=parameters;
		DATASTRUCT.notes=notes;
	
		if ~isempty(dig_in)		
			DATASTRUCT.digin.labels=cat(1,dig_in.channels(:).native_order)';
			DATASTRUCT.digin.data=dig_in.data';
			DATASTRUCT.digin.fs=parameters.board_dig_in_sample_rate;
			DATASTRUCT.digin.t=dig_in.t(:);
		end

		if ~isempty(dig_out) 
			DATASTRUCT.digout.labels=cat(1,dig_out.channels(:),native_order)';
			DATASTRUCT.digout.data=dig_out.data';
			DATASTRUCT.digout.t=dig_out.t(:);
		end
		
		if ~isempty(aux_input)	
			DATASTRUCT.aux.labels=cat(1,aux_input.channels(:).native_order)';
			DATASTRUCT.aux.ports=cat(1,aux_input.channels(:).port_prefix)';
			DATASTRUCT.aux.data=aux_input.data';
			DATASTRUCT.aux.fs=parameters.aux_input_sample_rate;
			DATASTRUCT.aux.t=aux_input.t(:);
		end

		if ~isempty(adc)
		    DATASTRUCT.adc.labels=cat(1,adc.channels(:).native_order)';
		    DATASTRUCT.adc.data=adc.data';
			DATASTRUCT.adc.fs=parameters.board_adc_sample_rate;
			DATASTRUCT.adc.t=adc.t(:);
	   	end

		% pack most of the data into a misc structure, the rest can be unpacked into
		% T AMPS DATA and AUX as before
	

	case '.continuous'

		% open ephys extraction code here
		% add option for open ephys?

	otherwise
		error('Could not recognize file type %s ', ext);
end
