function [BIRDID,RECID,MICTRACE,MICSOURCE,MICPORT,PORTS,TTLTRACE,TTLSOURCE,DATENUM]=frontend_fileparse(FILENAME,DELIM,FMT,DATEFMT)
%frontend_fileparse
%
%	frontend_fileparse(FILENAME,DELIM,FMT)
%
%	FILENAME
%
%	DELIM
%
%	FMT
%
%	

% B=subject name
% I=recording ID
% M=microphone
% T=ttl
% D=date

if nargin<4 | isempty(DATEFMT), DATEFMT='yymmddHHMMSS'; end
if nargin<3 | isempty(FMT), FMT='bimdd'; end
if nargin<2 | isempty(DELIM), DELIM='\_'; end
if nargin<1 | isempty(FILENAME), error('Need filename to continue!'); end

BIRDID=[];
RECID=[];
MICTRACE=[];
MICSOURCE='';
MICPORT='';

TTLTRACE=[];
TTLSOURCE='';

DATENUM=[];
PORTS='';

port_labels='abcd';

birdtoken=strfind(lower(FMT),'b');
rectoken=strfind(lower(FMT),'i');
mictoken=strfind(lower(FMT),'m');
ttltoken=strfind(lower(FMT),'t');
datetoken=strfind(lower(FMT),'d');
porttoken=strfind(lower(FMT),'p');

[path,filename,ext]=fileparts(FILENAME);

tokens=regexpi(filename,DELIM,'split');

% first token should be bird number

if ~isempty(birdtoken)
	BIRDID=tokens{birdtoken};
else
	BIRDID='birdname';
end

% second should be recording tag (normally nucleus)

if ~isempty(rectoken)
	RECID=tokens{rectoken};
else
	RECID='RECID';
end

% third should be mic trace

if ~isempty(mictoken)
	string=tokens{mictoken};
	[mictokens,startpoint,endpoint]=regexpi(string,'\d+','match');
	
	if isempty(mictokens)
		error('Could not find mic trace at token %g for file %s', mictoken,FILENAME);
	end

	MICTRACE=str2num(mictokens{1});

	if length(string)>endpoint
		tmp=string(endpoint+1:end);

		if strcmp(lower(tmp(1)),'m')
			MICSOURCE='m';
		elseif strcmp(lower(tmp(1:2)),'au') 
			MICSOURCE='a';
		elseif strcmp(lower(tmp(1:2)),'ad')
			MICSOURCE='c';
		else
			warning('Did not understand mic source %s setting to main.',tmp);
			MICSOURCE='m';
		end
	else
		warning('No microphone source given, set to main');
		MICSOURCE='m';
	end

	% need a port if we're using ephys channels for the mic

	[port,startpoint,endpoint]=regexpi(string,'port','match');

	if strcmp(MICSOURCE,'m')
		if ~isempty(port)
			MICPORT=string(endpoint+1);
		else
			warning('Could not find microphone port, setting to A');
			MICPORT='A';
		end
	end
end


if ~isempty(ttltoken)
	string=tokens{ttltoken};
	[ttltokens,startpoint,endpoint]=regexpi(string,'\d+','match');

	if isempty(ttltokens)
		error('Could not find ttl trace at token %g for file %s', ttltoken,FILENAME);
	end

	TTLTRACE=str2num(ttltokens{1});	

	if length(string)>endpoint
		tmp=string(endpoint+1:end);

		if strcmp(lower(tmp(1)),'m')
			TTLSOURCE='m';
		elseif strcmp(lower(tmp(1:2)),'au') 
			TTLSOURCE='a';
		elseif strcmp(lower(tmp(1:2)),'ad')
			TTLSOURCE='c';
		elseif strcmp(lower(tmp(1)),'d')
			TTLSOURCE='d';
		else
			warning('Did not understand TTL source %s setting to digital in.',tmp);
			TTLSOURCE='d';
		end
	else
		warning('No TTL source given, set to digital in');
		TTLSOURCE='d';
	end
end

if ~isempty(porttoken)

	for i=1:length(port_labels)
		if ~isempty(strfind(tokens{porttoken},port_labels(i)))
			PORTS=[ PORTS port_labels(i) ];
		end
	end
	
	% if isempty assume all ports
	
	PORTS=port_labels;

end	

% parse the date

if ~isempty(datetoken)
	DATENUM=datenum([tokens{datetoken}],'yymmddHHMMSS');
end

