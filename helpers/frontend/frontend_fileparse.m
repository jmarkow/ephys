function [BIRDID,RECID,MICTRACE,TTLTRACE,DATENUM]=frontend_fileparse(FILENAME,DELIM,FMT,DATEFMT)
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
if nargin<3 | isempty(FMT), FMT='bimtdd'; end
if nargin<2 | isempty(DELIM), DELIM='\_'; end
if nargin<1 | isempty(FILENAME), error('Need filename to continue!'); end

BIRDID=[];
RECID=[];
MICTRACE=[];
TTLTRACE=[];
DATENUM=[];

birdtoken=strfind(lower(FMT),'b');
rectoken=strfind(lower(FMT),'i');
mictoken=strfind(lower(FMT),'m');
ttltoken=strfind(lower(FMT),'t');
datetoken=strfind(lower(FMT),'d');

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
	mictokens=regexpi(tokens{mictoken},'\d+','match');
	
	if isempty(mictokens)
		error('Could not find mic trace at token %g for file %s', mictoken,FILENAME);
	end

	MICTRACE=str2num(mictokens{1});
end


if ~isempty(ttltoken)
	ttltokens=regexpi(tokens{ttltoken},'\d+','match');

	if isempty(ttltokens)
		error('Could not find ttl trace at token %g for file %s', ttltoken,FILENAME);
	end

	TTLTRACE=str2num(ttltokens{1});
end

% parse the date

DATENUM=datenum([tokens{datetoken}],'yymmddHHMMSS');

