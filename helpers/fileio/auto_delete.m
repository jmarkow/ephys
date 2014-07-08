function auto_delete(DIR,DATE,FILT)
%auto_deletes files older than DATE from DIR
%
%	auto_delete(DIR,DATE,FILT)
%
%	DIR
%
%	DATE
%
%	FILT
%
%


% get the DIR listing

if ~isempty(FILT)
	listing=dir(fullfile(DIR,['*.' FILT]));
else
	listing=dir(DIR);
end

dates={listing.date};
filenames={listing.name};

% make sure recycle is set to off

oldstate=recycle;
newstate=recycle('off');



for i=1:length(dates)

	% for each date check to see if it is > than
	% the threshold

	% days elapsed

	delapsed=daysdif(datenum(dates{i}),datenum(now));
	
	if delapsed>DATE
		disp(['Deleting ' fullfile(DIR,filenames{i})]);
		delete(fullfile(DIR,filenames{i}));
	end

end

newstate=recycle(oldstate);

