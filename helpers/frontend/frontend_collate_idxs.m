function [EXTRACTION_PTS]=frontend_collate_idxs(DETECTION,PAD)
% takes a vector where 1 indicates a hit and 0 nohit, and organizes
% into a series of extraction points suitable for extracting data
%
%
%

if nargin<2, PAD=1; end
if nargin<1, error('Need detection vector to continue.'); end

% find the rising and falling edges

idx=1:length(DETECTION)-1;

rising=find((DETECTION(idx)==0)&(DETECTION(idx+1)==1));
falling=find((DETECTION(idx)==1)&(DETECTION(idx+1)==0));

% now for each rising edge, find a matching falling edge

EXTRACTION_PTS=[];

if isempty(DETECTION) | (isempty(rising)&isempty(falling))
	return;
end

for i=1:length(rising)

	left_edge=rising(i);
	right_edge=min(falling(falling>rising(i)));

	% check for orphan rising edge

	if isempty(right_edge) & i==length(rising)
		right_edge=length(DETECTION);
	elseif isempty(right_edge) 
		error('Could not find a matching falling edge...');
	end

	EXTRACTION_PTS=[ EXTRACTION_PTS; left_edge right_edge ];

end

% is the first falling edge an orphan?

if ~isempty(falling)
	left_edge=max(rising(rising<falling(1)));

	if isempty(left_edge)
		left_edge=1;
	end

	right_edge=falling(1);

	EXTRACTION_PTS=[ EXTRACTION_PTS; left_edge right_edge ];
end 

% sort by left_edge

[~,idx]=sort(EXTRACTION_PTS(:,1));
EXTRACTION_PTS=EXTRACTION_PTS(idx,:);

% now account for the pads

npairs=size(EXTRACTION_PTS,1);
EXTRACTION_PTS=EXTRACTION_PTS+repmat([-PAD PAD],[npairs 1]);

EXTRACTION_PTS(EXTRACTION_PTS<1)=1;
EXTRACTION_PTS(EXTRACTION_PTS>length(DETECTION))=length(DETECTION);

% now find points where right_edge overlaps with left_edge and merge

exit_flag=0;

while exit_flag==0

	exit_flag=1;

	for i=1:npairs-1

		right_edge=EXTRACTION_PTS(i,2);
		left_edge=EXTRACTION_PTS(i+1,1);

		if right_edge>left_edge
			
			% merge operation
			
			EXTRACTION_PTS(i,2)=EXTRACTION_PTS(i+1,2);
			EXTRACTION_PTS(i+1,:)=[];

			npairs=npairs-1;
			exit_flag=0;
			break;
		end
	end
end
