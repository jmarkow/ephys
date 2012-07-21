function [BIN BINEDGES]=prctile_bin(DATA,BININTERVAL)
%
%
%
%
%

% takes a vector DATA and bins according to prctile

% e.g. if the BININTERVAL is 10, then 

% iterate through each bin until we reach the end

curr_perc=0;

BIN=zeros(size(DATA));
COUNTER=1;

while curr_perc<100

	edge_perc=curr_perc+BININTERVAL;

	if edge_perc>100, edge_perc=100; end

	if COUNTER==1
		left_edge=-inf;
	else
		left_edge=prctile(DATA,curr_perc);
	end

	if edge_perc>=100
		right_edge=inf;
	else
		right_edge=prctile(DATA,edge_perc);
	end

	BIN(find(DATA>left_edge&DATA<=right_edge))=COUNTER;

	BINEDGES(:,COUNTER)=[curr_perc edge_perc];

	curr_perc=edge_perc;
	COUNTER=COUNTER+1;

end
