function [N,BINVEC]=pretty_histogram(DATA,BINS)
%
%
%
%
%
%

[n,x]=histc(DATA,BINS);

N=[];
BINVEC=[];

for j=1:length(BINS)-1
	BINVEC=[BINVEC BINS(j:j+1)];
	N=[N;repmat(n(j),2,1)];
end

