function [xvec,yvec]=ephys_stairplot(N,BINS)
%%%% generates a stair plot given histogram data
%
%
%
%


% define a polygon to plot as a patch

yvec=[];
xvec=[];

for i=1:length(BINS)-1
	xvec=[xvec BINS(i:i+1)];
	yvec=[yvec repmat(N(i),[1 2])];
end



