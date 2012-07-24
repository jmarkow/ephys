function OUTPUT=im_reformat(IM,ROWS)
%
% OUTPUT=im_reformat(IM,ROWS) 
%
% Takes a matrix of intensity values and reformats so that 
% the image spans across ROWS
%
% IM
% matrix of intensity values
%
% ROWS
% number of rows for new image
%

% collect the size of the original image

[m,n]=size(IM);

% calculate the width of the new image

width=ceil(n/ROWS);

% preallocate the new image

OUTPUT=zeros(m*ROWS,width);

row_counter=0;
width_counter=0;
bar_size=5;
% place a segment of the original image in each row

for i=1:ROWS

	if width_counter+width>n
		width=n-width_counter;
	end

	OUTPUT(row_counter+1:row_counter+m,1:width)=IM(:,width_counter+1:width_counter+width);
	OUTPUT(row_counter+m-(bar_size-1):row_counter+m,1:width)=zeros(bar_size,width);
	row_counter=row_counter+m;
	width_counter=width_counter+width;

end



