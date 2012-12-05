function [filter_sound]=filter_sound_sam(filter_sound)

% [filter_sound]=filter_sound_sam(filter_sound)
% Written by Sigal Saar

xv=zeros(1,3);
yv=zeros(1,3);

  for i=1:length(filter_sound)
        xv(1) = xv(2);
        xv(2) = xv(3);
        xv(3) = filter_sound(i) /1.020353514; % GAIN;
        yv(1) = yv(2);
        yv(2) = yv(3);
        yv(3) =   (xv(1) + xv(3)) - 2 * xv(2) + ( -0.9605029194 * yv(1)) + (  1.9597070338 * yv(2));
        filter_sound(i) = yv(3);
  end
