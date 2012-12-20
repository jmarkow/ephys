function [CONSENSUS]=chirp_consensus(SONOGRAM,DX,ANGLES);
%
%

if nargin<2
	ANGLES=-pi/4:pi/16:pi/4;
end

if nargin<1
	error('ephysPipeline:chirpconsensus:nodata','Need the displacement to continue');
end

nangles= length(ANGLES);
Nfft=1024;
tScalef=1.0;

CONSENSUS=zeros(size(DX));
sonophase=angle(SONOGRAM);

parfor i=1:nangles

	theta = ANGLES(i);

	% contour computation

	s=-1*(imag(DX* (cos(theta)+(1j* tScalef*sin(theta))))<0)+(imag(DX*(cos(theta)+(1j* (tScalef)*sin(theta))))>0);
	[gx gy]=gradient(s);
	
	% for complex value simply add phase of sonogram using contours as mask

	contourmask=gy<0;

	% synthesize complex valued contour image

	complexcontour=contourmask.*exp(1j*(contourmask.*sonophase));

	% average over angles

	CONSENSUS=CONSENSUS+complexcontour./nangles;

end
