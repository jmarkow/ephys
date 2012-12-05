
function [consensus]=chirp_consensus(signal,SAMPLING,ARThreshold,angles,Nfft,dx,sono,i,ANGLE_CLASS);

	% sono is vestigial

%% matlab script zf_contours.m

%This is the single channel soundfile - keep this file small  - unless you have huge
%system memory.

%signal=wavread('finch_original.wav');
%SAMPLING=44100; %sampling rate of the signal
FREQUENCYLIMIT=12; %Upper Range of plot in kHz
%% Parameters  
%NtScale = length(TScale);
Nangle = length(angles);
 %Tested only for files shorter than one second. Memory issues may arise for longer files.

 sigma=1;
 tScale=sigma;
 sigmacount=1;
 Nfft=1024;

 for angle_variable=1:Nangle
	 [sigmacount angle_variable];

	 theta = angles(angle_variable);

	 tScalef=1.0;

	 %if(ANGLE_CLASS==1) theta=0;
	 %else
	 %   theta=pi/2;
	 %end



	 s=-1*(imag(dx* (cos(theta)+(1j* tScalef*sin(theta))))<0)+(imag(dx*(cos(theta)+(1j* (tScalef)*sin(theta))))>0);


	 %s=sign(imag(dx));


	 [gx gy]=gradient(s);

	 if(ANGLE_CLASS==1)
		 BW=gy<0; %depending on the chirplet trasform this sign should be opposite
	 else
		 BW=gx<0; %depending on the chirplet trasform this sign should be opposite
	 end



	 CC = bwconncomp(BW); %build a strucure containing a separate entry for each contour (connected component of BW).
	 cc_pix=regionprops(CC,'Area'); %this is the length of each contour
	 weightv=[];

	 for i = 1:length(cc_pix),
		 %indices=CC.PixelIdxList(i);
		 %weightv(i)=sum(sum(sono(indices{1})));
		 weightv(i)=cc_pix(i).Area;
	 end

	 %figure(1);
	 %hist(weightv,50);
	 %drawnow

	 a=find(weightv>prctile(weightv,ARThreshold));
	 %a=find(weightv>=median(weightv)+mad(weightv,1));

	 tempv=zeros(size(dx));

	 for(ik=1:length(a)),
		 ind=CC.PixelIdxList(a(ik));
		 tempv(ind{1})=1;
	 end          
	 BWallk{angle_variable}=sparse(tempv);

 end

 consensus=sparse(zeros(size(dx)));

 for angle_variable=1:Nangle,

	 consensus=consensus+BWallk{angle_variable};

 end






