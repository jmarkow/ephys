function [COEFFS]=get_geoemetric_coefficients(WAVEFORMS)
%
%
%
%
%
%

[samples,trials]=size(WAVEFORMS);

% take positive and negative peak

pospeak=max(WAVEFORMS);
negpeak=min(WAVEFORMS);

% positive and negative energy

% total energy

totalenergy=sum(WAVEFORMS.^2);

% positive and negative gradients

idx=1:samples;
neo=zeros(trials,1);
posgrad=zeros(trials,1);
neggrad=zeros(trials,1);
width=zeros(trials,1);
posenergy=zeros(trials,1);
negenergy=zeros(trials,1);

for i=1:trials

	% neo operator

	currwave=WAVEFORMS(:,i);

	[val loc]=min(currwave);
	negenergy(i)=sum(currwave(currwave<0).^2);

	if loc>1 && loc<samples-1
		neo(i)=currwave(loc)^2-currwave(loc-1)*currwave(loc+1);
	else
		neo(i)=NaN;
	end

	for j=loc:-1:1
		if currwave(j)>.25*val
			break;
		end
	end

	cross1=j;
	neggrad(i)=mean(diff(currwave(cross1:loc)));

	% get the width

	for j=loc:1:length(currwave)
		if currwave(j)>.25*val
			break;
		end
	end

	cross2=j;
	width(i)=abs(cross1-cross2);
	posgrad(i)=mean(diff(currwave(loc:cross2)));

	% get the pos gradient

	posenergy(i)=sum(currwave(currwave>0).^2);

end

% width seems unreliable, along with pos and neg grad

COEFFS=[pospeak' negpeak' log(posenergy) log(negenergy) log(totalenergy') neo width posgrad neggrad];


