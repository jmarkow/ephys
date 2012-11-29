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

	% right gradient

	maskwave=currwave;

	% mask everthing except for the slope to the min point
	% start at the minimum peak and walk back until we cross zero

	for j=loc:-1:1
		if currwave(j)>0
			break;
		end
    end
    
    cross1=j;
	maskwave(idx<=j|idx>=loc)=0;
	midpoint=min(find(maskwave<=.5*val));
    
    % get the angle at the midpoint

	gradangle=[0;atand(diff(maskwave))];
	
	if isempty(midpoint)
		neggrad(i)=NaN;
	else
		neggrad(i)=gradangle(midpoint);
    end
    
    % get the width
    
    for j=loc:1:length(currwave)
        if currwave(j)>0
            break;
        end
    end
    
    cross2=j;
    width(i)=abs(cross1-cross2);

	% get the pos gradient

	[val loc]=max(currwave);
    
    posenergy(i)=sum(currwave(currwave>0).^2);

	maskwave=currwave;

	for j=loc:-1:1
		if currwave(j)<.4*val
			break;
		end
	end

	%width(i)=abs(j-cross1);
	maskwave(idx<=j|idx>=loc)=0;
	midpoint=min(find(maskwave>=.5*val));

	gradangle=[0;atand(diff(maskwave))];
	
	if isempty(midpoint)
		posgrad(i)=NaN;
	else
		posgrad(i)=gradangle(midpoint);
    end
    
end

% width seems unreliable, along with pos and neg grad

COEFFS=[pospeak' negpeak' log(posenergy) log(negenergy) log(totalenergy') neo width posgrad neggrad];


