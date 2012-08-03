function IFRVEC=ephys_ifr(SPIKETIMES,TIMEVECTOR,FS)
%computes IFR per the Hahnloser et al. 2001 definition
%

%
%

% provide spikes in samples
% supply the length of the full vector to estimate over

SPIKETIMES=SPIKETIMES(:);

% shouldn't be any spike samples before 1, also add end of timevector

if length(SPIKETIMES==1)<1
	SPIKETIMES=[1;SPIKETIMES];
end

if length(SPIKETIMES==length(TIMEVECTOR))<1
	SPIKETIMES=[SPIKETIMES;length(TIMEVECTOR)];
end

IFRVEC=zeros(TIMEVECTOR,1);

for i=1:length(SPIKETIMES)-1
	IFRVEC(SPIKETIMES(i)+1:SPIKETIMES(i+1))=1/((SPIKETIMES(i+1)-SPIKETIMES(i))/FS);
end
