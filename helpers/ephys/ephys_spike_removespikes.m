function SPIKELESS=ephys_spike_removespikes(DATA,SPIKES)
%
%
%
%
%


nspikes=length(SPIKES.times);
edges=round(SPIKES.frame*SPIKES.fs);

SPIKELESS=[];
leftedge=1;

for i=1:nspikes

	% left edge of next spike

	rightedge=SPIKES.times(i)-edges(1);

	if leftedge<1
		continue;
	end

	SPIKELESS=[SPIKELESS;DATA(leftedge:rightedge)];

	leftedge=SPIKES.times(i)+edges(2);

end

