function ISI=get_isi(CLUSTERSPIKES)
%Computes the ISI vector from clust_spike_vec
%
%

% get the total number of spikes

nspikes=sum(cellfun(@length,CLUSTERSPIKES));
ntrials=length(CLUSTERSPIKES);

ISI=zeros(nspikes-ntrials,1);
counter=1;

for i=1:ntrials
	tmpisi=diff(CLUSTERSPIKES{i});
	ISI(counter:counter+length(tmpisi)-1)=tmpisi;

	counter=counter+length(tmpisi);
end
