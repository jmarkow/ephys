function ISI=get_spike_correlogram(SPIKES1,SPIKES2,varargin)
%Computes the ISI vector from clust_spike_vec
%
%

% get the total number of spikes

nparams=length(varargin);

if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs');
end

fs=25e3;
fig_num=[];
xres=.001;
maxlag=.2;
type='auto';

color=[0 0 0];

% just in case add the hot colormap at the end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'fs'
			fs=varargin{i+1};
		case 'fig_num'
			fig_num=varargin{i+1};
		case 'color'
			color=varargin{i+1};
		case 'xres'
			xres=varargin{i+1};
		case 'maxlag'
			maxlag=varargin{i+1};
		case 'type'
			type=varargin{i+1};
	end
end


if length(SPIKES1)~=length(SPIKES2)
	error('ephysPipeline:spikecorrel:unequaltrials','Unequal number of trials');
end

if iscell(SPIKES1)
    ntrials=length(SPIKES1);
else
    ntrials=1;
end

% center on each spike, check up to maxlag using nbins

windowbins=-maxlag:xres:maxlag;

% should have time in secs here

density=zeros(1,length(windowbins));
counter=0;
for i=1:ntrials
	
    if ntrials>1
        spiketimes1=SPIKES1{i};
        spiketimes2=SPIKES2{i};
    else
        spiketimes1=SPIKES1;
        spiketimes2=SPIKES2;
    end
        
	% center on each spike, subtract time from second 
	% vector of spiketimes and bin

	for j=1:length(spiketimes1)
		
		spiketmp=spiketimes2;

		if lower(type(1))=='a'
			spiketmp(j)=[];
		end

		spikediff=spiketmp-spiketimes1(j);

		density=density+histc(spikediff,windowbins);
		counter=counter+1;
	end
end

if isempty(fig_num)
	fig_num=figure();
end

h=bar(windowbins*1e3,(density./counter)*(1/xres),'histc');
set(h,'EdgeColor','none','FaceColor',color);

