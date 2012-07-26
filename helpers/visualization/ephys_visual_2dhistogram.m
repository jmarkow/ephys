function [DENSITY,TIME,V]=ephys_visual_2dhistogram(SPIKEWINDOWS,varargin)
%generates a 2D histogram for spike visualization
%
%
%
%
%

nparams=length(varargin);

if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs');
end


SR=50e3;
fig_num=[];
patch_color=[1 .6 0];
y_res=200;
noise_p2p=[];

for i=1:2:nparams
	switch lower(varargin{i})
		case 'sr'
			SR=varargin{i+1};
		case 'fig_num'
			fig_num=varargin{i+1};
		case 'noise_p2p'
			noise_p2p=varargin{i+1};
		case 'y_res'
			y_res=varargin{i+1};
	end
end

[samples,trials]=size(SPIKEWINDOWS);

TIME=([1:samples]./SR)*1e3;
timevec_mat=[1:samples];
coordmat=[];

for i=1:trials

	currwin=SPIKEWINDOWS(:,i);

	% form matrix with spike x y pairings

	coordmat_tmp=[timevec_mat(:) currwin(:)];
	coordmat=[coordmat;coordmat_tmp];

end

% construct 2D histogram
% for voltmin and voltmax cover 80% of voltage values

voltmin=prctile(coordmat(:,2),1)-10;
voltmax=prctile(coordmat(:,2),99)+10;

edges{1}=.5:1:samples+.5;
edges{2}=linspace(voltmin,voltmax,y_res);
DENSITY=hist3(coordmat,'Edges',edges);

V=edges{2};


