function fig_num=fisher_cluster(SPIKEWINDOWS,varargin)
%cluster statistics, include Fisher projection and other quality metrics
%
%
%

%
% perhaps include plot of stability across trials (maybe threshold?)
%

nparams=length(varargin);

if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs');
end

fig_num=[];

colors=[...	
	1 .6445 0;... % orange	
	1 0 1;... % magenta
	0 1 1;... % cyan	
	1 0 0;... % red				
	0 0 1;... % blue		
	.5 0 .5; ... % purple		
	.6445 .1641 .1641; ... % brown
	1 1 0;... % yellow
	.1953 .8008 .1953;... % lime-green
	.1328 .5430 .1328;... % forest green
	0 0 .5;... % navy
	0 .5 .5;... % teal
	.5430 .2695 .0742;... % saddle-brown
	1 .2695 0;... % orange red
	]; 

% just in case add the hot colormap at the end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'fig_num'
			fig_num=varargin{i+1};
		case 'colors'
			color=varargin{i+1};
		case 'maxlag'
			maxlag=varargin{i+1};
		case 'xres'
			xres=varargin{i+1};
	end
end

colorappend=colormap('hot');
colors=[colors;colorappend];

% fisher LDA for each cluster
% now compare all combinations

K=length(SPIKEWINDOWS);

clustercombos=nchoosek(1:K,2);
ncombos=size(clustercombos,1);

if isempty(fig_num)
	fig_num=figure('Visible','on','position',[0 0 150*K 150*K]);
else
	fig_num=figure('Visible','off','position',[0 0 150*K 150*K]);
end

for i=1:ncombos

	% compute the discriminant

	cluster1=SPIKEWINDOWS{clustercombos(i,1)}';
	cluster2=SPIKEWINDOWS{clustercombos(i,2)}';

	[density1,density2,xi]=fisher_projection(cluster1,cluster2);

	xcoord=[xi fliplr(xi)];
	ycoord1=[zeros(size(xi)) fliplr(density1)];
	ycoord2=[zeros(size(xi)) fliplr(density2)];

	% on a Mac transparency is completely borked :(

	subaxis(K,K,clustercombos(i,2),clustercombos(i,1),'margin',.1,'spacingvert',0,'spacinghor',0);

	patch(xcoord,ycoord1,1,'facecolor',colors(clustercombos(i,1),:),...
		'edgecolor','none');
	hold on
	patch(xcoord,ycoord2,1,'facecolor',colors(clustercombos(i,2),:),...
		'edgecolor','none');

	set(gca,'layer','top','xtick',[],'ytick',[],'linewidth',2);
	axis tight;

end


