function FIGNUM=pretty_polar(DATA,LABELS,varargin)
%wrapper for rose plot for circular histogramming
%
%
%
%

[samples,dims]=size(DATA);

if dims>3
	warning('ephysPipeline:sphereplot:toomanydims','Dims>3, truncating to 3');
	DATA=DATA(:,1:3);
end

surfcolor=colormap('hsv');
x_label='X';
y_label='Y';
z_label='Z';

nparams=length(varargin);
fignum=[];

if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs');
end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'x_label'
			x_label=varargin{i+1};
		case 'y_label'
			y_label=varargin{i+1};
		case 'z_label'
			z_label=varargin{i+1};
		case 'fignum'
			fignum=varargin{i+1};
		case 'surfcolor'
			surfcolor=varargin{i+1};
	end
end

if isempty(fignum)
	FIGNUM=figure();
else
	FIGNUM=fignum;
end

% for each label, get the max euclidean distance for those points, which we'll
% use to scale each sphere

uniqlabels=unique(LABELS)

for i=1:length(uniqlabels)
 
	grouppoints=find(LABELS==uniqlabels(i));
	groupmat=DATA(grouppoints,:)

	npoints=length(grouppoints);

	% simple mean should be fine...

	groupmean=mean(groupmat,1)
	
	%groupdist=squareform(pdist(grouppoints,'euclidean'));

	if length(grouppoints)<2
		continue;
	end

	dist=zeros(1,npoints);

	for j=1:npoints
		dist(j)=sqrt(sum((groupmat(j,:)-groupmean).^2));
	end

	r=max(dist)/2;

	% get sphere coordinates

	[x,y,z]=sphere;

	surf(x*r+groupmean(1),y*r+groupmean(2),z*r+groupmean(3));
	hold on;

	% could use simply mean and scale by max pairwise distance,
	% or assume that the first point is the center (i.e. how far
	% do we get from the first point?)

end
