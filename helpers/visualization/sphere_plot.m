function FIGNUM=sphere_plot(DATA,LABELS,varargin)
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

surfcolor=colormap('lines');
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
	dist=zeros(npoints,3);

	% take max distance in aech direction

	for j=1:npoints
		dist(j,:)=abs(groupmat(j,:)-groupmean);
	end

	r=max(dist);

	% get sphere coordinates

	[x,y,z]=sphere(15);

	h=surfc(x*r(1)+groupmean(1),y*r(2)+groupmean(2),z*r(3)+groupmean(3));
	set(h,'FaceColor',surfcolor(i,:),'FaceAlpha',.95,'EdgeColor',[0 0 0]);
	hold on;

	% could use simply mean and scale by max pairwise distance,
	% or assume that the first point is the center (i.e. how far
	% do we get from the first point?)

end

xlimits=xlim();
ylimits=ylim();
zlimits=zlim();

for i=1:length(uniqlabels)

	grouppoints=find(LABELS==uniqlabels(i));
	groupmat=DATA(grouppoints,:)

	npoints=length(grouppoints);

	% simple mean should be fine...

	groupmean=mean(groupmat,1)
	
	%groupdist=squareform(pdist(grouppoints,'euclidean'));

	dist=zeros(npoints,3);

	% take max distance in aech direction

	for j=1:npoints
		dist(j,:)=abs(groupmat(j,:)-groupmean);
	end

	r=max(dist);

	% get sphere coordinates

	[x,y,z]=sphere(15);

	h=surf(x*r(1)+groupmean(1),y*r(2)+groupmean(2),ones(size(z)).*zlimits(1));
	set(h,'FaceColor',[.4 .4 .4],'EdgeColor','none','FaceAlpha',.15);

	h=surf(ones(size(x)).*xlimits(2),y*r(2)+groupmean(2),z*r(3)+groupmean(3));
	set(h,'FaceColor',[.4 .4 .4],'EdgeColor','none','FaceAlpha',.15);

	h=surf(x*r(1)+groupmean(1),ones(size(y)).*ylimits(2),z*r(3)+groupmean(3));
	set(h,'FaceColor',[.4 .4 .4],'EdgeColor','none','FaceAlpha',.15);

end


view(-52.5,40);
box on
set(gca,'linewidth',2,'FontSize',15,'FontName','Helvetica');

xlabel('PC1','FontSize',15,'FontName','Helvetica');
ylabel('PC2','FontSize',15,'FontName','Helvetica');
zlabel('PC3','FontSize',15,'FontName','Helvetica');
