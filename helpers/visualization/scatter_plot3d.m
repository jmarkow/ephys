function FIGNUM=scatter_plot3d(DATA,LABELS,varargin)
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
shadow=0;
colors=colormap('lines');

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
		case 'shadow'
			shadow=varargin{i+1};
	end
end

if isempty(fignum)
	FIGNUM=figure();
else
	FIGNUM=fignum;
end

% for each label, get the max euclidean distance for those points, which we'll
% use to scale each sphere

uniqlabels=unique(LABELS);
r=.5;

for i=1:length(uniqlabels)
 
	grouppoints=find(LABELS==uniqlabels(i));
	[x,y,z]=sphere(40);

	for j=1:length(grouppoints)
		colors(i,:)
		mesh(x*r+DATA(grouppoints(j),1),...
			y*r+DATA(grouppoints(j),2),...
			z*r+DATA(grouppoints(j),3),...
			'FaceColor',colors(i,:),'EdgeColor','none','FaceAlpha',.5);

		hold on;
	end


end

xlimits=xlim()
ylimits=ylim()
zlimits=zlim()

[points dims]=size(DATA);

if shadow
	for i=1:points

		sh1=[xlimits(2),DATA(i,2),DATA(i,3)];
		sh2=[DATA(i,1),ylimits(2),DATA(i,3)];
		sh3=[DATA(i,1),DATA(i,2),zlimits(1)];

		scatter3(sh1(1),sh1(2),sh1(3),100,[.4 .4 .4],'filled');
		scatter3(sh2(1),sh2(2),sh2(3),100,[.4 .4 .4],'filled');
		scatter3(sh3(1),sh3(2),sh3(3),100,[.4 .4 .4],'filled');

	end
end
