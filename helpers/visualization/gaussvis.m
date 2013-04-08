function fig_num=gaussvis(MODEL,DATA,varargin)
%
%
%
%
%


% plot the first two dimensions

nparams=length(varargin);

fig_num=[];
colorsmap=colormap('copper');
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

if mod(nparams,2)>0
	error('ephysPipeline:argChk','Parameters must be specified as parameter/value pairs!');
end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'fig_num'
			fig_num=varargin{i+1};
	end
end

if ~isempty(fig_num)
	set(0,'CurrentFigure',fig_num);
	cla;
else
	fig_num=figure();
end

NCLUST=size(MODEL.mu,1);
NDIM=size(MODEL.mu,2);

% get the dimensions we're plotting

if NDIM<2
	error('gaussvis:toofewdimensions','Need more than 1 dimensions for plotting');
end

dimpairs=nchoosek(1:NDIM,2)

rows=max(dimpairs(:,1));
columns=max(dimpairs(:,2));

for ii=1:size(dimpairs,1)

	cx=dimpairs(ii,1);
	cy=dimpairs(ii,2);

	if size(dimpairs,1)>1
		subaxis(rows,columns,cy,cx,1,1,...
			'margin',.1,'spacingvert',.05,'spacinghor',.05);
	end

	if iscell(DATA)
		for i=1:length(DATA)

			currpoints=DATA{i};
			scatter(DATA{i}(:,cx),DATA{i}(:,cy),20,colors(i,:));
			axis tight;
			hold on;
		end
	else
		scatter(DATA(:,cx),DATA(:,cy),20);
		axis tight;
		hold on;
	end


	xlimits=xlim();
	ylimits=ylim();

	for i=1:NCLUST

		% evaluate x and y over a grid

		[x1,y1]=meshgrid(xlimits(1):.3:xlimits(2),ylimits(1):.3:ylimits(2));

		% get the mahal distance between each point in the grid and our cluster

		invcovm=inv(MODEL.sigma([cx cy],[cx cy],i));

		gridpoints=[x1(:) y1(:)];

		mahal=[];
		for j=1:size(gridpoints,1)
			mahal(j)=(gridpoints(j,:)-MODEL.mu(i,[cx cy]))*invcovm*(gridpoints(j,:)-MODEL.mu(i,[cx cy]))';
		end

		mahal=reshape(mahal,size(x1));

		% levels

		levels=chi2inv([.6827 .9545],2);

		contour(x1,y1,mahal,levels,'linewidth',1.5);
		colormap(colorsmap);
		hold on;
	end

	set(gca,'linewidth',2,'FontSize',18,'FontName','Helvetica','TickDir','out');
end
