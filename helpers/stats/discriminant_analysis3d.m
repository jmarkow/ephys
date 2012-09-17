function discriminant_analysis3d(TRAINDATA,CLASS,varargin)
%Simple wrapper for MATLAB's classify function, only supports two classes
%at the moment
%
%
%
%

% if the testdata is empty, create a mesh spanning the x and y range
nparams=length(varargin);
spike_fs=100e3;
type='linear';
xres=100;
yres=100;
zres=100;
colors=[ 0 1 1 ; 1 .6445 0 ];

if mod(nparams,2)>0
	error('ephysPipeline:argChk','Parameters must be specified as parameter/value pairs!');
end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'type'
			type=varargin{i+1};
		case 'labels'
			labels=varargin{i+1};
		case 'colors'
			colors=varargin{i+1};
	end
end

xrange=[ min(TRAINDATA(:,1)) max(TRAINDATA(:,1)) ];
yrange=[ min(TRAINDATA(:,2)) max(TRAINDATA(:,2)) ];
zrange=[ min(TRAINDATA(:,3)) max(TRAINDATA(:,3)) ];

%if size(TESTDATA,2)~=size(TRAINDATA,2)
%	error('ephysPipeline:discriminantanalysis:unequaldims',...
%		'Testing and training data do not have the same number of columns');
%end

[C,err,P,logp,coeff]=classify(TRAINDATA,TRAINDATA,CLASS,type);

% specify the classification boundary function

[X,Y,Z]=meshgrid(linspace(xrange(1),xrange(2),xres),linspace(yrange(1),yrange(2),yres),...
	linspace(zrange(1),zrange(2),zres));

switch lower(type)

	case 'linear'

		k=coeff(1,2).const;
		l=coeff(1,2).linear;

		f=@(x,y,z) k + [x y z]*l;

	case 'quadratic'

		k=coeff(1,2).const;
		l=coeff(1,2).linear;
		q=coeff(1,2).quadratic;

		f=@(x,y,z) k + [x y z]*l + sum([x y z].*([x y z]*q),2);

	otherwise
		
		error('ephysPipeline:discriminantanalysis','Did not understand classification type');

end


v=f(X(:),Y(:),Z(:));
v=reshape(v,size(X));

classes=unique(CLASS);

figure();

colors={'g','b','r'};
for i=1:length(classes)
	points=find(CLASS==classes(i));
	plot3(TRAINDATA(points,1),TRAINDATA(points,2),TRAINDATA(points,3),'.','color',colors{i},...
		'markersize',25);
	hold on;
end
p=patch(isosurface(X,Y,Z,v,-3),'FaceColor','r','edgecolor','none');
camlight;
%gscatter(TRAINDATA(:,1),TRAINDATA(:,2),CLASS,'.',20,'off');




