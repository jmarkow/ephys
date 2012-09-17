function [C]=discriminant_analysis(TESTDATA,TRAINDATA,CLASS,varargin)
%Simple wrapper for MATLAB's classify function, only supports two classes
%at the moment
%
%
%
%

% if the testdata is empty, create a mesh spanning the x and y range

nparams=length(varargin);
spike_fs=100e3;
type='quadratic';
xres=600;
yres=600;
colors=[ 0 1 1 ; 1 .6445 0 ];
x_label='Norm. waveform score';
y_label='ISI score';
figtitle='QDA';

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

if isempty(TESTDATA)
	[X Y]=meshgrid(linspace(xrange(1),xrange(2),xres),...
		linspace(yrange(1),yrange(2),yres));
	TESTDATA=[X(:) Y(:)];
end

if size(TESTDATA,2)~=size(TRAINDATA,2)
	error('ephysPipeline:discriminantanalysis:unequaldims',...
		'Testing and training data do not have the same number of columns');
end


[C,err,P,logp,coeff]=classify(TRAINDATA(CLASS==2,:),TRAINDATA,CLASS,type);

threshold=quantile(P(:,1),.95);
n=(2*threshold)/(1+2*threshold);
p=1/(1+2*threshold);
priors=[ n p ];

[C,err,P,logp,coeff]=classify(TESTDATA,TRAINDATA,CLASS,type,priors);

err

% finally classify the testing data

% print the error rate (confusion matrix)
% specify the classification boundary function

switch lower(type)

	case 'linear'

		k=coeff(1,2).const;
		l=coeff(1,2).linear;

		f=@(x,y) k + l(1)*x + l(2)*y;

	case 'quadratic'

		k=coeff(1,2).const;
		l=coeff(1,2).linear;
		q=coeff(1,2).quadratic;

		f=@(x,y) k + l(1)*x + l(2)*y + q(1,1)*x^2 + (q(1,2)+q(2,1))*x.*y + q(2,2)*y.^2;

	otherwise
		
		error('ephysPipeline:discriminantanalysis','Did not understand classification type');

end


% figure suitable for publication

figure();
points=find(CLASS==1);
scatter(TRAINDATA(points,1),TRAINDATA(points,2),100,[ 1 0 1 ],'x');
hold on

points=find(CLASS==2);
scatter(TRAINDATA(points,1),TRAINDATA(points,2),100,[ 0 1 0 ],'filled');

boundaryline=ezplot(f,[ xrange yrange ]);
set(boundaryline,'Color',[.4 .4 .4],'linewidth',2);

axis tight

xlabel(x_label);
ylabel(y_label);
title('')
prettify_axis(gca,'FontSize',12,'TickLength',[.025 .025],'linewidth',2,'FontName','Helvetica');
prettify_axislabels(gca,'FontSize',15,'FontName','Helvetica');

figure();
points=find(C==1);
scatter(TESTDATA(points,1),TESTDATA(points,2),100,[ 1 0 1 ],'x');
hold on

points=find(C==2);
scatter(TESTDATA(points,1),TESTDATA(points,2),100,[ 0 1 0 ],'filled');

boundaryline=ezplot(f,[ xrange yrange ]);
set(boundaryline,'Color',[.4 .4 .4],'linewidth',2);

axis tight

xlabel(x_label);
ylabel(y_label);
title('')
prettify_axis(gca,'FontSize',12,'TickLength',[.025 .025],'linewidth',2,'FontName','Helvetica');
prettify_axislabels(gca,'FontSize',15,'FontName','Helvetica');


