function [C,FP,FN]=discriminant_analysis(TESTDATA,TRAINDATA,CLASS,varargin)
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
		case 'x_label'
			x_label=varargin{i+1};
		case 'y_label'
			y_label=varargin{i+1};
	end
end

xrange=[ min([TRAINDATA(:,1);TESTDATA(:,1)]) max([TRAINDATA(:,1);TESTDATA(:,1)]) ];
yrange=[ min([TRAINDATA(:,2);TESTDATA(:,2)]) max([TRAINDATA(:,2);TESTDATA(:,2)]) ];

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

threshold=quantile(P(:,1),.73);
n=(2*threshold)/(1+2*threshold);
p=1/(1+2*threshold);
priors=[ n p ];
%priors=[ .5 .5 ];

[C,err,P,logp,coeff]=classify(TRAINDATA,TRAINDATA,CLASS,type,priors);

% false positives

FP=sum(C==2&CLASS==1);

% false negatives

FN=sum(C==1&CLASS==2);

FP+FN

[C,err,P,logp,coeff]=classify(TESTDATA,TRAINDATA,CLASS,type,priors);

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
subplot(1,2,1);
points=find(CLASS==1);
scatter(TRAINDATA(points,1),TRAINDATA(points,2),100,[ .5 .5 .5 ],'x');
hold on

points=find(CLASS==2);
scatter(TRAINDATA(points,1),TRAINDATA(points,2),60,[ 0 0 0 ]);

boundaryline=ezplot(f,[ xrange yrange ]);
set(boundaryline,'Color','k','linewidth',2);

axis([ xrange yrange ]);

xlabel(x_label);
ylabel(y_label);
title('Training labels')
prettify_axis(gca,'FontSize',15,'TickLength',[.025 .025],'linewidth',2,'FontName','Helvetica');
prettify_axislabels(gca,'FontSize',20,'FontName','Helvetica');

subplot(1,2,2);
points=find(C==1);
scatter(TESTDATA(points,1),TESTDATA(points,2),100,[ .5 .5 .5 ],'x');
hold on

points=find(C==2);
scatter(TESTDATA(points,1),TESTDATA(points,2),60,[ 0 0 0 ]);

boundaryline=ezplot(f,[ xrange yrange ]);
set(boundaryline,'Color','k','linewidth',2);


axis([ xrange yrange ]);
xlabel(x_label);
ylabel(y_label);
title('Testing labels')
prettify_axis(gca,'FontSize',15,'TickLength',[.025 .025],'linewidth',2,'FontName','Helvetica');
prettify_axislabels(gca,'FontSize',20,'FontName','Helvetica');


