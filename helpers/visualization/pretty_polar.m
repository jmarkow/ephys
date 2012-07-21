function FIGNUM=pretty_polar(ANGLES,NBINS,varargin)
%wrapper for rose plot for circular histogramming
%
%
%
%

filled=1;
linewidth=1;
labels={'0','$\pi$/2','$\pm\pi$','-$\pi$/2'};
nparams=length(varargin);
fig_title=[];
x_label=[];
fignum=[];

if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs');
end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'filled'
			filled=varargin{i+1};
		case 'linewidth'
			linewidth=varargin{i+1};
		case 'labels'
			labels=varargin{i+1};
		case 'fig_title'
			fig_title=varargin{i+1};
		case 'x_label'
			x_label=varargin{i+1};
		case 'y_label'
			y_label=varargin{i+1};
		case 'fignum'
			fignum=varargin{i+1};
	end
end

if isempty(fignum)
	FIGNUM=figure();
else
	FIGNUM=fignum;
end

[n,x]=rose(ANGLES,NBINS);
p=polar(n,x);

if filled
	x=get(p,'XData');
	y=get(p,'YData');
	g=patch(x,y,'y','FaceColor',[1 0 0],'EdgeColor',[.2 .2 .2],'linewidth',linewidth);
else
	set(p,'color',facecolor,'linewidth',linewidth);
end

h=findall(gcf,'type','line');
h(h==p)=[];
delete(h);

hidden_text=findall(gca,'type','text');

angles=0:30:330;
k=0;
nodel=[];
todel=[];
for ang=angles
	hObj=findall(hidden_text,'string',num2str(ang));

	switch ang
		case 0
			nodel=[nodel hObj];
			set(hObj,'string',labels{1},'HorizontalAlignment','Left','FontSize',12,'FontName','Helvetica','Interpreter','Latex');
		case 90
			nodel=[nodel hObj];
			set(hObj,'string',labels{2},'VerticalAlignment','Bottom','FontSize',12,'FontName','Helvetica','Interpreter','Latex');
		case 180
			nodel=[nodel hObj];
			set(hObj,'string',labels{3},'HorizontalAlignment','Right','FontSize',12,'FontName','Helvetica','Interpreter','Latex');
		case 270 
			nodel=[nodel hObj];
			set(hObj,'string',labels{4},'VerticalAlignment','Top','FontSize',12,'FontName','Helvetica','Interpreter','Latex');
		otherwise
			todel=[todel hObj];
	end

end

% add xlabel and title to the NO DELETE list

delete(todel);
h=findall(gca,'type','text');

x_label_id=[];
title_id=[];

if ~isempty(x_label)
	xlabel(x_label,'FontName','Helvetica','FontSize',13);
	x_label_id=findall(h,'string',x_label);
end

if ~isempty(fig_title)
	title(fig_title,'FontName','Helvetica','FontSize',13);
	title_id=findall(h,'string',fig_title);
	title_id=setdiff(title_id,x_label_id); % in case we have overlapping title and xlabel
end

nodel=[nodel x_label_id title_id];

for i=1:length(h)
	if ~any(h(i)==nodel)
		delete(h(i))
	end
end


% maybe also delete other labels

