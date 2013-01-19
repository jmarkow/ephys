function prettify_axis(FIGHANDLE,varargin)
%changes axis to look sensible
%
%
%
%

nparams=length(varargin);

if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs');
end

linewidth=1.5;
fontsize=11;
font='Helvetica';
ticklength=[.02 .02];

for i=1:2:nparams
	switch lower(varargin{i})
		case 'linewidth'
			linewidth=varargin{i+1};
		case 'fontsize'
			fontsize=varargin{i+1};
		case 'font'
			font=varargin{i+1};
		case 'ticklength'
			ticklength=varargin{i+1};

	end
end

% tweak the axes

all_axis=findall(FIGHANDLE,'Type','Axes');
set(all_axis,'Color','w','XColor','k','YColor','k','ticklength',ticklength,...
	'FontSize',fontsize,'FontName',font,'Linewidth',linewidth,'tickdir','out','layer','top');
