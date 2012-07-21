function prettify_axislabels(FIGHANDLE,varargin)
%makes axis labels sensible
%
%
%
%
nparams=length(varargin);

if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs');
end

fontsize=11;
font='Helvetica';

for i=1:2:nparams
	switch lower(varargin{i})
		case 'fontsize'
			fontsize=varargin{i+1};
		case 'font'
			font=varargin{i+1};

	end
end

% tweak the axes

all_text=findall(FIGHANDLE,'Type','Text');
set(all_text,'FontSize',fontsize,'FontName',font);

