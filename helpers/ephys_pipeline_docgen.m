function ephys_pipeline_docgen()
%
%

% run one directory up from from the primary directory

ignoredirectories{1}='external';
ignoredirectories{2}='doc';

m2html('mfiles','ephys','htmldir','ephys/doc','recursive','on','global','on',...
	'template','blue','index','menu','graph','on','ignoredDir',ignoredirectories)
