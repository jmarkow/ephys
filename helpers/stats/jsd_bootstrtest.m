function [H,P,JSD,NULL]=jsd_bootstrtest(P,Q,varargin)
%
%
%
%
%
%

alpha=.01;

binmethod='root'; % root method for choosing number of bins, other
		  % options include 'friedman' or 'scott'
nbins=100; % only set if you are using the constant bin method
trials=10e3;
tail='left';

nparams=length(varargin);

if mod(nparams,2)>0
	error('ephysPipeline:argChk','Parameters must be specified as parameter/value pairs!');
end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'binmethod'
			binmethod=varargin{i+1};
		case 'nbins'
			nbins=varargin{i+1};
		case 'trials'
			trials=varargin{i+1};
		case 'tail'
			tail=varargin{i+1};
	end
end

if ~isvector(P) || ~isvector(Q)
	error('ephysPipeline:jsdbootstr:inputnotvector','Input must be vectors.');
end

n_p=length(P);
n_q=length(Q);
NULL=zeros(trials,1);

parfor i=1:trials

	% case resample the first distribution

	sample1=P(randsample(n_p,n_p,1));
	sample2=P(randsample(n_p,n_p,1));

	NULL(i)=kld(sample1,sample2,'binmethod',binmethod,'nbins',nbins,'jsd',1);

end

JSD=kld(P,Q,'binmethod',binmethod,'nbins',nbins,'jsd',1);

switch lower(tail(1))
	
	case 'l'

		P=1-(sum(JSD>=NULL)/trials);
		H=P<=alpha;

	case 'r'

		P=1-(sum(JSD<=NULL)/trials);
		H=P<=alpha;

	case 'b'

	otherwise

end
