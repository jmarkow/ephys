function [H,P,JSD,NULL]=jsd_bootstrtest(P,Q,varargin)
%
%
%
%
%
%

alpha=.01;

binmethod='c'; % root method for choosing number of bins, other
		  % options include 'friedman' or 'scott'
nbins=75; % only set if you are using the constant bin method
trials=3e3;
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

JSD=kld(P,Q,'binmethod',binmethod,'nbins',nbins,'jsd',0);

n_p=length(P);
n_q=length(Q);
NULL=zeros(trials,1);

samplepool=[P(:);Q(:)];
npool=length(samplepool);

for i=1:trials

	% case resample the first distribution

	sample1=samplepool(randsample(npool,n_p,1));
	sample2=samplepool(randsample(npool,n_q,1));
	NULL(i)=kld(sample1,sample2,'binmethod',binmethod,'nbins',nbins,'jsd',0);


end

switch lower(tail(1))
	
	case 'l'
		
		P=sum(JSD>=NULL)/trials;
		H=P<=alpha;

	case 'r'
		
		P=sum(JSD<=NULL)/trials;
		H=P<=alpha;

	case 'b'

	otherwise

end
