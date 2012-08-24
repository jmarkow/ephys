function DIV=kld(P,Q,varargin)
%computes the Kullback-Leibler divergence for two populations P and Q
%, also includes the ability to compute Jensen-Shannon divergence
%
%


binmethod='root'; % root method for choosing number of bins, other
		  % options include 'friedman' or 'scott'
jsd=1; % setting jsd to 1 will compute Jensen Shannon divergence
       % instead of kld

nbins=100; % only set if you are using the constant bin method

nparams=length(varargin);

if mod(nparams,2)>0
	error('ephysPipeline:argChk','Parameters must be specified as parameter/value pairs!');
end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'jsd'
			jsd=varargin{i+1};
		case 'binmethod'
			binmethod=varargin{i+1};
		case 'nbins'
			nbins=varargin{i+1};
	end
end

if ~isvector(P) || ~isvector(Q)
	error('ephysPipeline:kld:notvectorinputs','Need to supply two vectors as P and Q');
end

P=P(:);
Q=Q(:);

n_p=length(P);
n_q=length(Q);

n_max=max([n_p n_q]);

switch lower(binmethod(1))

	case 'r'

		bins_1=sqrt(n_p);
		bins_2=sqrt(n_q);

	case 'f'

		bins_1=2*(iqr(P)/(n_p^(1/3)));
		bins_2=2*(iqr(Q)/(n_q^(1/3)));

	case 's'



	case 'c'

		bins_1=nbins;
		bins_2=nbins;
		
	otherwise
		
		error('ephysPipeline:kld:binmethod','Did not understand bin method %s',binmethod);

end

bins=round(mean([bins_1 bins_2]));

binedges=linspace(1,max([P;Q]),bins);

[density_1]=histc(P,binedges);
[density_2]=histc(Q,binedges);

pi=density_1./sum(density_1);
qi=density_2./sum(density_2);

pi=pi+eps;
qi=qi+eps;

if ~jsd
	DIV=sum(pi.*(log2(pi)-log2(qi)));
else
	m=(pi+qi)./2;

	d_pm=sum(pi.*(log2(pi)-log2(m)));
	d_qm=sum(qi.*(log2(qi)-log2(m)));

	DIV=.5*d_pm+.5*d_qm;
end
