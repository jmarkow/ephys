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
histmax=[];
histmin=[];
bins=[];

nparams=length(varargin);

if mod(nparams,2)>0
	error('ephysPipeline:argChk','Parameters must be specified as parameter/value pairs!');
end

for i=1:2:nparams
	switch lower(varargin{i})
        	case 'bins'
            		bins=varargin{i+1};
		case 'jsd'
			jsd=varargin{i+1};
		case 'binmethod'
			binmethod=varargin{i+1};
		case 'nbins'
			nbins=varargin{i+1};
		case 'histmax'
			histmax=varargin{i+1};
		case 'histmin'
			histmin=varargin{i+1};
	end
end

if ~isvector(P) || ~isvector(Q)
	error('ephysPipeline:kld:notvectorinputs','Need to supply two vectors as P and Q');
end

P=P(:);
Q=Q(:);

n_p=length(P);
n_q=length(Q);

if ~isempty(histmax)
	n_p=sum(P<=histmax);
	n_q=sum(Q<=histmax);
end

n_max=max([n_p n_q]);

if isempty(bins)
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
	
	if isempty(histmax)
		histmax=max([P;Q]);
	end

	if isempty(histmin)
		histmin=min([P;Q]);
	end

	binedges=linspace(histmin,histmax,bins);

else

	binedges=bins;
end


[density_1]=histc(P,binedges);
[density_2]=histc(Q,binedges);

%density_1=density_1+ones(size(density_1));
%density_2=density_2+ones(size(density_2));

p_i=density_1./sum(density_1);
q_i=density_2./sum(density_2);

if ~jsd
	DIV=sum(p_i.*(log(p_i./q_i)));
else
	m=(p_i+q_i).*.5;

	d_pm=p_i.*(log(p_i./m));
	d_qm=q_i.*(log(q_i./m));

	d_pm(isnan(d_pm))=0;
	d_qm(isnan(d_qm))=0;

	d_pm=sum(d_pm);
	d_qm=sum(d_qm);

	DIV=.5*d_pm+.5*d_qm;

	%DIV=sum(abs(cumsum(p_i)-cumsum(q_i)));


end
