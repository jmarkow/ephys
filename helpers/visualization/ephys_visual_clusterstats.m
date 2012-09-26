function fig_num=ephys_visual_clusterstats(DATAPOINTS,CLUSTERSTATS,varargin)
%cluster statistics, include Fisher projection and other quality metrics
%
%
%

%
% perhaps include plot of stability across trials (maybe threshold?)
%

nparams=length(varargin);

if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs');
end

fs=25e3;
fig_num=[];
noise_p2p=[];
y_res=200;
spike_fs=50e3;

colors=[...	
	1 .6445 0;... % orange	
	1 0 1;... % magenta
	0 1 1;... % cyan	
	1 0 0;... % red				
	0 0 1;... % blue		
	.5 0 .5; ... % purple		
	.6445 .1641 .1641; ... % brown
	1 1 0;... % yellow
	.1953 .8008 .1953;... % lime-green
	.1328 .5430 .1328;... % forest green
	0 0 .5;... % navy
	0 .5 .5;... % teal
	.5430 .2695 .0742;... % saddle-brown
	1 .2695 0;... % orange red
	]; 


for i=1:2:nparams
	switch lower(varargin{i})
		case 'fs'
			fs=varargin{i+1};
		case 'spike_fs'
			spike_fs=varargin{i+1};
		case 'fig_num'
			fig_num=varargin{i+1};
		case 'patch_color'
			patch_color=varargin{i+1};
		case 'noise_p2p'
			noise_p2p=varargin{i+1};
		case 'y_res'
			y_res=varargin{i+1};
		case 'isi_method'
			isi_method=varargin{i+1};
	end
end

% TODO: spike autocorrelation and cross-correlations functions, simply xcorr of binned spikes-mean(lambda)

% fisher LDA for each cluster

clustermeans=CLUSTERSTATS.mu;
clustercovs=CLUSTERSTATS.Sigma;

[K,M]=size(clustermeans);

% now compare all combinations

clustercombos=nchoosek(1:K,2);
ncombos=size(clustercombos,1);


for i=1:ncombos

	% compute the discriminant

	cluster1=DATAPOINTS{clustercombos(i,1)};
	cluster2=DATAPOINTS{clustercombos(i,2)};

	mu1=mean(cluster1)';
	mu2=mean(cluster2)';

	cov1=cov(cluster1,1);
	cov2=cov(cluster2,1);

	% regularize the covariances

	cov1=cov1+eye(size(cov1)).*1e-10;
	cov2=cov2+eye(size(cov2)).*1e-10;

	w_proj=(cov1+cov2)\(mu2-mu1);

	% project both sets of points

	points1=cluster1*w_proj;
	points2=cluster2*w_proj;

	% between class separation over within class separation

	separation(i)=(w_proj'*(mu1-mu2))^2/(w_proj'*(cov1+cov2)*w_proj)

	separation_test=(abs(mean(points1)-mean(points2))^2)/(var(points1)+var(points2))

	% define mean and covariance for samples
	% plot the histogram

	xi=linspace(min([points1;points2]),max([points1;points2]),...
		sqrt(length(points1)+length(points2)));

	[density1]=ksdensity(points1,xi);
	[density2]=ksdensity(points2,xi);

	density1=density1./sum(density1);
	density2=density2./sum(density2);

	xcoord=[xi fliplr(xi)];
	ycoord1=[zeros(size(xi)) fliplr(density1)];
	ycoord2=[zeros(size(xi)) fliplr(density2)];

	% on a Mac transparency is completely borked :(

	subaxis(ncombos,1,1,i,'margin',.12,'spacingvert',0,'paddingbottom',.05);

	patch(xcoord,ycoord1,1,'facecolor',colors(clustercombos(i,1),:),...
		'edgecolor','none','linewidth',2);
	hold on
	patch(xcoord,ycoord2,1,'facecolor',colors(clustercombos(i,2),:),...
		'edgecolor','none','linewidth',2);

	set(gca,'layer','top')

end




%linkaxes(ax,'x');

