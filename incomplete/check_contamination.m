function CONTAMINATION=check_contamination(SPIKEDATA,MODEL)
%
%
%
%

nclust=size(MODEL.mu,1);
npoints=1e4; % point from each cluster
CONTAMINATION=ones(1,nclust).*NaN;

% or generate points from each cluster
%

if any(MODEL.garbage)
	garbage=1;
	nclust=nclust+1;
	alldata=cat(1,SPIKEDATA{:});
	datarange=range(alldata);
	datamin=min(alldata);
	datamax=max(alldata);
	p=prod(datarange);
	P=(1/p);
end

total_points=nclust*npoints;
labels=zeros(total_points,1);
D=size(SPIKEDATA{1},2);
points=zeros(total_points,D);

for i=1:nclust
	startpoint=((i-1)*npoints);
	labels(startpoint+1:startpoint+npoints)=i;

	% generate using uniform density spanning entire space

	if i==nclust & garbage
		for j=1:D
			points(startpoint+1:startpoint+npoints,j)=(rand(npoints,1).*datamax(j))-datamin(j);
		end
	else
		points(startpoint+1:startpoint+npoints,:)=mvnrnd(MODEL.mu(i,:),MODEL.sigma(:,:,i),npoints);
	end
end

for i=1:nclust

	if i==nclust & garbage
		prob(:,i)=P.*ones(size(points,1),1);
	else
		prob(:,i)=mvnpdf(points,MODEL.mu(i,:),MODEL.sigma(:,:,i));
	end
end

prob=prob.*repmat(MODEL.mixing,[total_points 1]);
[~,classification]=max(prob,[],2);

for i=1:size(MODEL.mu,1)

	% how many false positives and false negatives
	%
	fp=sum((classification==i)&(labels~=i));

	% false negatives
	%
	fn=sum((classification~=i)&(labels==i));
	CONTAMINATION(i)=(fp/npoints)+(fn/npoints);
end
