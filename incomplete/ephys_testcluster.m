function LABELS=ephys_testcluster(SPIKEDATA,varargin)
%
%
%
%
%

merge=.22;

[coeff score]=princomp(SPIKEDATA);
SPIKEDATA=score(:,1:4);
LABELS=kmeans(SPIKEDATA,30);

% try merging now

idx=LABELS;

clusters=1:length(unique(LABELS));

if merge>0
	
	% use exp(-d), where d is Battacharyya distance to find clusters for potential merging

	count=1;

	while count>0 && length(clusters)>1

		clustpairs=nchoosek(1:length(clusters),2);
		clustdist=zeros(size(clustpairs,1),1);
		count=0;

		for i=1:size(clustpairs,1)

			c1=clustpairs(i,1);
			c2=clustpairs(i,2);

			m1=mean(SPIKEDATA(LABELS==c1,:))';
			m2=mean(SPIKEDATA(LABELS==c2,:))';
			cov1=cov(SPIKEDATA(LABELS==c1,:));
			cov2=cov(SPIKEDATA(LABELS==c2,:));

			cov1=cov1+ones(size(cov1)).*eps;
			cov2=cov2+ones(size(cov2)).*eps;

			mcov=(cov1+cov2)/2;

			dist=(1/8)*(m1-m2)'*inv(mcov)*(m1-m2)+...
				(1/2)*log((det(mcov))/(sqrt(det(cov1)*det(cov2))));	

			clustdist(i)=exp(-dist);

			fprintf('Exp(-d) between %g and %g:\t%.3f\n',c1,c2,clustdist(i));

		end

		count=sum(clustdist>merge);

		if count>0
	
			[val loc]=max(clustdist);

			merge1=clustpairs(loc(1),1);
			merge2=clustpairs(loc(1),2);

			fprintf('Merging clusters %g and %g\n',merge1,merge2);
			
			if merge1<merge2
				LABELS(LABELS==merge2)=merge1;
			else
				LABELS(LABELS==merge1)=merge2;
			end

			% ensure the labeling is contiguous again

			idx=LABELS;
			clusters=unique(idx);
			LABELS=zeros(size(LABELS));

			for i=1:length(clusters)
				LABELS(idx==clusters(i))=i;
			end

			clusters=unique(LABELS);

		end

	end

end


