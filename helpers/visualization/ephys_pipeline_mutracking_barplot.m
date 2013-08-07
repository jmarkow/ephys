function ephys_pipeline_mutracking_barplot(CONTROL,DOWN,UP,varargin)
%
%
%
%
%
%

% creates a scatter plot of multi-unit activity against time with a smooth average

% take the data, square and take the average

intan_fs=25e3;
freq_range=[500 3e3];
plotcolor=[.7 .7 .7];
windowsize=15;
mads=30;
grpid=[];
grpcolors=[0 0 0;...
	0 0 1;...
	1 0 0];
startid=1; % which datenum indicates feedback start?

nparams=length(varargin);

if mod(nparams,2)>0
	error('ephysPipeline:argChk','Parameters must be specified as parameter/value pairs!');
end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'grpid'
			grpid=varargin{i+1};
		case 'startid'
			startid=varargin{i+1};
	end
end

% bandpass filter to recover the multi-unit data

[b,a]=ellip(5,.2,40,freq_range./(intan_fs/2),'bandpass');

mucontrol=mean(CONTROL);
stdcontrol=std(CONTROL);

zcontrol=(CONTROL-mucontrol)./(stdcontrol);
zdown=(DOWN-mucontrol)./(stdcontrol);
zup=(UP-mucontrol)./(stdcontrol);

zcontrol_mu=mean(zcontrol);
zcontrol_std=std(zcontrol);

zdown_mu=mean(zdown);
zdown_std=std(zdown);

zup_mu=mean(zup);
zup_std=std(zup);

minval=min([CONTROL DOWN UP]);
maxval=max([CONTROL DOWN UP]);

bins=[minval:.2:maxval inf]'; 

x=histc(CONTROL,bins);
y=histc(UP,bins);
z=histc(DOWN,bins);

x=x./sum(x);
y=y./sum(y);
z=z./sum(z);
figure();

xvec=[];
yvec=[];
zvec=[];
binvec=[];

for i=1:length(bins)-1
	binvec=[binvec;bins(i:i+1)];
	xvec=[xvec;repmat(x(i),2,1)];
	yvec=[yvec;repmat(y(i),2,1)];
	zvec=[zvec;repmat(z(i),2,1)];
end

ax(1)=stairs(binvec,xvec,'k-','linewidth',3);hold on;
ax(2)=stairs(binvec,yvec,'r-','linewidth',3);
ax(3)=stairs(binvec,zvec,'b-','linewidth',3);
xlim([11.5 17]);

box off;
axis tight;
ylabel('P','FontSize',20,'FontName','Helvetica');
xlabel('MU activity','FontSize',20,'FontName','Helvetica');
set(gca,'layer','top','linewidth',2,'ticklength',[.025 .025],...
	'tickdir','out','fontsize',20,'fontname','helvetica');

L=legend(ax,'No FB','Escape up','Escape down');
legend boxoff;
set(L,'FontSize',20,'FontName','Helvetica','Location','NorthEast');
%%%%% bar chart

% plot means and 3*sem

barheights=[zcontrol_mu zdown_mu zup_mu];
barerr=[ zcontrol_std zdown_std zup_std ];
figure();

for i=1:3
	x=i-1+.75;
	h=bar(x,barheights(i),1);
	hold on;
	plot([x;x],[barheights(i)-barerr(i);barheights(i)+barerr(i)],'k-','linewidth',1.5)
	set(h,'FaceColor',grpcolors(i,:),'EdgeColor','none');
end
xlim([0 3.5]);


box off;
set(gca,'layer','top','linewidth',2,'ticklength',[.025 .025],...
	'tickdir','out','fontsize',20,'fontname','helvetica','xtick',[]);
ylabel('MU activity (z-score)');



