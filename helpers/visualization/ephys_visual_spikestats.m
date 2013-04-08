function fig_num=ephys_visual_spikestats(SPIKEWINDOWS,SPIKEISI,varargin)
%spike statistics figure, include 2D histogram and ISI distribution
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
patch_color=[1 .6 0];
noise_p2p=[];
y_res=200;
spike_fs=50e3;
note=[];
channelboundary=[];

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
		case 'note'
			note=varargin{i+1};
		case 'channelboundary'
			channelboundary=varargin{i+1};
	end
end

% TODO: spike autocorrelation and cross-correlations functions, simply xcorr of binned spikes-mean(lambda)


% patch coordinates

[samples,trials]=size(SPIKEWINDOWS);

% isi bin edges (msec)
% also plot 2D histogram

isipoints=[0:.1:100];

% need the upper/lower edges for the 2D histogram

voltmin=inf;
voltmax=-inf;

timevec=([1:samples]./spike_fs)*1e3;
timevec_mat=[1:samples]';
coordmat=[];

% reshape takes elements columnwise, so should simply have samples*trials values

spikevalues=reshape(SPIKEWINDOWS,[samples*trials 1]);

% simply repeat the time vector

samplemat=repmat(timevec_mat,[trials 1]);
coordmat=[samplemat spikevalues];

% construct 2D histogram
% for voltmin and voltmax cover 99% of voltage values

voltmin=prctile(coordmat(:,2),1)-10;
voltmax=prctile(coordmat(:,2),99)+10;

edges{1}=.5:1:samples+.5;
edges{2}=linspace(voltmin,voltmax,y_res);

if isempty(fig_num)
	fig_num=figure('Visible','on');
end

subplot(3,1,1);
plot(timevec,SPIKEWINDOWS,'m-');
ylimits=ylim();
if ~isempty(channelboundary)
	hold on;
	plot([timevec(channelboundary);timevec(channelboundary)],...
		[ylimits(1).*ones(size(channelboundary));ylimits(2).*ones(size(channelboundary))],'b--','linewidth',1.5);
end
ylabel('Voltage (µVolts)','FontName','Helvetica','FontSize',13);
set(gca,'layer','top');
box off
axis tight;

if ~isempty(noise_p2p)
	mean_waveform=mean(SPIKEWINDOWS,2);
	peaktopeak=max(mean_waveform)-min(mean_waveform);
	title([' SNR:  ' num2str(peaktopeak/noise_p2p)],'FontName','Helvetica','FontSize',13)
end

ylimits=ylim();
yticks=[ ceil(ylimits(1)/10)*10 floor(ylimits(2)/10)*10 ];

if yticks(2)<=yticks(1)
	yticks(2)=yticks(1)+1;
end

set(gca,'YTick',yticks);
prettify_axis(gca,'FontSize',12,'FontName','Helvetica');
prettify_axislabels(gca,'FontSize',15,'FontName','Helvetica');
set(gca,'xcolor',get(gcf,'color'));

yticks=[ceil(voltmin/10)*10 floor(voltmax/10)*10];

if yticks(2)<=yticks(1)
	yticks(2)=yticks(1)+1;
end

subplot(3,1,2);
density=hist3(coordmat,'Edges',edges);
imagesc(timevec(1:end-1),edges{2},density(1:length(timevec)-1,:)');
ylimits=ylim();
if ~isempty(channelboundary)
	hold on;
	plot([timevec(channelboundary);timevec(channelboundary)],...
		[ylimits(1).*ones(size(channelboundary));ylimits(2).*ones(size(channelboundary))],'y--','linewidth',1.5);
end

set(gca,'YTick',yticks);
colormap(hot);

%box off

xlabel('Time (ms)','FontName','Helvetica','FontSize',13);
prettify_axis(gca,'FontSize',12,'FontName','Helvetica');
prettify_axislabels(gca,'FontSize',15,'FontName','Helvetica');
axis xy
box off
axis tight

subplot(3,1,3);

if isempty(SPIKEISI);
	warning('ephysPipeline:visualspikestats:emptyspikeisi','ISI vector is empty, skipping ISI density plot.');
	return;
end

[density,xi]=ksdensity((SPIKEISI/fs)*1e3,isipoints,'support','positive');
%density=density./sum(density); % normalization UNNECESSARY w/ ksdensity

%density=histc((SPIKEISI/fs)*1e3,isipoints);
%h=bar(isipoints,density,'histc');

% get percentage of ISI values below 1 msec

violations=sum((SPIKEISI/fs)<.001)/length(SPIKEISI);

% round off to percent

violations=round(1000*violations)/10;

semilogx(xi,density,'r-','linewidth',3);
hline=findobj(gca,'type','line');
set(hline,'clipping','off');
box off
%set(h,'FaceColor',[.7 .7 .7],'EdgeColor','k','LineWidth',1.5);
xlabel({['ISI (ms), ' num2str(violations) '% < 1 ms '];[note]});
ylabel('P(ISI)');

ylimits(1)=min(density);
ylimits(2)=ceil(max(density)*1e2)/1e2;

if ylimits(1)<ylimits(2)
	set(gca,'YLim',ylimits,'YTick',ylimits);
end

prettify_axis(gca,'FontSize',12,'FontName','Helvetica');
prettify_axislabels(gca,'FontSize',15,'FontName','Helvetica');
set(gca,'layer','top');
xlim([xi(1) xi(end)]);
set(gca,'XTick',[ .1 1 10 100 ],'XTickLabel',[ .1 1 10 100 ],'XMinorTick','off');

%linkaxes(ax,'x');
