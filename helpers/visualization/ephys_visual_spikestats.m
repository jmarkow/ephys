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

SR=25e3;
fig_num=[];
patch_color=[1 .6 0];
noise_p2p=[];
y_res=200;
spike_SR=50e3;

for i=1:2:nparams
	switch lower(varargin{i})
		case 'sr'
			SR=varargin{i+1};
		case 'spike_sr'
			spike_SR=varargin{i+1};
		case 'fig_num'
			fig_num=varargin{i+1};
		case 'patch_color'
			patch_color=varargin{i+1};
		case 'noise_p2p'
			noise_p2p=varargin{i+1};
		case 'y_res'
			y_res=varargin{i+1};
	end
end

if isempty(fig_num)
	fig_num=figure('Visible','off');
end

% patch coordinates

[samples,trials]=size(SPIKEWINDOWS);

% isi bin edges (msec)
% also plot 2D histogram

isipoints=[0:.1:10];

subplot(2,1,1);

% need the upper/lower edges for the 2D histogram

voltmin=inf;
voltmax=-inf;

timevec=([1:samples]./spike_SR)*1e3;
timevec_mat=[1:samples]';
coordmat=[];

% slows down to a CRAWL with many spikes!
% TODO re-implement with reshape and repmat

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
density=hist3(coordmat,'Edges',edges);
imagesc(timevec(1:end-1),edges{2},density(1:length(timevec)-1,:)');
colormap(hot);

%box off

if ~isempty(noise_p2p)
	mean_waveform=mean(SPIKEWINDOWS,2);
	peaktopeak=max(mean_waveform)-min(mean_waveform);
	title([' SNR:  ' num2str(peaktopeak/noise_p2p)],'FontName','Helvetica','FontSize',13)
end


xlabel('Time (ms)','FontName','Helvetica','FontSize',13);
ylabel('Voltage (in $\mu$V)','FontName','Helvetica','FontSize',13,'Interpreter','Latex');
set(gca,'YDir','Normal','tickdir','out','FontSize',10,'FontName','Helvetica');
box off
axis tight

subplot(2,1,2);
density=histc((SPIKEISI/SR)*1e3,isipoints);
h=bar(isipoints,density,'histc');
box off
set(h,'FaceColor',[.7 .7 .7],'EdgeColor','k','LineWidth',1.5);
set(gca,'TickDir','out','FontSize',10,'FontName','Helvetica');
xlabel('ISI (ms)','FontName','Helvetica','FontSize',13);
ylabel('N','FontName','Helvetica','FontSize',13);
axis tight

%linkaxes(ax,'x');
