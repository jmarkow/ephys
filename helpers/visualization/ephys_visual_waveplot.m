function ephys_visual_waveplot(WINDOWS,varargin)
%for cluster quality visualization, plots mean/median waveforms against variance
%
%	ephys_visual_waveplot(WINDOWS,varargin)
%
%	WINDOWS
%	samples x trials matrix of waveforms (loaded from sua_channels x.mat, stored in clusterwindows)
%
%	the following may be passed as parameter/value pairs:
%	
%		fs
%		sampling frequency of spikes (normally twice Intan sampling rate, spikes are interpolated by default)
%
%		snr
%		snr to display as a title
%
%
%
%

nparams=length(varargin);
if mod(nparams,2)>0
	error('ephysPipeline:argChk','Parameters must be specified as parameter/value pairs!');
end

%%%

fs=50e3; % default interpolate fs
snr=[];

for i=1:2:nparams
	switch lower(varargin{i})
		case 'fs'
			fs=varargin{i+1};
		case 'snr'
			snr=varargin{i+1};
	end
end

[samples,trials]=size(WINDOWS);
time=([1:samples]./fs)*1e3;

spike_mean=mean(WINDOWS,2);
spike_upper=prctile(WINDOWS,75,2);
spike_lower=prctile(WINDOWS,25,2);

spikefig=figure();

plot(time,spike_mean,'k-','linewidth',1.3);
hold on
plot(time,spike_upper,'m--');
plot(time,spike_lower,'m--');
xlabel({'Time (in ms)'},'FontSize',18,'FontName','Helvetica','interpreter','latex');
ylabel('Voltage (in $\mu$V)','interpreter','latex','FontSize',18,'FontName','Helvetica');
set(gca,'tickdir','out','FontSize',15,'FontName','Helvetica','linewidth',1.25,'TickLength',[.025 .025]);
if ~isempty(snr)
	title(['SNR  ' num2str(snr)],'FontSize',18,'FontName','Helvetica');
end
axis tight
box off


