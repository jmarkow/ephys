function greg_proc_data()
%
%
%
%
%
%

root_dir='/Volumes/MRJBOD/workspace/other/greg/Examples of Acute Data/';

listing=dir(root_dir);
proc_dirs={};

for i=1:length(listing)
	if listing(i).name(1)~='.' && listing(i).isdir
		proc_dirs{end+1}=fullfile(root_dir,listing(i).name);
	end
end

% specify TDT or Intan
% specify CAR channels and spike channels

ref_channel_tdt=[];
ref_channel_intan=[];

intan_fs=25e3;
tdt_fs=24414;

% subtract CAR, filter if necessary then use sensible threshold (e.g. 7*std)
% and compute SNR 

for i=1:length(proc_dirs)

	disp([proc_dirs{i}]);

	savedir=fullfile(proc_dirs{i},'analysis');
	
	proceed=1;

	if ~exist(savedir,'dir');
		mkdir(savedir);
	else
		response=[];
		while isempty(response)
			response=input('Directory already processed, skip (y or n)?  ','s');
			switch lower(response(1))
				case 'y'
					proceed=0;
				case 'n'
					proceed=1;
				otherwise
					response=[];
			end
		end

	end

	if ~proceed
		continue;
	end


	if ~exist(fullfile(proc_dirs{i},'spikes.txt'),'file')
		disp('Did not find spikes file...');
		continue;
	end

	if ~exist(fullfile(proc_dirs{i},'CAR.txt'),'file')
		disp('Did not find CAR file...');
		continue;
	end

	if ~exist(fullfile(proc_dirs{i},'info.txt'),'file')
		disp('Did not find info file...');
		birdid='Pretty bird';
		area='';
		misc='';
	else
		fid=fopen(fullfile(proc_dirs{i},'info.txt'),'r');
		birdid=fgetl(fid);
		area=fgetl(fid);
		misc=fgetl(fid);

		if ~ischar(misc)
			misc='';
		end

		fclose(fid);

	end

	% read in the spike channels

	fid=fopen(fullfile(proc_dirs{i},'spikes.txt'),'r');
	spike_channels=fscanf(fid,'%d');
	fclose(fid);

	% read in CAR channels

	fid=fopen(fullfile(proc_dirs{i},'CAR.txt'),'r');
	car_channels=fscanf(fid,'%d');
	fclose(fid);

	% are we working with Intan or TDT

	intan_files=dir(fullfile(proc_dirs{i},'*.int'));
	tdt_files=dir(fullfile(proc_dirs{i},'*.mat'));

	if ~isempty(intan_files) && ~isempty(tdt_files)
		error('Found both Intan and TDT files!');
	end


	if ~isempty(intan_files)

		[t,amps,data]=read_intan_data_cli(fullfile(proc_dirs{i},intan_files(1).name));

		samples=length(t);

		car_mat=zeros(samples,length(car_channels));
		spike_data=zeros(samples,length(spike_channels));

		for j=1:length(car_channels)
			car_mat(:,j)=data(:,find(amps==car_channels(j)));
		end

		for j=1:length(spike_channels)
			spike_data(:,j)=data(:,find(amps==spike_channels(j)));;
		end

		fs=intan_fs;

	else

		% odd file format...
		load(fullfile(proc_dirs{i},tdt_files(1).name));

		samples=length(ChRawDat{1}.rawdat);

		car_mat=zeros(samples,length(car_channels));
		spike_data=zeros(samples,length(spike_channels));

		for j=1:length(car_channels)
			car_mat(:,j)=ChRawDat{car_channels(j)}.rawdat;
		end

		for j=1:length(spike_channels)
			spike_data(:,j)=ChRawDat{spike_channels(j)}.rawdat;
		end

		fs=tdt_fs;

		% process each spike channel now
	end





	[b,a]=butter(2,[10e3]/(fs/2),'low');

	car=mean(car_mat,2);
	clear car_mat;

	% re-filter to be safe

	for j=1:length(spike_channels)


	
		curr_data=filtfilt(b,a,spike_data(:,j)-car(:));

		if exist(fullfile(proc_dirs{i},'avoid.txt'),'file')
			fid=fopen(fullfile(proc_dirs{i},'avoid.txt'),'r');
			chop_idxs=fscanf(fid,'%e');
			curr_data(chop_idxs(1):chop_idxs(2))=[];
			fclose(fid);
		end

		% subtract CAR, choose sensible threshold, collect spike times
		% then delete the spikes and compute SNR

		threshold_file=fullfile(proc_dirs{i},['threshold_' num2str(spike_channels(j)) '.txt' ]);

		if exist(threshold_file,'file')
			fid=fopen(threshold_file,'r');
			threshold=fscanf(fid,'%e');
			fclose(fid);
		else
			data_display=figure();
			plot([1:length(curr_data)]./fs,curr_data);
			title('Raw data');
			response=[];
			
			while isempty(response)
				response=input('Please enter a threshold:  ','s');
				threshold=str2num(response);
			end
			
			fid=fopen(threshold_file,'w');
			fprintf(fid,'%f',threshold);
			fclose(fid);

			close([data_display]);

		end	


		% threshold, may want to use Quiroga's threshold, let's test 
		% with this

		noise_events=find(curr_data>1000);

		if ~isempty(noise_events)
			curr_data(noise_events(1):noise_events(end))=[];
		end

		%threshold=14*std(curr_data);

		disp('Detecting spikes..');
		disp(['Threshold ' num2str(threshold)]);

		[spikes_pp spikes_pb]=ephys_spike_detect(curr_data,threshold,...
			'sr',fs,'visualize','n','window',[.001 .002],'align','max');
		%=spikes_pp.abs.windows;

		[samples,nspikes]=size(spikes_pp.abs.windows)

		if nspikes<2
			continue;
		end

		window=floor((samples-1)/2); % how many samples to the left and right of the spike
		% should we delete?

		adjust=0;

		% compute the peak to peak of the mean waveform

		mean_waveform=mean(spikes_pp.abs.windows,2);
		peaktopeak=max(mean_waveform)-min(mean_waveform);
		raw_data=curr_data;

		for k=1:nspikes
			spike_point=spikes_pp.abs.times(k)-adjust;

			if spike_point<1
				break;
			end

			curr_data(spike_point-window:spike_point+window)=[];
			adjust=adjust+((window*2)+1);
		end

		% rms

		rms=sqrt(mean(curr_data.^2));

		% kipke critical value, assuming noise is white then 6*std
		% accounts for 99.97% of the distribution, thus is a good 
		% approximation of noise p2p

		kipke=6*std(curr_data);

		snr=mean(peaktopeak)./rms;
		snr_kipke=mean(peaktopeak)./kipke;

		

		greg_spike_plot(spikes_pp.abs.windows,raw_data,threshold,...
			fs,curr_data,birdid,area,misc,snr_kipke,spike_channels(j),savedir);

		% write out the stats to a text file

	end


end

end

function greg_spike_plot(waveforms,rawdata,threshold,fs,rawdata_nospike,birdid,area,notes,snr,channel,savedir)
%
%
%

% plot:
%
% 1, raw data with threshold overlaid
% 2, mean waveform with variance shown
% 3, raw data without the spikes

time=[1:length(rawdata)]./fs;

% we may want to cycle through two-five second windows and generate figures...

% 5 second segments

segments=1:5*fs:length(rawdata);

if length(segments)-1>15
	segments=segments(1:16);
end


rawfig=figure('Visible','off');
plot(time,rawdata,'k-');
hold on
plot(time,ones(size(time)).*threshold,'k--','color','r','linewidth',1.3);
plot(time,ones(size(time)).*-threshold,'k--','color','r','linewidth',1.3);
ylimits=ylim();
xlabel({'Time (in s)'},'FontSize',18,'FontName','Helvetica','interpreter','latex');
ylabel('Voltage (in $\mu$V)','FontSize',18,'FontName','Helvetica','interpreter','latex')
title({[ birdid ' ' area];[notes];['Channel ' num2str(channel) ' threshold ' num2str(threshold)]},'FontSize',18,'FontName','Helvetica');
box off
set(gca,'tickdir','out','FontSize',13,'FontName','Helvetica','linewidth',1.25,'TickLength',[.025 .025]);
axis tight

set(gcf,'PaperPositionMode','auto');
multi_fig_save(rawfig,savedir,['channel_' num2str(channel) '_rawdata'],'eps');
close([rawfig]);

for i=2:length(segments)

	disp(['Plotting segment ' num2str(i-1) ]);

	segmentfig=figure('Visible','off');
	plot(time(segments(i-1):segments(i)),rawdata(segments(i-1):segments(i)),'k-');

	xlabel({'Time (in s)'},'FontSize',18,'FontName','Helvetica','interpreter','latex');
	ylabel('Voltage (in $\mu$V)','FontSize',18,'FontName','Helvetica','interpreter','latex')
	title({[ birdid ' ' area];[notes];['Channel ' num2str(channel) ' threshold ' num2str(threshold)]},...
		'FontSize',18,'FontName','Helvetica');
	box off
	set(gca,'tickdir','out','FontSize',13,'FontName','Helvetica','linewidth',1.25,'TickLength',[.025 .025]);
	axis tight
	ylim(ylimits)

	set(segmentfig,'PaperPositionMode','auto');
	multi_fig_save(segmentfig,savedir,['channel_' num2str(channel) '_segment_' num2str(i-1) '_rawdata'],'eps');
	close([segmentfig]);

end

wavefig=figure('Visible','off');

[samples,nspikes]=size(waveforms);
time_waveforms=[1:samples]./fs*1e3;

% plot the median and 25th-->75th percentile

spike_median=median(waveforms,2);
spike_lower=prctile(waveforms',25);
spike_upper=prctile(waveforms',75);

patch_x=time_waveforms;
patch_x=[patch_x,fliplr(patch_x)];
patch_y=[spike_lower,fliplr(spike_upper)];
patch_z=-eps*ones(size(patch_x));
patch_color=[1 .6 0];

%patch(patch_x,patch_y,patch_z,'facecolor',patch_color,'edgecolor','none')
plot(time_waveforms,spike_median,'k-','linewidth',1.3);
hold on
plot(time_waveforms,spike_upper,'m--');
plot(time_waveforms,spike_lower,'m--');
xlabel({'Time (in ms)'},'FontSize',18,'FontName','Helvetica','interpreter','latex');
ylabel('Voltage (in $\mu$V)','interpreter','latex','FontSize',18,'FontName','Helvetica');
set(gca,'tickdir','out','FontSize',15,'FontName','Helvetica','linewidth',1.25,'TickLength',[.025 .025]);
title(['SNR  ' num2str(snr)],'FontSize',18,'FontName','Helvetica');
axis tight
box off

set(gcf,'PaperPositionMode','auto');
multi_fig_save(wavefig,savedir,['channel_' num2str(channel) '_spike_waveforms'],'eps');
close([wavefig]);


noisefig=figure('Visible','off');

time=[1:length(rawdata_nospike)]./fs;
plot(time,rawdata_nospike,'k-');

% noise estimate

peaktopeaknoise=6*std(rawdata_nospike);
hold on
plot(time,(peaktopeaknoise/2).*ones(size(time)),'k--','color','r','linewidth',1.3);
plot(time,(-peaktopeaknoise/2).*ones(size(time)),'k--','color','r','linewidth',1.3);

xlabel({'Time (in s)'},'FontSize',18,'FontName','Helvetica','interpreter','latex');
ylabel('Voltage (in $\mu$V)','FontSize',18,'FontName','Helvetica','interpreter','latex')
title({birdid;'Data with spikes removed';[ 'peak-to-peak:  ' num2str(peaktopeaknoise)]},'FontSize',18,'FontName','Helvetica');
set(gca,'tickdir','out','FontSize',13,'FontName','Helvetica','linewidth',1.25,'TickLength',[.025 .025]);
box off
xlim([time(1) time(end)]);
ylim(ylimits);

set(gcf,'PaperPositionMode','auto');
multi_fig_save(noisefig,savedir,['channel_' num2str(channel) '_noise_estimate'],'eps');
close([noisefig]);



end
