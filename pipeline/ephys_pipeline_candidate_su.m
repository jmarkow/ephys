function SNR=ephys_candidate_su(EPHYS_DATA,CHANNELS,varargin)
%ephys_candidate_su.m approximates the SNR of channels by
%assuming a low threshold (Quiroga et al. 2004) and computing 
%the mean peak-to-peak of the spikes and comparing it to the 
%the noise peak-to-peak (Ludwig et al. 2009), SNR>1.1 is considered
%a "good" or discriminable unit.  If channels with a sufficient SNR
%are found, text files are written out with the channel number and
%the trial number with putative units
%
%	ephys_candidate_su(EPHYS_DATA,CHANNELS,varargin)
%
%	EPHYS_DATA
%	sample x trial x channel matrix (single) that contains song-aligned
%	voltage traces (stored in extracted_data.mat or aggregated_data.mat)
%
%	CHANNELS
%	vector of channel labels (stored in extracted_data.mat or aggregated_data.mat)
%
%	The following may be 
%
%
%	the following may be specified as parameter/value pairs:
%
%		fs
%		sampling rate for aligned data (25e3, default Intan)
%
%		freq_range
%		filter corner frequencies (default: 600)
%
%		snr_threshold
%		SNR threshold for "good unit" (default: 1.1)
%		
%		savedir
%		directory to write out results to
%

nparams=length(varargin);

if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs');
end

fs=25e3;
noise='none'; % common-average reference for noise removal
car_exclude=[];
filtering='y'; % if defined then filtering filter the traces
savedir=pwd;
freq_range=[800]; % frequency range for filtering
filt_type='high';
snr_threshold=8;
exclude_channels=[];
channels=CHANNELS;

for i=1:2:nparams
	switch lower(varargin{i})
		case 'fs'
			fs=varargin{i+1};
		case 'noise'
			noise=varargin{i+1};
		case 'filtering'
			filtering=varargin{i+1};
		case 'savedir'
			savedir=varargin{i+1};
		case 'snr_threshold'
			snr_threshold=varargin{i+1};
		case 'car_exclude'
			car_exclude=varargin{i+1};
	end
end

[samples,ntrials,nchannels]=size(EPHYS_DATA);

proc_data=ephys_denoise_signal(EPHYS_DATA,CHANNELS,channels,'method',noise,'car_exclude',car_exclude);
proc_data=ephys_condition_signal(proc_data,'s','freq_range',freq_range,'filt_type',filt_type);

clear EPHYS_DATA;

% collect spike windows and estimate SNR

disp('Computing SNR');

for i=1:length(channels)

	disp(['Electrode ' num2str(CHANNELS(i))])

	% collect spikes

	tmp_spikewindows={};
	tmp_spiketimes={};
	snr=[];

	for j=1:ntrials

		%disp(['Trial ' num2str(j)]);

		curr_data=[];
		spike_pp=[];
		threshold=4*median(abs(proc_data(:,j,i))/.6745);

		spikes_pp=ephys_spike_detect(proc_data(:,j,i),threshold,'fs',fs,'visualize','n','interpolate',0);
	
		% need to adjust our windows if we've interpolated

		[samples,nspikes]=size(spikes_pp.abs.windows);
	
		if nspikes<2
			snr(j)=NaN;
			continue;
		end

		% if tetrode then check each channel for significant spike

		window=ceil((samples-1)/2); % how many samples to the left and right of the spike

		% should we delete?

		adjust=0;

		% compute the peak to peak of the mean waveform
	
		%mean_waveform=mean(spikes_pp.abs.windows,2);
		%peaktopeak=max(mean_waveform)-min(mean_waveform);
		curr_data=proc_data(:,j,i);

		for k=1:nspikes

			spike_point=spikes_pp.abs.times(k)-adjust;

			if spike_point-window<1
				snr(j)=NaN;
				break;
			end

			curr_data(spike_point-window:spike_point+window)=[];
			adjust=adjust+((window*2)+1);
		end

		noise_est=std(curr_data);
		peaktopeak=max(spikes_pp.abs.windows,[],2)-min(spikes_pp.abs.windows,[],2);
		fullsnr=peaktopeak(:)./noise_est;

		snr(j)=prctile(fullsnr,90); % check the top spikes at 90th percentile to 
					    % avoid using outliers (from max)

	end

	SNR(:,i)=snr; % for all the spikes, how does the upper end of the distribution look?

	if any(snr>snr_threshold)

		disp(['Found trials with snr > ' num2str(snr_threshold)]);
		fid=fopen(fullfile(savedir,['snr_channel_' num2str(channels(i)) '.txt']),'w');
		[val,loc]=find(snr>snr_threshold);

		fprintf(fid,'Trial\tSNR(mean)\n\n');
		for j=1:length(val)
			fprintf(fid,'%g\t%0.2f\n',loc(j),snr(loc(j)));
		end

		fclose(fid);
	end
end
