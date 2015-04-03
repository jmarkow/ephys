function [t,f,forcingfunction,all_sum_cross,signal_sono_mean,consensus_sono_mean,consensus_noweight_mean]=...
		lfp_chirp3(LFPDATA,SPIKETIMES,varargin)
%
%
%
%

CHIRPLET=0;

nparams=length(varargin);

if mod(nparams,2)>0
	error('ephysPipeline:argChk','Parameters must be specified as parameter/value pairs!');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETER COLLECTION %%%%%%%%%%%%%%%%%

fs=25e3;
freq_range=[100];
filt_order=5;
filt_name='butter';
proc_fs=500;
medfilt_scale=1.5;
debug=0;
ts=50:20:150; %timescales in milliseconds;
angles=-pi/4:pi/16:pi/4; 
nfft=512;
overlap=511;
consensus=0;
low_cutoff=5;
rand_lfpphase=0;
rand_spiketime=0;
contour_thresh=98; % prctile in contour length

%%%

for i=1:2:nparams
	switch lower(varargin{i})
		case 'fs'
			fs=varargin{i+1};
		case 'freq_range'
			freq_range=varargin{i+1};
		case 'filt_order'
			filt_order=varargin{i+1};
		case 'angles'
			angles=varargin{i+1};
		case 'nfft'
			nfft=varargin{i+1};
		case 'overlap'
			overlap=varargin{i+1};
		case 'ts'
			ts=varargin{i+1};
		case 'proc_fs'
			proc_fs=varargin{i+1};
		case 'debug'
			debug=varargin{i+1};
		case 'angles'
			angles=varargin{i+1};
		case 'consensus'
			consensus=varargin{i+1};
		case 'low_cutoff'
			low_cutoff=varargin{i+1};
		case 'rand_lfpphase'
			rand_lfpphase=varargin{i+1};
		case 'rand_spiketime'
			rand_spiketime=varargin{i+1};
		case 'contour_thresh'
			contour_thresh=varargin{i+1};
	end
end

% process the LFP and downsample
nfft
overlap
contour_thresh

[b,a]=butter(3,[100]/(fs/2),'low');
proc_data=filtfilt(b,a,double(LFPDATA));

downfact=fs/proc_fs;

if mod(downfact,1)>0
	error('ephysPipeline:spectcoherence:downsamplenotinteger','Need to downsample by integer');
end

disp(['Downsampling to ' num2str(proc_fs) ]);
proc_data=downsample(proc_data,downfact);

proc_data=ephys_condition_signal(proc_data,'l','freq_range',freq_range,'medfilt_scale',medfilt_scale,...
	'filt_name',filt_name,'filt_order',filt_order,'fs',proc_fs,'notch',0,'demean',1);

[nsamples,ntrials]=size(proc_data);

% bin spikes at the same frequency

binedges=[1:nsamples]./proc_fs;
binned_spikes=zeros(nsamples,ntrials,'double');

for i=1:size(proc_data,2)
	binned_spikes(:,i)=histc(SPIKETIMES{i},binedges);

	if debug
		figure(1);
		stairs(binned_spikes(:,i))
		drawnow;
		pause(.01);
	end
end

[t,f]=getspecgram_dim(nsamples,nfft,overlap,nfft,proc_fs);

%this next line is specific to this signal size

all_sum_cross=zeros(length(f),length(t));
spike_sono_mean=zeros(size(all_sum_cross));
signal_sono_mean=zeros(size(all_sum_cross));
consensus_sono_mean=zeros(size(all_sum_cross));
cross_sono_mean=zeros(size(all_sum_cross));
consensus_noweight_mean=zeros(size(all_sum_cross));

% freeze a random phase vector for all frequency bins over trials
% rand number length(f)

if rand_lfpphase

	% rand shift for all bins, 0-2pi

	phaseshift=rand(nsamples,1).*2*pi

	%  
end

for loopi=1:ntrials %loop over trials

	loopi

	signal=proc_data(:,loopi);
	spikes=double(binned_spikes(:,loopi));

	%for controls

	if rand_spiketime
		spikes=spikes(randperm(length(spikes))); % this is used to randomize the
	end
	
	%time of spikes

	
	if rand_lfpphase

		sigfft=fft(signal);
		amp=abs(sigfft);
		theta=angle(sigfft);
		z=theta+phaseshift; % include the same phaseshift for each frequency bin
		ar=amp.*cos(z)+1j*amp.*sin(z);
		signal=real(ifft(ar)); 

	end

	signalv{loopi}=signal; %record of LFP in different trials
	spikesv{loopi}=spikes;  %record of spikes in different trials

	spike_sono_tmp=zeros(size(all_sum_cross));
	signal_sono_tmp=zeros(size(all_sum_cross));
	cross_sono_tmp=zeros(size(all_sum_cross));
	consensus_sono_tmp=zeros(size(all_sum_cross));

	for tsv=1:length(ts)

		% setting first parameter>0 uses chirplet, 0 uses Gabor 

		[dxsub sonosub]=chirp2(0,ts(tsv),signal,proc_fs,overlap,nfft);

		% scramble angle here

		[dxsub_s sonosub_s]=chirp2(0,ts(tsv),spikes,proc_fs,overlap,nfft);
		cross_sono=sonosub.*conj(sonosub_s);
		%cross_sono_norm=sqrt(abs(sonosub).*abs(sonosub_s));

		signal_power=sonosub.*conj(sonosub);
		spike_power=sonosub_s.*conj(sonosub_s);

		% this can be ignored (can be 1,2)

		ANGLE_CLASS=1;

		%create a contour image combining information across all angles.

		consensus_tmp=zeros(size(all_sum_cross));
		consensus_noweight_tmp=zeros(size(all_sum_cross));

		if consensus

			for i=1:length(angles)


				a_consensus=chirp_consensus(signal,proc_fs,contour_thresh,angles(i),nfft,dxsub,0,i,ANGLE_CLASS);

				% contour length-weighted contours

				%timescale_consensus{tsv}=timescale_consensus{tsv}+(a_consensus.*cross_sono)./cross_sono_norm;
				%timescale_consensus{tsv}=timescale_consensus{tsv}+a_consensus.*sonosub; %buld average contour image
				%timescale_consensus{tsv}=timescale_consensus{tsv}+(sonosub);

				consensus_tmp=consensus_tmp+(cross_sono.*a_consensus)./length(angles);
				consensus_sono_tmp=consensus_sono_tmp+(sonosub.*a_consensus)./length(angles);
				consensus_noweight_tmp=consensus_noweight_tmp+a_consensus./length(angles);


			end

			cross_sono=consensus_tmp;

		end

		spike_sono_tmp=spike_sono_tmp+spike_power./length(ts);	
		signal_sono_tmp=signal_sono_tmp+signal_power./length(ts);
		cross_sono_tmp=cross_sono_tmp+cross_sono./length(ts);

		% merge across time-scales

		% debug figure

		if debug
			figure(1);
			imagesc(flipdim(abs(cross_sono_tmp),1));
			drawnow;
		end

	end

	spike_sono_mean=spike_sono_mean+spike_sono_tmp./ntrials;
	signal_sono_mean=signal_sono_mean+signal_sono_tmp./ntrials;
	cross_sono_mean=cross_sono_mean+cross_sono_tmp./ntrials;
	consensus_sono_mean=consensus_sono_mean+consensus_sono_tmp./ntrials;
	consensus_noweight_mean=consensus_noweight_mean+consensus_noweight_tmp./ntrials;

end

% all_sum_vector keeps track of image over each randomization

all_sum_cross=cross_sono_mean./sqrt(spike_sono_mean.*signal_sono_mean);

if debug
	figure(1);
	imagesc(flipdim(abs(all_sum_cross),1));
	drawnow;
end

forcing=zeros(length(f),length(t));

for loopi=1:ntrials
	parfor tsv=1:length(ts)
	
		signal=signalv{loopi};

		[dxsub sonosub]=chirp2(0,ts(tsv),signal,proc_fs,overlap,nfft);
		forcing=forcing+(abs(all_sum_cross).*sonosub)./(ntrials*length(ts)); %cross-spectrum weighted LFPs
	
		%in the last line, can also explore using the complex value rather than
		%abs.

	end
end


absv=abs(forcing);
preang = angle(forcing);

% phase shift every other line by pi

%for jj=1:2:length(f)
%	preang(jj,:) = preang(jj,:) + pi;
%end

%preang = mod(preang,2*pi);
sono = absv.*(cos(preang)+1i*sin(preang));

%sono(1:low_cutoff,:)=0; % high pass filter 10 Hz
forcingfunction=sum(real(sono));
