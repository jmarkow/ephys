function [SPECT_AVE]=ephys_visual_lfp_tf(EPHYS_DATA,HISTOGRAM,CHANNELS,varargin)
%generates a song-aligned average spectrogram of the LFP
%
%	[LFP_RASTER TIME LABEL HISTOGRAM]=ephys_lfp_amp(MIC_DATA,EPHYS_DATA,CHANNELS,varargin)
%
%	EPHYS_DATA
%	aligned ephys data generated by ephys_cluster or in extracted_data/aggregated_data.mat
%	(should be the variable ephys_data), the data should be a matrix of doubles that is 
%	samples x trials x channels
%
%	HISTOGRAM
%	contour histogram returned by ephys_visual_histogram.m (or loaded from histogram.mat)
%
%	CHANNELS
%	channel labels (i.e. the channel that corresponds to a given element in the cell array
%	ephys_data) from ephys_cluster.m or extracted_data/aggregated_data.mat
%
%	the following may be specified as parameter/value pairs:
%
%		car_exclude
%		electrodes to exclude from noise estimate
%
%		SR
%		data sampling rate (default: 25e3)
%
%		noise
%		noise rejection method ('car' for common average 'nn' for nearest neighbor, or 'none',
%		default: 'none')
%
%		freq_range
%		vector with two elements to specify the frequency range (one element specifies low pass, default: 300)
%
%		savedir
%		directory to store results (default: pwd)
%
%		hist_min_f
%		lowermost frequency to display for contour histogram (default: 1e3)
%
%		hist_max_f
%		uppermost frequency to display for contour histogram (default: 10e3)
%
%		lfp_min_f
%		lowermost frequency to display for contour histogram (default: 1)
%
%		lfp_max_f
%		uppermost frequency to display for contour histogram (default: 100)
%
%		scale
%		scale for LFP spectrogram amplitude (linear or log, default: log)
%
%		lfp_colors
%		colormap (string) for lfp data (default: jet)
%
%		hist_colors
%		colormap(string) for histogram (default: jet)
%
%		medfilt_scale
%		timescale for median filter (in ms, leave blank to skip, default: 1.5)
%
%		method
%		method for computing the LFP spectrogram ('raw' for single Hanning taper or 'mt' for 
%		multi-taper using Slepian functions, default: mt)
%
%		lfp_nfft
%		nfft for lfp spectrogram (default: 10e3)
%
%		lfp_n
%		window size for lfp spectrogram (default: 6250)
%
%		lfp_overlap
%		overlap for lfp spectrogram (default: 6000)
%
%		lfp_w
%		bandwidth for multi-taper (default: 2)
%
%		lfp_ntapers
%		override default number of tapers (leave blank for default, default: 2*w-1)
%
%		singletrials
%		number of single trial spectrograms to plot
%
% see also ephys_visual_sua.m,ephys_visual_lfp_amp.m,ephys_visual_mua.m	


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETER COLLECTION %%%%%%%%%%%%%%%%%

if nargin<2 | isempty(CHANNELS), CHANNELS=1:16; end

nparams=length(varargin);

if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs');
end

SR=25e3;
noise='none'; % common-average reference for noise removal
car_exclude=[];
savedir=pwd;
hist_colors='jet';
lfp_colors='jet';
lfp_min_f=1; % bring down to 1
lfp_max_f=100;
lfp_n=4000; % defined frequency resolution
lfp_overlap=3750;
lfp_nfft=8000; % superficial, makes the spectrogram smoother
lfp_w=2; % time bandwidth product if using multi-taper
lfp_ntapers=[]; % number of tapers, leave blank to use 2*w-1

hist_min_f=1;
hist_max_f=10e3;

figtitle=[];
freq_range=[300]; % frequency range for filtering
channels=CHANNELS;
method='mt'; % raw or mt for multi-taper
scale='log';
scalelabel='dB';
singletrials=5;
medfilt_scale=1.5; % median filter scale (in ms)

for i=1:2:nparams
	switch lower(varargin{i})
		case 'sr'
			SR=varargin{i+1};
		case 'savedir'
			savedir=varargin{i+1};
		case 'hist_colors'
			hist_colors=varargin{i+1};
		case 'car_exclude'
			car_exclude=varargin{i+1};
		case 'figtitle'
			figtitle=varargin{i+1};
		case 'freq_range'
			freq_range=varargin{i+1};
		case 'channels'
			channels=varargin{i+1};
		case 'method'
			method=varargin{i+1};
		case 'noise'
			noise=varargin{i+1};
		case 'lfp_ntapers'
			lfp_ntapers=varargin{i+1};
		case 'lfp_w'
			lfp_w=varargin{i+1};
		case 'scale'
			scale=varargin{i+1};
		case 'singletrials'
			singletrials=varargin{i+1};
		case 'lfp_min_f'
			lfp_min_f=varargin{i+1};
		case 'lfp_max_f'
			lfp_max_f=varargin{i+1};
		case 'lfp_colors'
			lfp_colors=varargin{i+1};
		case 'lfn_n'
			lfp_n=varargin{i+1};
		case 'lfp_nfft'
			lfp_nfft=varargin{i+1};
		case 'lfp_overlap'
			lfp_overlap=varargin{i+1};
	end
end


[nsamples,ntrials,nchannels]=size(EPHYS_DATA);
TIME=[1:nsamples]./SR;

% if the window size is greater than the length of the signal
% make the window size roughly the length of the signal/5
% and adjust the overlap accordingly


if lfp_n>nsamples
	difference=lfp_n-lfp_overlap;
	lfp_n=round(nsamples/5);
	lfp_overlap=lfp_n-difference;

	disp('Reset window size and overlap, longer than nsamples');
	disp(['Window size:  ' num2str(lfp_n) ' samples']);
	disp(['Overlap:  ' num2str(lfp_overlap) ' samples']);

end

% if nfft is empty set to nextpow2

if isempty(lfp_nfft)
	lfp_nfft=max([n 2^nextpow2(lfp_n)]);
else
	lfp_nfft=2^nextpow2(lfp_nfft);
end

% check if we're using multi-taper

if lower(method(1))=='m'
	resolution=lfp_w*1/(lfp_n/SR);
	
	if isempty(lfp_ntapers)
		lfp_ntapers=2*(lfp_w)-1;
	end
	[tapers,lambda]=dpss(lfp_n,lfp_w,lfp_ntapers);
else
	lfp_ntapers='';
	resolution=lfp_n/SR;
	tapers=[];
end

disp(['Resolution:  ' num2str(resolution)  ' Hz']);
disp(['NFFT:  ' num2str(lfp_nfft)]);


% get rows and columns


[t,f,lfp_startidx,lfp_stopidx]=getspecgram_dim(nsamples,lfp_n,lfp_overlap,lfp_nfft,SR,lfp_min_f,lfp_max_f);

rows=length(f);
columns=length(t);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SIGNAL CONDITIONING %%%%%%%%%%%%%%%%

% preallocate matrices

SPECT_AVE.t=t;
SPECT_AVE.f=f;
SPECT_AVE.image=zeros(length(f),length(t),length(channels),'single');

% denoise and condition the signal

proc_data=ephys_denoise_signal(EPHYS_DATA,CHANNELS,channels,'method',noise,'car_exclude',car_exclude);
clear EPHYS_DATA;
proc_data=double(ephys_condition_signal(proc_data,'l','freq_range',freq_range));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTTING CODE %%%%%%%%%%%%%%%%%%%%%%

% create output directory

[path,name,ext]=fileparts(savedir);
savedir=fullfile(savedir,'lfp_tf');

if ~exist(savedir,'dir')
	mkdir(fullfile(savedir));
end

% single trials go in a subdir

savefilename=[ name '_lfptf_electrode_'];

delete(fullfile(savedir,[savefilename '*']));

if exist(fullfile(savedir,'singletrials'));
	rmdir(fullfile(savedir,'singletrials'),'s');
end

for i=1:length(channels)

	savedir_st=fullfile(savedir,'singletrials',[ 'ch' num2str(channels(i))]);

	if ~exist(savedir_st,'dir')
		mkdir(savedir_st);
	end
	
	spect_ave_tmp=zeros(rows,columns);

	parfor j=1:ntrials
	
		currdata=proc_data(:,j,i);

		if lower(method(1))=='r'

			[spect]=spectrogram(currdata,lfp_n,lfp_overlap,lfp_nfft);
			spect=abs(spect);

		else

			% average over tapers

			spect=zeros(rows,columns);

			for k=1:lfp_ntapers

				spect_tmp=spectrogram(currdata,tapers(:,k),lfp_overlap,lfp_nfft);
				spect_tmp=spect_tmp.*conj(spect_tmp);

				% reduce instead to save memory

				spect=spect+spect_tmp./lfp_ntapers;

			end

		end
	
		% reduction to save memory, storing in a large matrix probably not viable for large n...
		
		spect_ave_tmp=spect_ave_tmp+spect./ntrials;

		if j<=singletrials
			
			singletrialfig=figure('Visible','off','position',[0 0 800 600]);

			% just plot a simple spectrogram and save


			switch lower(scale)
				case 'linear'
					imagesc(t,f(lfp_startidx:lfp_stopidx),spect(lfp_startidx:lfp_stopidx,:));
				case 'log'
					imagesc(t,f(lfp_startidx:lfp_stopidx),20*log10(spect(lfp_startidx:lfp_stopidx,:)+eps));
				otherwise
					error('Did not understand scale, must be log or linear!');
			end

			set(gca,'ydir','normal');

			axis tight;
			colormap(lfp_colors);
			set(gca,'tickdir','out','ytick',[round(lfp_min_f/10):1:round(lfp_max_f/10)]*10,...
				'linewidth',1.5,'ticklength',[.025 .025],'FontSize',16,'FontName','Helvetica');
			xlabel('Time (in s)','FontSize',13,'FontName','Helvetica');
			ylabel('Hz','FontSize',13,'FontName','Helvetica');
			box off;
			set(singletrialfig,'PaperPositionMode','auto')
			multi_fig_save(singletrialfig,savedir_st,['trial' num2str(j)],'eps,png');
			close([singletrialfig]);

			parsave(fullfile(savedir_st,[ 'trial' num2str(j)]),t,f,spect);

		end
		
	end

	SPECT_AVE.image(:,:,i)=single(spect_ave_tmp); % for saving
	spect_ave_plot.image=spect_ave_tmp./max(spect_ave_tmp(:)); % normalize so max is 0 db
	spect_ave_plot.t=SPECT_AVE.t;
	spect_ave_plot.f=SPECT_AVE.f;

	spect_fig=figure('visible','off','Units','Pixels','Position',[0 0 round(800*nsamples/SR) 800]);
	
	fig_title=['CH' num2str(channels(i)) ' NTAP' num2str(lfp_ntapers) ' RES ' num2str(resolution) ' Hz' ' NTRIALS' num2str(ntrials)];

	spect_fig=time_frequency_raster(HISTOGRAM,spect_ave_plot,'fig_num',spect_fig,'fig_title',fig_title,'scale',scale,...
		'scalelabel',scalelabel,'hist_min_f',hist_min_f,'hist_max_f',hist_max_f,'tf_min_f',lfp_min_f,'tf_max_f',lfp_max_f);

	set(spect_fig,'PaperPositionMode','auto');

	multi_fig_save(spect_fig,savedir,...
		[ savefilename num2str(channels(i)) ],'eps,png');
	close([spect_fig]);

end

save(fullfile(savedir,'lfp_tf_data.mat'),'CHANNELS','channels','SPECT_AVE');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARFOR SAVING %%%%%%%%%%%%%%%%%%%%%%

function parsave(FILE,t,f,spect)
	
save(FILE,'t','f','spect');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



