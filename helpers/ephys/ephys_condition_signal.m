function EPHYS_DATA=ephys_condition_signal(EPHYS_DATA,DATA_TYPE,varargin)
%conditions signals for further processing
%
%
%
%

% filter the data appropriately, single unit, multi unit or LFPs

% global defaults

fs=25e3;
medfilt_scale=1.5; % median filter scale (in ms)
filt_order=2;
sigma=.0025; % sigma for Gaussian smoothing kernel (in s)
filt_name='butter'; % default to Butterworth, use a Kaiser filter if we need
		    % sharp cutoffs, per Logothetis et al. 2001
ripple=1; % Kaiser params, ripple in dB
attenuation=.01; % Kaiser params, attenuation (not dB)
decomp_level=7;
wavelet_denoise=0;

% data type specific defaults

switch lower(DATA_TYPE)
	
	case 's'

		% single unit data

		freq_range=[800]; % 500 Hz high pass 8e3 low pass
		filt_order=2;
		filt_type='high';

		% do not de-mean,trend or median filter!

		filtering=1;
		demean=0;
		detrenddata=0;
		medfilt=0;
		rectify=0;
		smoothdata=0;


	case 'm'

		% multi unit data
		
		freq_range=[500 5e3]; %5-5k bandpass
		filt_order=2;
		filt_type='bandpass';

		filtering=1;
		demean=0;
		detrenddata=0;
		medfilt=0;
		rectify=1;
		smoothdata=1;

	case 'l'

		% lfp

		freq_range=[300]; % 300 Hz low pass
		filt_order=2;
		filt_type='low';

		filtering=1;
		demean=1;
		detrenddata=1;
		medfilt=0;
		rectify=0;
		smoothdata=0;


end

nparams=length(varargin);

if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs');
end

% overrides

for i=1:2:nparams
	switch lower(varargin{i})
		case 'fs'
			fs=varargin{i+1};
		case 'freq_range'
			freq_range=varargin{i+1};
		case 'filt_type'
			filt_type=varargin{i+1};
		case 'filt_order'
			filt_order=varargin{i+1};
		case 'filt_name'
			filt_name=varargin{i+1};
		case 'demean'
			demean=varargin{i+1};
		case 'detrenddata'
			detrenddata=varargin{i+1};
		case 'rectify'
			rectify=varargin{i+1};
		case 'smoothdata'
			smoothdata=varargin{i+1};
		case 'medfilt'
			medfilt=varargin{i+1};
		case 'medfilt_scale'
			medfilt_scale=varargin{i+1};
		case 'decomp_level'
			decomp_level=varargin{i+1};
		case 'wavelet_denoise'
			wavelet_denoise=varargin{i+1};
	end
end

if length(freq_range)>1
	filt_type='bandpass';
end

% if median filtering, this comes first

% make sure EPHYS_DATA is a double, otherwise filtering will
% bail on us


[nsamples,ntrials,nchannels]=size(EPHYS_DATA);

% cast data type if necessary

if ~isa(EPHYS_DATA,'double')
	EPHYS_DATA=double(EPHYS_DATA);
end

% order of operations:
% 1) median filter to remove spikes if necessary (usually LFP only)
% 2) demean/detrend (usually LFP only)
% 3) filter to desired frequency range
% 4) rectify (usually MU only)
% 5) smooth with a Gaussian kernel (usually MU only)

% usually only median filter the signal if we're computing fields, a timescale of ~1.5 ms 
% seems to do a good job of removing spikes

if wavelet_denoise

	disp('Wavelet denoising');

	correction=.8*sqrt(2*log(nsamples))

	% minimax threshold is a bit more conservative
	% multiply by scalar <1 if noise reduction is too severe (altered spike waveforms)

	%correction=1*(.3936+.1829*log2(nsamples));
	
	% could fold the filtering into this step as well

	for i=1:nchannels
		for j=1:ntrials
			tmp=EPHYS_DATA(:,j,i);
			[c,l]=wavedec(tmp,6,'sym7');

			% estimate sigma from cd_1
			
			lidxs=[0;cumsum(l(1:end-1))];

			for k=2:length(l)-1
				sigma=median(abs((c(lidxs(k):lidxs(k+1)))))/.6745;
				idxs=find(abs(c(lidxs(k):lidxs(k+1)))<=correction*sigma);
				c(idxs+lidxs(k)-1)=0;
			end

			EPHYS_DATA(:,j,i)=waverec(c,l,'sym7');
		end

	end

end

if medfilt	

	disp(['Median filtering, timescale:  ' num2str(medfilt_scale) ' ms']);

	medfilt_order=round(medfilt_scale/1e3*fs);

	for i=1:nchannels
		parfor j=1:ntrials
			EPHYS_DATA(:,j,i)=medfilt1(EPHYS_DATA(:,j,i),medfilt_order);
		end
	end
end

% demean

if demean

	disp('Demeaning data...');

	for i=1:nchannels
		for j=1:ntrials
			EPHYS_DATA(:,j,i)=EPHYS_DATA(:,j,i)-mean(EPHYS_DATA(:,j,i));
		end
	end
end

% detrending

if detrenddata

	disp('Detrending (removing linear trend across each trial)');

	for i=1:nchannels
		for j=1:ntrials
			EPHYS_DATA(:,j,i)=detrend(EPHYS_DATA(:,j,i));
		end
	end
end

% filtering

if filtering

	switch lower(filt_name(1))

		case 'b'	

			disp('Butterworth filter');
			disp(['Filter order ' num2str(filt_order)]);
			disp(['Frequency range ' num2str(freq_range)]);
			disp(['Filter type ' filt_type]);

			[b,a]=butter(filt_order,[freq_range]/(fs/2),filt_type);

			% zero-phase filter here

			for i=1:nchannels
				EPHYS_DATA(:,:,i)=filtfilt(b,a,squeeze(EPHYS_DATA(:,:,i)));
			end
	
		case 'e'
			disp('Elliptic filter');


			[b,a]=ellip(2,.1,40,[freq_range]/(fs/2));


			for i=1:nchannels
				EPHYS_DATA(:,:,i)=filtfilt(b,a,squeeze(EPHYS_DATA(:,:,i)));
			end

		case 'k'
		

			disp('Kaiser filter');
			disp(['Frequency cutoffs ' num2str(freq_range)]);
			disp(['Ripple ' num2str(ripple) ' dB']);
			disp(['Attenuation ' num2str(20*log10(attenuation)) ' dB']);
		
			switch lower(filt_type(1))
				case 'l'
					mags=[1 0];
					dev=[ 10^(ripple/20)-1 attenuation ];
				case 'h'
					mags=[0 1];
					dev=[ attenuation 10^(ripple/20)-1 ];
				case 'b'
					mags=[0 1 0];
					dev=[ attenuation 10^(ripple/20)-1 attenuation];
				otherwise
					error('Did not understand filter type!');
			end

			[M,Wn,beta,typ]=kaiserord(freq_range,mags,dev,fs);
			b=fir1(M,Wn,typ,kaiser(M+1,beta),'noscale');

			disp(['Filter order ' num2str(length(b)) ]);

			for i=1:nchannels
				EPHYS_DATA(:,:,i)=filtfilt(b,1,squeeze(EPHYS_DATA(:,:,i)));
			end
	
		case 'w'

			disp('Wavelet filter');
			disp(['Decomposition level ' num2str(decomp_level)]);
			disp(['Frequency cutoff ' num2str(fs/(2^decomp_level))]);

			for i=1:nchannels
				
				for j=1:ntrials
					
					[c,l]=wavedec(EPHYS_DATA(:,j,i),decomp_level,'sym7');

					% set the approximation coefficients to zero

					c(1:l(1))=0;

					% reconstruct

					EPHYS_DATA(:,j,i)=waverec(c,l,'sym7');
				end
			end



		otherwise
		
			error('Did not understand filter name (must be kaiser or butter)');
		
	end

end

% rectification

if rectify

	disp('Rectifying data (squaring)');

	for i=1:nchannels
		EPHYS_DATA(:,:,i)=EPHYS_DATA(:,:,i).^2;
	end

end

% finally, smoothing

if smoothdata
	
	disp(['Smoothing data with Gaussian kernel, sigma=' num2str(sigma*1e3) ' ms']);
	
	edges=[-3*sigma:1/fs:3*sigma];
	kernel=(1/(sigma*sqrt(2*pi)))*exp((-(edges-0).^2)./(2*sigma^2));
	kernel=kernel./sum(kernel); % normalize to sum to 1

	for i=1:nchannels
		parfor j=1:ntrials
			EPHYS_DATA(:,j,i)=conv(EPHYS_DATA(:,j,i),kernel,'same');
		end
	end

end





