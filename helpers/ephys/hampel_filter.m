function [REJECTS]=hampel_filter(INPUTDATA,varargin)
%MAD-based filter for removing noise
%
%
%
%

% need to add input options

[samples,trials]=size(INPUTDATA);

estimate_window=30; % length of sliding window to estimate RMS

if estimate_window>trials, estimate_window=trials-1; end
% 1.4826 converges to STD if the distribution is normal

hampel_factor=4;  % decrease to make outlier identification more aggressive

nparams=length(varargin);

if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs');
end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'estimate_window'
			estimate_window=varargin{i+1};
		case 'hampel_factor'
			hampel_factor=varargin{i+1};
	end
end

disp(['Hampel factor ' num2str(hampel_factor)]);
REJECTS=[];

for i=1:trials-estimate_window

	% use a sliding hampel filter to detect outliers

	window_rms=sqrt(mean(INPUTDATA(:,i:i+estimate_window).^2)); % RMS within the window
	window_rms_abs_dev=abs(window_rms-median(window_rms));	% MAD RMS in the window
	outlier_cutoff=hampel_factor*median(window_rms_abs_dev); % cutoff for "outliers"

	%median_rms=median(window_rms)
	%std_rms=std(window_rms)
	%[noise_sample,noise_trial]=find(abs(window_rms)>median_rms+std_rms)

	[outlier_trial]=find(window_rms_abs_dev>outlier_cutoff); % if the MAD>cutoff, CUT IT OUT!
	REJECTS=[REJECTS;outlier_trial'+i-1]; % collect the rejected trials
end

REJECTS=unique(REJECTS); % don't be redundant

disp('The following trials were rejected as noise trials...');
disp([num2str(REJECTS')]);

end%
