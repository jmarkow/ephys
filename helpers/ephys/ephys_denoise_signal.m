function [DATA,CAR]=ephys_denoise_signal(EPHYS_DATA,CHIN,CHOUT,varargin)
%takes ephys data from Intan recordings and re-references
%
%
%
%

% how to denoise the raw signal, CAR or nearest neighbor

% CHOUT is the channel mapping for the output, CHIN for the input

% i.e. if EPHYS_DATA contains channels 1,2,4,5,6 then CHIN is 1:6
% and if you only want to return channels 1 and 2 set CHOUT to 1:2


if nargin<3 | isempty(CHOUT), CHOUT=CHIN; end

nparams=length(varargin);

if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs');
end

%%%

car_exclude=[];
method='none';
car_trim=40;

for i=1:2:nparams
	switch lower(varargin{i})
		case 'car_exclude'
			car_exclude=varargin{i+1};
		case 'method'
			method=varargin{i+1};
		case 'car_trim'
			car_trim=varargin{i+1};
		otherwise

	end
end

CAR=[];

% intan nearest neighbor mapping

channel_map=[4 3 2 1 16 15 14 13 5 6 7 8 9 10 11 12];
channel_ref=[3 4 3 2 1 16 15 14 6 5 6 7 8 9 10 11];

exclude_channels=[];

for i=1:length(car_exclude)
	exclude_channels(i)=find(CHIN==car_exclude(i));
end

car_electrodes=setdiff(1:length(CHIN),exclude_channels); % which electrodes are good for CAR?
ndims_ephys=ndims(EPHYS_DATA);

[samples,ntrials,nchannels]=size(EPHYS_DATA);

% map each channel appropriately

chmap=[];
for i=1:length(CHOUT)
	idx=find(CHOUT(i)==CHIN);
	if ~isempty(idx)
		chmap=[chmap idx];
	end
end

DATA=zeros(samples,ntrials,length(chmap));

switch lower(method)

	case 'car'

		% trimmed mean to avoid subtracting in artifacts and spikes

		disp(['Using electrodes ' num2str(CHIN(car_electrodes)) ' for CAR']);
		disp(['Trimmed mean prctile ' num2str(car_trim)]);

		CAR=trimmean(EPHYS_DATA(:,:,car_electrodes),car_trim,'round',3);

		for i=1:length(chmap)
			DATA(:,:,i)=EPHYS_DATA(:,:,chmap(i))-CAR;
		end

	case 'nn'

		for i=1:length(chmap)
			DATA(:,:,i)=EPHYS_DATA(:,:,chmap(i))-EPHYS_DATA(:,:,find(channel_map==chmap(i)));
		end

	otherwise

		for i=1:length(chmap)
			DATA(:,:,i)=EPHYS_DATA(:,:,chmap(i));
		end

end



