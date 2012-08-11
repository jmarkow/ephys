function ephys_multi_mua(HISTOGRAMS,MUA,varargin)
%THIS FUNCTION IS CURRENTLY DEPRECATED, USE AT YOUR OWN RISK!
%
%
%

% basic idea:
% 1) pass MUA and HISTOGRAMS as cell arrays, one cell per MUA raster and HISTOGRAM
% 2) use the histograms to find the approximate onset to align all histograms to the first
%    (so it normally makes sense to compute larger histograms to have room to maneuver)
% 3) adjust histogram imask and MUA according to the offset of the histogram onset
% 4) plot MUAs on the left, HISTOGRAMS on the right, allowing you to visualize the change
%    in the histogram and the MUA over multiple days
%

smooth_window=.0025; % smoothing window in secs
SR=40e3;
dir=pwd;
min_f=1e3;
max_f=10e3;
hist_colors=hot;
mua_colors=1-bone;
figtitle='';
onset_thresh=6e-4;

onset_f=[3e3 7e3];

nparams=length(varargin);

for i=1:2:nparams
	switch lower(varargin{i})
		case 'sr'
			SR=varargin{i+1};
		case 'dir'
			dir=varargin{i+1};
		case 'min_f'
			min_f=varargin{i+1};
		case 'max_f'
			max_f=varargin{i+1};
		case 'hist_colors'
			hist_colors=varargin{i+1};
		case 'mua_colors'
			mua_colors=varargin{i+1};
		case 'onset_f'
			onset_f=varargin{i+1};
		case 'onset_thresh'
			onset_thresh=varargin{i+1};
	end
end


if length(MUA)~=length(HISTOGRAMS)
	error('Did not pass in the same number of histograms and rasters, bailling!');
end

nplots=length(MUA)*2; %need one raster and one histogram per row

% first we need to determine the offset, align everything to the histogram with the
% earlier onset

for i=1:length(HISTOGRAMS)
	
	sdi{i}=HISTOGRAMS{i}.imask./sum(HISTOGRAMS{i}.imask(:));
	sdi_t{i}=HISTOGRAMS{i}.t;
	mua{i}=MUA{i}.image;
	mua_t{i}=MUA{i}.time;

	startidx=max([find(HISTOGRAMS{i}.f<=onset_f(1))]);
	stopidx=min([find(HISTOGRAMS{i}.f>=onset_f(2))]);

	songpower{i}=smooth(sum(sdi{i}(startidx:stopidx,:)),20);
	onsetvec=find(songpower{i}>onset_thresh);

	figure();plot(songpower{i})

	%songpower{i}=double(songpower{i}>onset_thresh);

	onsetvec(onsetvec<5)=[];

	onset(i)=min(onsetvec);
	duration(i)=length(songpower{i});

end

[min_onset ref_hist]=min(onset)
min_duration=min(duration);

% find min onset and then align everything to it

for i=setdiff(1:length(HISTOGRAMS),ref_hist)
	
	%[c,lags]=xcorr(songpower{i},songpower{ref_hist});
	%[maxc,delay]=max(c);
	%align_point=lags(delay)

	align_point=onset(i)-min_onset

	align_t=sdi_t{i}(align_point);

	% shifted version
	
	sdi_t{i}=sdi_t{i}-align_t/2;

	% now shift the MUA

	mua_t{i}=mua_t{i}-align_t/2;

end


% chop to the minimum duration, then we'll shift MUA and histograms
% histograms first, save axes for linking

multi_fig=figure('Visible','off');
counter=1;

% chop the rasters to the max or 400 trials

max_trials=-inf;

for i=1:length(MUA)

	[ntrials,samples]=size(MUA{i}.image);

	if ntrials>max_trials
		max_trials=ntrials;
	end

end

if max_trials>350
	max_trials=350;
end

for i=1:length(MUA)
	
	[ntrials,samples]=size(mua{i});

	% pad at the bottom if necessary so all images are the right aspect ratio

	if ntrials<max_trials
		mua{i}=[mua{i};zeros(max_trials-ntrials,samples)];
	end

end

for i=2:2:6
	
	ax(i)=subaxis(nplots/2,2,i,'margin',.1,'spacingvert',0.05,'spacinghoriz',.1);

	startidx=max([find(HISTOGRAMS{counter}.f<=min_f)]);
	stopidx=min([find(HISTOGRAMS{counter}.f>=max_f)]);

	imagesc(sdi_t{counter},HISTOGRAMS{counter}.f(startidx:stopidx),sdi{counter}(startidx:stopidx,:));
	set(ax(i),'ydir','normal')
	colormap(hist_colors);
	freezeColors;

	counter=counter+1;

end

counter=1;

for i=1:2:5

	% to size the MUA correctly just pad the smaller rasters with zeros near the top?
	% then eliminate the axes and done

	ax(i)=subaxis(nplots/2,2,i,'margin',.1,'spacingvert',0.05,'spacinghoriz',.1);
	imagesc(mua_t{counter},1:max_trials,mua{counter}(1:max_trials,:));
	%xlim([0 max_duration./SR]);
	colormap(mua_colors);
	freezeColors;

	counter=counter+1;

end

linkaxes(ax,'x');

set(multi_fig,'visible','on');

% use subaxis to align everything appropriately

% use 1e3-7e3 power to find onset (see Barish syllable duration)
% 
% align everything to the first crossing
% 
% chop off the boundaries if necessary

% this will also be useful to compute summary statistics...
