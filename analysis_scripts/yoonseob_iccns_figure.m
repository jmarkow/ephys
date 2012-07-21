function yoonseob_iccns_figure()
%
%
%
%
%

root_dir='/Volumes/MRJBOD/lockbox/canary/mic/2012/04/05/211558/mat/rm7_sleepplayback/mat/full_stimulus';
mua_mat=fullfile(root_dir,'mua','mua.mat');
%lfp_mat=fullfile(root_dir,'lfp_amp','300','lfp_amp_freqrange_300.mat');
stimulus='~/Desktop/Dropbox/code/research/matlab/egardner/trunk/ephys/playback/rm7/stimuli/data_20120323T151029_chunk_1_ammodulatednoise (rm7, bos).wav';
stim_mat=fullfile(root_dir,'extracted_data.mat');

load(mua_mat,'MUA');
%load(lfp_mat,'LFP_RASTER');
load(stim_mat,'mic_data');
[sound_y,sound_fs]=wavread(stimulus);

mua_fs=1./diff(MUA.t);
mua_fs=mua_fs(1);

%lfp_fs=1./diff(LFP_RASTER.t);
%lfp_fs=lfp_fs(1);

mua_channel=MUA.image(:,:,5);
mua_channel(34,:)=[];
summua=sum(mua_channel);
%summua=zscore(summua);
summua=summua-mean(summua);
%summua=detrend(summua);

sound_y=sound_y-mean(sound_y);
%sound_y=zscore(sound_y);
sound_y=detrend(sound_y);

%avelfp=mean(LFP_RASTER.image(:,:,5));
%avelfp=avelfp-mean(avelfp);
%avelfp=detrend(avelfp);

% multi-taper estimate of mua signal

[muavec,muafreqs]=mt_spectrum(summua,mua_fs,'freq_range',[1 100]);
[soundvec,soundfreqs]=mt_spectrum(sound_y,sound_fs,'freq_range',[1 100]);

%[muavec,muafreqs]=ave_fft(summua(:),mua_fs,'freq_range',[5 100]);
%[lfpvec,lfpfreqs]=ave_fft(LFP_RASTER.image(:,:,5)',lfp_fs,'freq_range',[5 100]);
%[soundvec,soundfreqs]=ave_fft(sound_y(:),sound_fs,'freq_range',[5 100]);
%[micvec,micfreqs]=ave_fft(mic_data,25e3,'freq_range',[5 100]);

figure();

ax(1)=subaxis(6,1,1,1,1,2,'spacingvert',0);
plot([1:length(mic_data(:,1))]./25e3,mic_data(:,1),'k');
ylabel('Osc.','FontSize',11,'FontName','Helvetica');
axis tight;
set(gca,'YTick',[],'XTick',[]);
ax(2)=subaxis(6,1,1,3,1,2,'spacingvert',0);
plot(MUA.t,summua./std(summua),'k');
ylabel('MU activity (zscore)','FontSize',11,'FontName','Helvetica');
xlabel('Time (sec)');
linkaxes(ax,'x');
axis tight

multi_fig_save(gcf,pwd,'signal_and_mu','eps,png');

%

figure();
plot(muafreqs,muavec,'r','linewidth',1.25);box off
axis tight;
hold on
plot(soundfreqs,soundvec,'m-.','linewidth',1.25);
%ylim([-1e-6 5e-6]);
xlabel('FS (Hz)','FontSize',20,'FontName','Helvetica');
set(gca,'TickDir','out');

multi_fig_save(gcf,pwd,'mu_mt_spectrum','eps,png');

%




