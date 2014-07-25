function IS_LEGACY=check_legacy(FILENAME)
%checks to see if base data file uses legacy pipeline format,
%ATM simply checks to see if the variable mic_data is present
%
%

tmp=whos('-file',FILENAME);
flag1=ismember('mic_data',{tmp(:).name});
flag2=ismember('audio_extraction',{tmp(:).name});
flag3=ismember('template_features',{tmp(:).name});

IS_LEGACY=flag1|flag2|flag3;

