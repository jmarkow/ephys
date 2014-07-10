function IS_LEGACY=check_legacy(FILENAME)
%checks to see if base data file uses legacy pipeline format,
%ATM simply checks to see if the variable mic_data is present
%
%

tmp=whos('-file',FILENAME);
IS_LEGACY=ismember('mic_data',{tmp(:).name});

