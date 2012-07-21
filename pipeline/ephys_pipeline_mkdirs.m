function ephys_pipeline_mkdirs()
%ephys_pipe_mkdirs.m creates the directory structure for the ephys pipeline:
%
%pwd/staging/processed
%pwd/staging/unprocessed
%pwd/data/intan_data
%
%the Intan files should be stored in pwd/staging/unprocessed

staging=fullfile(pwd,'staging');
data=fullfile(pwd,'data');

mkdir(fullfile(staging,'unprocessed'));
mkdir(fullfile(staging,'processed'));
mkdir(fullfile(data,'intan_data'));
