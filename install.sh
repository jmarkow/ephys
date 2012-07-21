#!/bin/bash

# RUN WITH SUDO!

# this install script will configure everything up to the point of installing
# the MCR libraries, which you will need to install for MATLAB 2011a 

# two directories we need to symlink are pipeline/bash and pipeline/bash/binscripts

# destination for the symlink


DEST=/usr/local/bin

if [ ! -d $DEST ]; then
	mkdir $DEST
fi

echo 'Creating symlinks in ' $DEST

# source 1

BASE=$PWD

SOURCE=$BASE/pipeline/bash

cd -- "$SOURCE"

find . -type f -maxdepth 1 -exec ln -sf -- "$SOURCE"/{} "$DEST"/{} \;

# source 2

SOURCE=$BASE/pipeline/bash/binscripts

cd -- "$SOURCE"

find . -type f -maxdepth 1 -exec ln -sf -- "$SOURCE"/{} "$DEST"/{} \;

# inform user to add /usr/local/bin/ to PATH

echo 'To finish the installation:  '
echo '1) Add /usr/local/bin/ to your PATH in ~/.profile'
echo '2) Change /usr/local/bin/ephys_pipeline_wrapper.cfg accordingly'
echo '3) Install the Matlab 2011b MCR and add the environment variables to ~/.profile'
echo '4) Run ephys_pipeline_wrapper.sh from the command line to start'
echo '5) Add ephys and all subdirectories to your MATLAB path'
