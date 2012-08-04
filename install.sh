#!/bin/bash

# REQUIREMENTS:  Mac OS X 10.6+, Matlab 2010a or later
# MATLAB Toolboxes:  Bioinformatics (SVM), Statistics (GMM), Wavelet (Spike sorting)

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
COMPILED=$BASE/pipeline/compiled/

if [ -f $DEST/ephys_pipeline_dirs.cfg ];then
	
	source $DEST/ephys_pipeline_dirs.cfg
	
	echo "Local directory:  " $LOCAL
	
	for NETDIR in ${NETWORK[@]}; do
		echo "Network directory:  " $NETDIR
	done

	echo "Compiled directory:  " $COMPILED

	echo -n "Overwrite old directories [y/n] (compiled will be overwritten for y or n)?  "
	read response

	COMPILED=$BASE/pipeline/compiled/

	case "$response" in
	"y" | "Y" ) 
		echo "Replacing directories..."

		echo -n "Enter network directory with Intan files and press [ENTER] (leave blank to set later or if you want to sync multiple directories):  "
		read networkdir

		# read out user directories

		echo -n "Enter local directory where you would like to store the Intan files and press [ENTER] (leave blank to set later):  "
		read localdir

		echo "LOCAL=$localdir" > $DEST/ephys_pipeline_dirs.cfg
		echo "COMPILED=$COMPILED" >> $DEST/ephys_pipeline_dirs.cfg
		echo "# set NETWORK to an array to sync multiple directories" >> $DEST/ephys_pipeline_dirs.cfg
		echo "# e.g. NETWORK[1]=/path/to/dir" >> $DEST/ephys_pipeline_dirs.cfg
		echo "# NETWORK[2]=/path/to/other/dir" >> $DEST/ephys_pipeline_dirs.cfg
		echo "NETWORK=$networkdir" >> $DEST/ephys_pipeline_dirs.cfg
		;;
	* ) 
		echo "Not replacing directories..."
		
		# always update the compiled directory

		sed -i '' "s|\(COMPILED\=\).*|\1$COMPILED|" $DEST/ephys_pipeline_dirs.cfg
		COUNTER=1
		;;
	esac
else

	echo "Replacing directories..."

	echo -n "Enter network directory with Intan files and press [ENTER] (leave blank to set later or if you want to sync multiple directories):  "
	read networkdir

	# read out user directories

	echo -n "Enter local directory where you would like to store the Intan files and press [ENTER] (leave blank to set later):  "
	read localdir

	echo "LOCAL=$localdir" > $DEST/ephys_pipeline_dirs.cfg
	echo "COMPILED=$COMPILED" >> $DEST/ephys_pipeline_dirs.cfg
	echo "# set NETWORK to an array to sync multiple directories" >> $DEST/ephys_pipeline_dirs.cfg
	echo "# e.g. NETWORK[1]=/path/to/dir" >> $DEST/ephys_pipeline_dirs.cfg
	echo "# NETWORK[2]=/path/to/other/dir" >> $DEST/ephys_pipeline_dirs.cfg
	echo "NETWORK=$networkdir" >> $DEST/ephys_pipeline_dirs.cfg
fi

SOURCE=$BASE/pipeline/bash

cd -- "$SOURCE"

# changed to interactive to not overwrite important files

find . -type f -maxdepth 1 -exec ln -sf -- "$SOURCE"/{} "$DEST"/{} \;

# source 2

SOURCE=$BASE/pipeline/bash/binscripts

cd -- "$SOURCE"

find . -type f -maxdepth 1 -exec ln -sf -- "$SOURCE"/{} "$DEST"/{} \;

# query the user for updating the settings files, also for updating the extras

echo -n "Do you want to replace your settings [y/n] (this will overwrite your old settings, enter y if installing for the first time)?  "
read response

# don't use symlinks for settings so that they're preserved through updates

case "$response" in
	"y" | "Y" ) 
		echo "Replacing settings..."

		rm /usr/local/bin/ephys_pipeline.cfg
		rm /usr/local/bin/ephys_pipeline_wrapper.cfg

		cp $BASE/pipeline/bash/settings/*.cfg /usr/local/bin

		;;
	* ) 
		echo "Not replacing settings..."
		;;
esac

echo -n "Do you want to overwrite the extras [y/n] (this will overwrite data_sync.sh and screenrc, enter y if installing for the first time)? "
read response

case "$response" in
	"y" | "Y" ) 
		echo "Replacing extras..."
		SOURCE=$BASE/pipeline/bash/extras
		cd -- "$SOURCE"
		find . -type f -maxdepth 1 -exec ln -sf -- "$SOURCE"/{} "$DEST"/{} \;
		;;
	* ) 
		echo "Not replacing extras..."
		;;
esac

# inform user to add /usr/local/bin/ to PATH

echo 'To finish the installation:  '
echo '1) Add /usr/local/bin/ to your PATH in ~/.profile'
echo '2) Change any settings in ephys_pipeline_wrapper.cfg or ephys_pipeline.cfg if they are not to your liking (in ' $DEST ')'
echo '3) Install the Matlab 2011b MCR and add the environment variables to ~/.profile'
echo '4) Run ephys_pipeline_wrapper.sh from the command line to start'
echo '5) Add ephys and all subdirectories to your MATLAB path'
