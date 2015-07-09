#!/bin/bash
#
# Hunt for the .mat files in submit to the queue
#

OLD_IFS=$IFS
IFS=$'\n'

# check if the features have already been computed

COUNTER=1

# should write this out to a logfile

# pull out the list of templates

# get the template list

# construct a list of files to process for each template

while true; do

	# get the list of birds

	# re-source in case we've changed the configuration

	source ephys_pipeline_dirs.cfg

	# sort birdlist in reverse chronological order (most recently modified first)

	FILELIST=( `find $NETFILEDUMP -type f -mtime +5m` )

	for FILE in ${FILELIST[@]}; do
		echo 'Moving ' $FILE ' to ' $LOCALFILEDUMP
		mv $FILE $LOCALFILEDUMP
	done

	date 
	sleep $(($INTERVAL_FILECOPY*60))
done

IFS=$OLD_IFS


