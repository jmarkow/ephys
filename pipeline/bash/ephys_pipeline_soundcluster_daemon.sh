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

if [ $# -lt 1 ]; then
	1=&1
fi

while true; do

	# get the list of birds

	# re-source in case we've changed the configuration

	source ephys_pipeline_wrapper.cfg

	BIRDLIST=( `ls -1 $ROOTDIR` )

	for BIRD in ${BIRDLIST[@]}; do 

		if [ ! -d $ROOTDIR/$BIRD/templates ]; then
			continue;
		fi

		FILELIST=( `find $ROOTDIR/$BIRD -name "*songdet*_chunk_*score.mat" | grep -v '/\.' | egrep -v 'ttl_score.mat'` )

		# check for files that have not been processed for this template

		# the files processed should have a dotfile in a subdir that matches
		# the template directory

		# all the files we have spectral features for

		TEMPLATELIST=( `find $ROOTDIR/$BIRD/templates -type f -name "template_data.mat" -maxdepth 2` )

		for TEMPLATE in ${TEMPLATELIST[@]}; do

			#source ephys_pipeline_wrapper.cfg

			# extract the toplevel directory from the template listing

			TEMPLATEDIRNAME=( `dirname $TEMPLATE` )
			TEMPLATETOPDIR=( `basename $TEMPLATEDIRNAME` )

			if [ ! -e "$TEMPLATEDIRNAME/classify_data.mat" ]; then
				echo 'No classify_data.mat in ' $TEMPLATEDIRNAME >> $OUTPUT
				continue;
			fi

			# exclude dot files what we want to check

			let COUNTER=1

			for FILE in ${FILELIST[@]}; do

				FILENAME=( `basename $FILE` )
				DIRNAME=( `dirname $FILE` )

				RESULTSDIR=$DIRNAME/../$TEMPLATETOPDIR

				#echo $RESULTSDIR/.${FILENAME}

				if [ ! -e $RESULTSDIR/.$FILENAME ]; then

					echo "Spawning subshell " $COUNTER >> $1

					# specify template_data file, spectral feature data file then classification object file

					nice -n $NICELVL $EXEC_CLUSTER "$TEMPLATE" "$FILE" "$TEMPLATEDIRNAME/classify_data.mat" $CONFIG >> $1 &

					let COUNTER+=1

					if [ $COUNTER -gt $SUBSHELLS_CLUSTER ]; then
						echo 'Waiting for subshells to finish' >> $1
						wait
						let COUNTER=1
					fi

				fi

			done
		done

		echo 'Checking TTL templates' >> $1

		if [ ! -d $ROOTDIR/$BIRD/templates_ttl ]; then
			continue;
		fi

		FILELIST=( `find $ROOTDIR/$BIRD -name "*songdet*_chunk_*ttl_score.mat" | grep -v '/\.'` )

		# check for files that have not been processed for this template

		# the files processed should have a dotfile in a subdir that matches
		# the template directory

		# all the files we have spectral features for

		TEMPLATELIST=( `find $ROOTDIR/$BIRD/templates_ttl -type f -name "template_data.mat" -maxdepth 2` )

		for TEMPLATE in ${TEMPLATELIST[@]}; do

			#source ephys_pipeline_wrapper.cfg

			# extract the toplevel directory from the template listing

			TEMPLATEDIRNAME=( `dirname $TEMPLATE` )
			TEMPLATETOPDIR=( `basename $TEMPLATEDIRNAME` )

			if [ ! -e "$TEMPLATEDIRNAME/classify_data.mat" ]; then
				echo 'No classify_data.mat in ' $TEMPLATEDIRNAME >> $OUTPUT
				continue;
			fi

			# exclude dot files what we want to check

			let COUNTER=1

			for FILE in ${FILELIST[@]}; do

				FILENAME=( `basename $FILE` )
				DIRNAME=( `dirname $FILE` )

				RESULTSDIR=$DIRNAME/../$TEMPLATETOPDIR

				#echo $RESULTSDIR/.${FILENAME}

				if [ ! -e $RESULTSDIR/.$FILENAME ]; then

					echo "Spawning subshell " $COUNTER >> $1

					# specify template_data file, spectral feature data file then classification object file

					nice -n $NICELVL $EXEC_CLUSTER "$TEMPLATE" "$FILE" "$TEMPLATEDIRNAME/classify_data.mat" $CONFIG >> $1 &

					let COUNTER+=1

					if [ $COUNTER -gt $SUBSHELLS_CLUSTER ]; then
						echo 'Waiting for subshells to finish' >> $1
						wait
						let COUNTER=1
					fi

				fi

			done
		done
	done


	date >> $1
	sleep $INTERVAL_CLUSTER

done

IFS=$OLD_IFS

