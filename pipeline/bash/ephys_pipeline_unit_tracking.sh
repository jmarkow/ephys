#!/bin/bash
#
# Hunt for the .mat files in submit to the queue
#

OLD_IFS=$IFS
IFS=$'\n'

if [ $# -lt 1 ]; then
	OUTPUT=/dev/null
else
	OUTPUT=$1
fi

while true; do

	# get the list of birds

	# re-source in case we've changed the configuration

	source ephys_pipeline_wrapper.cfg

	# need boundary points defined and spike dir

	BIRDLIST=( `ls -1 $ROOTDIR` )

	for BIRD in ${BIRDLIST[@]}; do 

		# only work with files at least 24 hours old (to prevent excessive duplication of effort)

		FILELIST=( `find $ROOTDIR/$BIRD -type f -name "*candidate_unit_ch*" -cmin +600 ` )

		for FILE in ${FILELIST[@]}; do

			FILENAME=( `basename $FILE` )
			DIRNAME=( `dirname $FILE` )
			CLUSTER=( `echo $FILENAME | sed -n 's/.*cl\([0-9][0-9]\{0,\}\).*/\1/p' ` )
			CHANNEL=( `echo $FILENAME | sed -n 's/.*ch\([0-9][0-9]\{0,\}\).*/\1/p' ` )

			CANDIDATEFILE='sua_channels '$CHANNEL'.mat'

			DONEFILE=$DIRNAME/.$CANDIDATEFILE
			CANDIDATEFILE=$DIRNAME/$CANDIDATEFILE

			# this has to be run one shell at a time

			if [ ! -e $DONEFILE ] && [ -e $CANDIDATEFILE ]; then


				# SED out cluster and channel pass to ephys_pipeline_sua_track and
				# you are done

				echo $CANDIDATEFILE >> $OUTPUT

				$EXEC_TRACK "$CANDIDATEFILE" "$CLUSTER" "$CONFIG" >> $OUTPUT &

				wait
			fi
		done
	done

	wait
	date >> $OUTPUT
	sleep $(($INTERVAL_TRACK*60))
done

IFS=$OLD_IFS

