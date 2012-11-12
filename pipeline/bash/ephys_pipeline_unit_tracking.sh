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

		if [ ! -d $ROOTDIR/$BIRD/spike_templates ]; then
			continue;
		fi

		TEMPLATEBASELIST=( `find $ROOTDIR/$BIRD/spike_templates -type d -maxdepth 1` )
		
		for TEMPLATEBASE in ${TEMPLATEBASELIST[@]}; do

			# get the TEMPLATES
			# check for the CANDIDATE files, parse cluster and channel
			# if CANDIDATE matches CHANNEL, then pass CHANNEL AND CLUSTER
			# create a MATCHING DOTFILE to indicate that you're done!

			TEMPLATELIST=( `find $TEMPLATEBASE -type f -name "sua_*.mat" -maxdepth 2` )

			for TEMPLATE in ${TEMPLATELIST[@]}; do

				# extract the toplevel directory from the template listing

				TEMPLATEDIRNAME=( `dirname $TEMPLATE` )
				TEMPLATETOPDIR=( `basename $TEMPLATEDIRNAME` )

				# exclude dot files what we want to check

				if [ ! -e "$TEMPLATEDIRNAME/template.cfg" ]; then
					echo 'No template.cfg in ' $TEMPLATEDIRNAME >> $OUTPUT
					continue;
				fi

				# read in channel name

				TMPCHANNEL=( `cat $TEMPLATEDIRNAME/template.cfg | \
					grep "channel" | sed -n 's/.* \([0-9][0-9]\{0,\}\).*/\1/p' `)
				TMPCLUSTER=( `cat $TEMPLATEDIRNAME/template.cfg | \
					grep "cluster" | sed -n 's/.* \([0-9][0-9]\{0,\}\).*/\1/p' `)

				# need to read in the tempalte config and check for candidate units
				# on the same channel

				FILELIST=( `find $ROOTDIR/$BIRD -name "*candidate_unit_ch${TMPCHANNEL}*" | grep -v '/\.'` )

				let COUNTER=1

				for FILE in ${FILELIST[@]}; do

					FILENAME=( `basename $FILE` )
					DIRNAME=( `dirname $FILE` )
					CLUSTER=( `echo $FILENAME | sed -n 's/.*cl\([0-9][0-9]\{0,\}\).*/\1/p' `)

					CANDIDATEFILE='sua_channels '$TMPCHANNEL'.mat'

					# check for CELL ID before filename .CELLID_filename

					DONEFILE=$DIRNAME/.${TEMPLATETOPDIR}_$CANDIDATEFILE
					CANDIDATEFILE=$DIRNAME/$CANDIDATEFILE
					
					if [ ! -e $DONEFILE ]; then

						echo "Spawning subshell " $COUNTER >> $OUTPUT

						# SED out cluster and channel pass to ephys_pipeline_sua_track and
						# you're done

						$EXEC_TRACK "$TEMPLATE" "${TEMPLATEDIRNAME}/template.cfg" \
							"$CANDIDATEFILE" "$CLUSTER" "$CONFIG" >> $OUTPUT &

						let COUNTER+=1

						if [ $COUNTER -gt $SUBSHELLS_TRACK ]; then
							echo 'Waiting for subshells to finish' >> $OUTPUT
							wait
							let COUNTER=1
						fi

					fi

				done
			done
		done

	done

	wait
	date >> $OUTPUT
	sleep $(($INTERVAL_TRACK*60))
done

IFS=$OLD_IFS

