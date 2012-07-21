#!/bin/bash
#
# Hunt for the .mat files in submit to the queue
#

OLD_IFS=$IFS
IFS=$'\n'

# check if the features have already been computed

COUNTER=1

# should write this out to a logfile

if [ $# -lt 1 ]; then
	1=&1
fi

while true; do

	# re-source in case we've changed the configuration

	source ephys_pipeline_wrapper.cfg


	FILELIST=( `find $ROOTDIR -maxdepth 5 -name "*songdet1*_chunk_*.mat"`)

	let COUNTER=1

	for FILE in ${FILELIST[@]}; do

		BASEDIRECTORY=( `dirname $FILE `)
		FILENAME=(`basename $FILE .mat`)

		if [ ! -e "$BASEDIRECTORY/syllable_data/${FILENAME}_score.mat" ]; then

			echo $BASEDIRECTORY/syllable_data/${FILENAME}_score.mat >> $1
			NEWFILELIST[$COUNTER]=$FILE
			let COUNTER+=1

		fi	

	done

	let COUNTER-=1
	
	#echo $COUNTER ' files to compute features for' >> $1

	i="1"

	while [ $i -le ${#NEWFILELIST[@]} ]; do
		
		for ((n=1;n<=$SUBSHELLS_SMSCORE;n+=1)); do

			# insert MATLAB command here

			echo 'Spawning thread ' $n >> $1
			echo ${NEWFILELIST[$i]} >> $1
			
			$EXEC_SMSCORE "${NEWFILELIST[$i]}" $CONFIG >> $1 &

			let i+=1

			if [ ! $i -lt ${#NEWFILELIST[@]} ]; then
				break
			fi

		done

		echo 'Waiting for subshells to complete...' >> $1
		
		wait

	done

	unset NEWFILELIST

	date >> $1
	sleep $INTERVAL_SMSCORE
done

	
IFS=$OLD_IFS
