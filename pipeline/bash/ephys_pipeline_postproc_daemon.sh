#!/bin/bash

# script to be evoked by CRON

if [ $# -lt 1 ]; then
	1=&1
fi

T="$(date +%s)"

while true; do

	# re-source in case we've changed the configuration

	source ephys_pipeline_wrapper.cfg

	ELAPSED="$(($(date +%s)-T))"

	echo 'Seconds elapsed since last run: ' $ELAPSED >> $1


	T="$(date +%s)"

	# what has changed since the last time we ran?

	# set mindepth to 5 to prevent neuron lists from being sucked in...

	UNITLIST=( `find $ROOTDIR -type f -name "*.postproc" -ctime -${ELAPSED}s -mindepth 5` )

	# wait to see what's changed since last time

	# now in each directory compute multi-unit rasters with standalone script

	
	let COUNTER=1

	for UNITFILE in ${UNITLIST[@]}; do

		echo 'Spawning subshell ' $COUNTER  >> $1
		echo $DIR >> $1

		#PROCDIR=( `dirname $DIR`)

		$EXEC_POSTPROC "$UNITFILE" "$CONFIG" >> $1 &	
		
		let ARRAYCOUNT=$COUNTER-1

		PIDS[$ARRAYCOUNT]=$!
		let COUNTER+=1
	
		if [ $COUNTER -gt $SUBSHELLS_POSTPROC ]; then

			echo 'Waiting for subshell to finish' >> $1

			let OLDCOUNT=$COUNTER

			while true; do

				let COUNTER=$OLDCOUNT
				let PIDCOUNT=0

				for pid in ${PIDS[@]}; do

					kill -0 $pid > /dev/null 2>&1

					if [ "$?" -eq 1 ]; then
						let COUNTER-=1
						PIDS=(${PIDS[@]:0:$PIDCOUNT} ${PIDS[@]:$(($PIDCOUNT + 1))})
					else
						let PIDCOUNT+=1
					fi
				done

				if [ $COUNTER -le $SUBSHELLS_POSTPROC ]; then
					break;
				fi

				sleep .5

			done
		fi

	done

	wait
	date >> $1
	sleep $(($INTERVAL_POSTPROC*60))

done





