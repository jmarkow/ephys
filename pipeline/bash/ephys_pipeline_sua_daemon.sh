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
	
	DIRLIST=( `find $ROOTDIR -type f -name ".sua_signal" -ctime -${ELAPSED}s ` )

	# wait to see what's changed since last time

	# now in each directory compute multi-unit rasters with standalone script
	
	let COUNTER=1

	for DIR in ${DIRLIST[@]}; do

		echo 'Spawning subshell ' $COUNTER  >> $1
		echo $DIR >> $1

		PROCDIR=( `dirname $DIR`)

		$EXEC_SUA "$PROCDIR" "$CONFIG" >> $1 &	
		
		let ARRAYCOUNT=$COUNTER-1

		PIDS[$ARRAYCOUNT]=$!
		let COUNTER+=1
	
		if [ $COUNTER -gt $SUBSHELLS_SUA ]; then

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

				if [ $COUNTER -le $SUBSHELLS_SUA ]; then
					break;
				fi

				sleep .5

			done
		fi

	done

	wait
	date >> $1
	sleep $(($INTERVAL_SUA*60))

done




