#!/bin/bash

source ephys_pipeline_wrapper.cfg

# script to be evoked by CRON

if [ $# -lt 1 ]; then
	LOGFILE=~/ephys_pipeline_aggregate.log
else
	LOGFILE=$1
fi

T="$(date +%s)"

while true; do

	# re-source in case we've changed the configuration

	source ephys_pipeline_wrapper.cfg

	ELAPSED="$(($(date +%s)-T))"

	echo 'Seconds elapsed since last run: ' $ELAPSED >> $LOGFILE

	T="$(date +%s)"

	# what has changed since the last time we ran?
	
	DIRLIST=( `find $ROOTDIR -type f -name ".songextraction" -ctime -${ELAPSED}s ` )

	# wait to see what's changed since last time
	# now in each directory compute multi-unit rasters with standalone script
	
	let COUNTER=1

	for DIR in ${DIRLIST[@]}; do

		echo 'Spawning subshell ' $COUNTER  >> $LOGFILE
		echo $DIR >> $LOGFILE

		PROCDIR=( `dirname $DIR`)

		echo $CONFIG >> $LOGFILE &

		nice -n $NICELVL $EXEC_AGGREGATE "$PROCDIR/mat" "ephys" "$CONFIG" "0" >> $LOGFILE &	
		
		let ARRAYCOUNT=$COUNTER-1

		PIDS[$ARRAYCOUNT]=$!
		let COUNTER+=1
	
		if [ $COUNTER -gt $SUBSHELLS_AGGREGATE ]; then

			echo 'Waiting for subshell to finish' >> $LOGFILE

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

				if [ $COUNTER -le $SUBSHELLS_AGGREGATE ]; then
					break;
				fi

				sleep .5

			done
		fi

	done

	DIRLIST=( `find $ROOTDIR -type d -name "sleep" -ctime +4h -mindepth 3 -maxdepth 4` )

	let COUNTER=1

	for DIR in ${DIRLIST[@]}; do

		echo 'Spawning subshell ' $COUNTER  >> $LOGFILE
		echo $DIR >> $LOGFILE

		echo $CONFIG >> $LOGFILE &

		# if we've already processed, skip

		if [ -e "$DIR/.sleep_extract" ]; then
			continue
		fi

		touch $DIR/.sleep_extract
		$EXEC_AGGREGATE "$DIR" "ephys" "$CONFIG" "1" >> $LOGFILE &
		
		let ARRAYCOUNT=$COUNTER-1

		PIDS[$ARRAYCOUNT]=$!
		let COUNTER+=1
	
		if [ $COUNTER -gt $SUBSHELLS_AGGREGATE ]; then

			echo 'Waiting for subshell to finish' >> $LOGFILE

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

				if [ $COUNTER -le $SUBSHELLS_AGGREGATE ]; then
					break;
				fi

				sleep .5

			done
		fi

	done

	wait
	date >> $LOGFILE
	sleep $(($INTERVAL_AGGREGATE*60))

done


