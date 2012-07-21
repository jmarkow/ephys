#!/bin/bash

# first let's list all birds, keep track of IDs,

# then find all syllables or motifs, we want a printout of the number of 
# extractions, finally perhaps a list of putative single units 

OLD_IFS=$IFS
IFS=$'\n'

#if [ $# -lt 1 ]; then
#1=&1
#fi

source /usr/local/bin/ephys_pipeline_wrapper.cfg

YEAR=( `date "+%Y"` )
MONTH=( `date "+%m"` )
DAY=( `date "+%d"` )

BIRDLIST=( `ls -1 $ROOTDIR 2>/dev/null` )

FORMAT="%-10s %-20s %-5d"
HEADER="%-10s %-20s %s\n"

printf "$HEADER" "Bird ID" "Extraction" "N"
printf "================================================\n\n"

for BIRD in ${BIRDLIST[@]}; do

	# concatenate year month and day to get today's listing
	# take all extractions from today to start

	EXTRACTIONS=( `find $ROOTDIR/$BIRD -type f -regex ".*$YEAR\-$MONTH\-$DAY.*" -name ".songextraction" 2>/dev/null` )

	for EXTRACT in ${EXTRACTIONS[@]}; do

		# take the base directory and the number of extractions
		# is simply the number of mat files

		EXTRACTDIR=( `dirname $EXTRACT` )

		EXTRACTNAME=( `basename $EXTRACTDIR` )

		# let's write everything out to a tmpfile then 
		# redirect to stdout when finished

		#echo $EXTRACTDIR

		FILENUM=( `ls -1 ${EXTRACTDIR}/mat/*.mat 2>/dev/null | wc -l | tr -d ' '` )

		printf "$FORMAT" "$BIRD" "$EXTRACTNAME" $FILENUM

		# list any potential single units

		SUALIST=( `find $EXTRACTDIR/mat -type f -name "sua_channel*" | sed -n 's/^.*channels \([0-9][0-9]*\).*$/\1/p' 2>/dev/null ` )

		printf "%s" "Units ("
		
		for SUA in ${SUALIST[@]}; do
			
			printf "%-5d" $SUA

		done

		printf ")\n"

	done
done


printf "\n\n================================================\n\n"
