#!/bin/bash

OLD_IFS=$IFS
IFS=$'\n'

#source intan_songdet.cfg

if [ $# == 1 ]; then
	EXTRACTLIST=( `find $PWD -type f -name ".songextraction"` )
elif [ $# == 2 ]; then
	EXTRACTLIST=( `find $PWD -type f -name ".songextraction" | grep "$2"` )
else
	echo 'Did not understand command line options'
	echo 'Standard usage: ' $0 ' mua|lfp|sua extraction_name '
	exit 1
fi

case "$1" in
	"sua" ) REPLACEFILE="mat/ephys/.sua_signal";;
	"mua" ) REPLACEFILE="mat/ephys/.mua_signal";;
	"lfp" ) REPLACEFILE="mat/ephys/.lfp_signal";;
	"extract" ) REPLACEFILE=".songextraction";;
	* ) 
		echo 'Did not understand argument, bailing...'
		echo 'Possible options are "sua" "mua" "lfp" or "extract"'
		echo 'Pass as the first argument when calling the script'
		exit 1;;
esac

for EXTRACT in ${EXTRACTLIST[@]}; do
	EXTRACTDIR=( `dirname $EXTRACT` )
	REPLACEDIR=( `dirname $EXTRACTDIR/REPLACEFILE` )

	if [ -d "$REPLACEDIR" ]; then
		echo 'Touching ' $EXTRACTDIR/$REPLACEFILE
		touch $EXTRACTDIR/$REPLACEFILE
	fi
done

IFS=$OLD_IFS

