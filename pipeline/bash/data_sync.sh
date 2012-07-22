#!/bin/sh

LOCKFILE=/tmp/synclock

source ephys_pipeline_wrapper.cfg

if [ -f $LOCKFILE ]; then
	exit
fi

touch $LOCKFILE

rsync --exclude "Thumbs.db" -av /path_to_data $ROOTDIR >> ~/.intan_sync.log>&1

rm -f $LOCKFILE

exit

