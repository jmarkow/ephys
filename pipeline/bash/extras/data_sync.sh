#!/bin/sh

LOCKFILE=/tmp/synclock

source /usr/local/bin/ephys_pipeline_dirs.cfg

if [ -f $LOCKFILE ]; then
	exit
fi

touch $LOCKFILE

for NETWORKDIR in ${NETWORK[@]}; do
	echo "Syncing " $NETWORKDIR >> ~/.intan_sync.log>&1
	rsync --exclude "Thumbs.db" -av $NETWORKDIR $LOCAL >> ~/.intan_sync.log>&1
done

for NETWORKDIR in ${ALTNETWORK[@]}; do
	echo "Syncing " $NETWORKDIR >> ~/.intan_sync.log>&1
	rsync --exclude "Thumbs.db" -av $NETWORKDIR $ALTLOCAL >> ~/.intan_sync.log>&1
done

#rsync -av /Volumes/mic_test/ /Users/jmarkow/Desktop/test/staging/unprocessed/ >> ~/.mic_sync.log>&1

rm -f $LOCKFILE

exit

