#!/bin/bash

source ephys_pipeline_wrapper.cfg

trap 'jobs -p | xargs kill'  EXIT

echo 'Starting smscore daemon: ' $EXEC_SMSCORE_D
echo 'Log:  ' $LOG_SMSCORE
 
$EXEC_SMSCORE_D $LOG_SMSCORE &
SMSCORE_PID=$!

echo 'Starting ephys daemons:  ' $EXEC_AGGREGATE_D
echo 'Starting ephys daemons:  ' $EXEC_SUA_D
echo 'Starting ephys daemons:  ' $EXEC_MUA_D
echo 'Starting ephys daemons:  ' $EXEC_LFP_D
echo 'Starting ephys daemons:  ' $EXEC_POSTPROC_D
echo 'Starting ephys daemons:  ' $EXEC_TRACK_D
echo 'Starting ephys daemons:  ' $EXEC_FILECOPY_D
echo 'Log:  ' $LOG_AGGREGATE
echo 'Log:  ' $LOG_SUA
echo 'Log:  ' $LOG_MUA
echo 'Log:  ' $LOG_LFP
echo 'Log:  ' $LOG_POSTPROC
echo 'Log:  ' $LOG_TRACK
echo 'Log:  ' $LOG_FILECOPY

$EXEC_AGGREGATE_D $LOG_AGGREGATE &
$EXEC_SUA_D $LOG_SUA &
$EXEC_MUA_D $LOG_MUA &
$EXEC_LFP_D $LOG_LFP &
$EXEC_POSTPROC_D $LOG_POSTPROC &
$EXEC_TRACK_D $LOG_TRACK &
$EXEC_FILECOPY_D $LOG_FILECOPY &

echo 'Starting cluster daemon:  ' $EXEC_CLUSTER_D
echo 'Log:  ' $LOG_CLUSTER

$EXEC_CLUSTER_D $LOG_CLUSTER &

read -p "Press any key to kill Intan daemons and exit... " -n1 -s
echo ''
