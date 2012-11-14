#!/bin/bash

source ephys_pipeline_wrapper.cfg

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
echo 'Log:  ' $LOG_AGGREGATE
echo 'Log:  ' $LOG_SUA
echo 'Log:  ' $LOG_MUA
echo 'Log:  ' $LOG_LFP
echo 'Log:  ' $LOG_POSTPROC
echo 'Log:  ' $LOG_TRACK

$EXEC_AGGREGATE_D $LOG_AGGREGATE &
AGGREGATE_PID=$! 

$EXEC_SUA_D $LOG_SUA &
SUA_PID=$!

$EXEC_MUA_D $LOG_MUA &
MUA_PID=$!

$EXEC_LFP_D $LOG_LFP &
LFP_PID=$!

$EXEC_POSTPROC_D $LOG_POSTPROC &
POSTPROC_PID=$!

$EXEC_TRACK_D $LOG_TRACK &
TRACK_PID=$!

echo 'Starting cluster daemon:  ' $EXEC_CLUSTER_D
echo 'Log:  ' $LOG_CLUSTER

$EXEC_CLUSTER_D $LOG_CLUSTER &
CLUSTER_PID=$!

read -p "Press any key to kill Intan daemons and exit... " -n1 -s
echo ''

kill -9 $SMSCORE_PID
kill -9 $AGGREGATE_PID
kill -9 $SUA_PID
kill -9 $MUA_PID
kill -9 $LFP_PID
kill -9 $CLUSTER_PID
kill -9 $POSTPROC_PID
kill -9 $TRACK_PID
