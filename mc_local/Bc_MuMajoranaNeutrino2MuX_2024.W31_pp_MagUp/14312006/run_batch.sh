#!/bin/bash

n1000=$1
nStart=0
nargs=$#
if (( nargs > 1 )); then
    nStart=$2
fi

for ((iJob=nStart; iJob<n1000+nStart; iJob++)); do
    idir=$(printf "%03d" $iJob)
    rm -rf $idir
    mkdir $idir
    cp template/* $idir/
    cd $idir/
    sed -i "s/10000001/1${idir}0001/g" prodConf_00012345_00006789_1.py
    echo "nohup source cmd_all.sh > log_all_timing.log 2>&1 &"
    nohup ./cmd_all.sh > log_all_timing.log 2>&1 &
    sleep 1
    cd ..
done
