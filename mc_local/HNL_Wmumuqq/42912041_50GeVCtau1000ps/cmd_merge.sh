#!/bin/bash

n1000=$1
nStart=0
nargs=$#
if (( nargs > 1 )); then
    nStart=$2
fi

for ((iJob=nStart; iJob<n1000+nStart; iJob++)); do
    idir=$(printf "%02d" $iJob)
    rm -rf merge$idir
    mkdir merge$idir
    cp template_merge/* merge$idir/
    cd merge$idir/
    sed -i "s/..\/01/..\/${idir}/g" lbexec*.yaml
    echo "nohup source cmd_all.sh > log.log 2>&1 &"
    nohup ./cmd_all.sh > log.log 2>&1 &
    sleep 1
    cd ..
done
