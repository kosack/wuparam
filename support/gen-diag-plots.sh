#!/bin/sh
# usage: gen-cut-plots.sh WUPARAMDIR RUNID N2ID CACHEDIR

wupdir=$1
runid=$2
n2id=$3
cachedir=$4
telid=$5

cat $wupdir/wuparam-plots.gpl \
    | sed "s:RUNID:$runid:g;s:CACHEDIR:$cachedir:g;s:N2ID:$n2id:g;s:TELID:$telid:g" \
    | gnuplot
