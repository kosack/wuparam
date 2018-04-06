#!/bin/sh
# usage: gen-cut-plots.sh WUPARAMDIR RUNID N2ID CACHEDIR

wupdir=$1
runid=$2

cat $wupdir/gamma-rate.gpl \
    | sed "s:RUNID:$runid:g" | gnuplot
