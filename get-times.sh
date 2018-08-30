#!/usr/bin/bash

# Usage: ./get-times.sh Log.log 

REPORTNAME=tool-times

grep 'Execution time' $1 | awk '{ print $3, $10/60, "minutes" }' > $REPORTNAME

grep 'Execution time' $1 | awk '{ print $3, $10/60, "minutes" }' | awk '{ sum+=$2 } END { print "Total", sum/60, "hours" }' >> $REPORTNAME

cat $REPORTNAME
