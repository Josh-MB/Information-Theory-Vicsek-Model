#!/bin/bash

numCpus=10
beginRun=1
endRun=10


for r in $(seq $beginRun 1 $endRun)
do
	echo $r
	while true; do
		if [ `ps aux | grep [vi]csek | wc -l` -lt $numCpus ]; then
			break
		fi
		sleep 5
	done
	nice -20 ./singleRun.sh $r &
done
