#!/bin/bash

#$1 run number
#$2 eta/S tuples

RID=$(printf "%04d" $1)
echo RID got $RID
# Skip steps
lruns=1000
pruns=20000
hruns=50000

# eta value and corresponding skip steps
params=(6.28 $lruns 6.00 $lruns 5.50 $lruns 5.00 $lruns 4.50 $lruns 4.00 $lruns 3.50 $lruns 3.00 $lruns 2.50 $lruns 2.00 $lruns 1.70 $pruns 1.60 $pruns 1.50 $pruns 1.40 $pruns 1.30 $pruns 1.20 $pruns 1.10 $pruns 1.00 $pruns 0.90 $hruns 0.80 $hruns 0.70 $hruns 0.60 $hruns 0.50 $hruns 0.40 $hruns 0.30 $hruns 0.20 $hruns 0.10 $hruns)

prog="./vicsek"
outputDir="results"
seedDir="seeds_states"
mkdir -p $outputDir
mkdir -p $seedDir

#Flags - comment to disable
doMI="--MI"
doTE="--TE"
doGTE="--GTE"
#ksg_per_timestep="--calc-metrics-per-timestep"
ksg_whole_series="--calc-metrics"
#shuffle="--shuffle"

# Options
N=1000
rho=0.25
TA=5000
v=0.1
KSG_neighbours=3 # Neighbours for KSG
topo_neighbours=6 # Neighbours for topo
B=16
imethod=0
umethod=0
hist_gte_dims=1
GK_consensus_avg=0
shuffle_dim=1
subset_size=0
subset_runs=0
prevEta=${params[0]}

for ((i=0; i<${#params[@]}; i+=2)); do
	echo "${params[i]}: ${params[i+1]}"
	eta=${params[i]}
	S=${params[i+1]}

	gen="$prog generate -N $N --rho $rho -v $v --imode 0 -S $S
	--imethod $imethod --umethod $umethod --topo-neighbours ${topo_neighbours}"

	if [[ $i -eq 0 ]]; then
		$gen --eta $eta --out-file $seedDir/generate\_$eta\_$RID.bin > $seedDir/generate\_$eta\_$RID.log 2>&1
	else
		$gen --eta $eta --in-file $seedDir/generate\_$prevEta\_$RID.bin --out-file $seedDir/generate\_$eta\_$RID.bin > $seedDir/generate\_$eta\_$RID.log 2>&1
	fi

	if [[ $? -ne 0 ]]; then
		echo "Failed generating at eta=$eta"
		exit $?
	fi

	all="$prog analyse $doMI $doTE $doGTE -U $TA --KSG-neighbours ${KSG_neighbours}
	--update-initial-state $ksg_per_timestep $ksg_whole_series 
	--GTE-hist-dimensions $hist_gte_dims $shuffle --shuffle-dimension $shuffle_dim
	--GTE-consensus-average $GK_consensus_avg --subset-size $subset_size 
	--subset-runs $subset_runs -B $B"

	$all --in-file $seedDir/generate\_$eta\_$RID.bin --out-file $outputDir/\_$eta\_$RID > $outputDir/all\_$eta\_$RID.log 2>&1

	if [[ $? -ne 0 ]]; then
		echo "Failed analysis at eta=$eta"
		exit $?
	fi

	prevEta=$eta
done

