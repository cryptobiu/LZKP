#!/usr/bin/bash

PORT=8000

NUM_TRIALS=10

git log -1 | head -1 > measurements.txt
date +"%Y/%m/%d %H:%M:%S.%N" >> measurements.txt


# Parameters are going in parallel, first entry of each parameter are corresponding
q=()	# e.g q=("31")
n=()	# e.g n=("256")
m=()	# e.g m=("1024")

N=()	# e.g. N=("2")
M=()	# Split with : for multiple values e.g. M=("40:80")
tau=()
echo "protocol 1" >> measurements.txt
for((p=0;p<${#q[@]};p++)); do
    MM=(${M[$i]//:/ })
    TT=(${tau[$i]//:/ })
    for((j=0;j<${#MM[@]};j++)); do

    	R=""
    	for((nt=0;nt<${NUM_TRIALS};nt++)); do
    		echo "../LZKP -p0 -i -q${q[$p]} -M${MM[$j]} -t${TT[$j]} -N${N[$p]} -n${n[$p]} -m${m[$p]} -a --port ${PORT}"
    		S=$(echo "../LZKP -p0 -i -q${q[$p]} -M${MM[$j]} -t${TT[$j]} -N${N[$p]} -n${n[$p]} -m${m[$p]} -a --port ${PORT}") # Single threaded
    		R="${R}${S} "
    		wait
    	done
    	echo ${R} >> measurements.txt

    	R=""
    	for((nt=0;nt<${NUM_TRIALS};nt++)); do
    		echo "../LZKP -p0 -i -q${q[$p]} -M${MM[$j]} -t${TT[$j]} -N${N[$p]} -n${n[$p]} -m${m[$p]} -a --port ${PORT} -x8"
    		S=$(echo "../LZKP -p0 -i -q${q[$p]} -M${MM[$j]} -t${TT[$j]} -N${N[$p]} -n${n[$p]} -m${m[$p]} -a --port ${PORT} -x8") # Multi threaded
    		R="${R}${S} "
    		wait
    	done
    	echo ${R} >> measurements.txt
    done
done


q=()
n=()
m=()

N=()
M=()
echo "protocol 2" >> measurements.txt
for((p=0;p<${#q[@]};p++)); do
    MM=(${M[$i]//:/ })
    for((j=0;j<${#MM[@]};j++)); do

    	R=""
    	for((nt=0;nt<${NUM_TRIALS};nt++)); do
    		echo "../LZKP -p1 -i -q${q[$p]} -M${MM[$j]} -t1 -N${N[$p]} -n${n[$p]} -m${m[$p]} -a --port ${PORT}"
    		S=$(echo "../LZKP -p1 -i -q${q[$p]} -M${MM[$j]} -t1 -N${N[$p]} -n${n[$p]} -m${m[$p]} -a --port ${PORT}") # Single threaded
    		R="${R}${S} "
    		wait
    	done
    	echo ${R} >> measurements.txt

    	R=""
    	for((nt=0;nt<${NUM_TRIALS};nt++)); do
    		echo "../LZKP -p1 -i -q${q[$p]} -M${MM[$j]} -t1 -N${N[$p]} -n${n[$p]} -m${m[$p]} -a --port ${PORT} -x8"
    		S=$(echo "../LZKP -p1 -i -q${q[$p]} -M${MM[$j]} -t1 -N${N[$p]} -n${n[$p]} -m${m[$p]} -a --port ${PORT} -x8") # Multi threaded
    		R="${R}${S} "
    		wait
    	done
    	echo ${R} >> measurements.txt
    done
done

date +"%Y/%m/%d %H:%M:%S.%N" >> measurements.txt