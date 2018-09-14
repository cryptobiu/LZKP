#!/usr/bin/bash

PORT=8000

NUM_TRIALS=10

git log -1 | head -1 > measurements.txt
date +"%Y/%m/%d %H:%M:%S.%N" >> measurements.txt


q=()
n=()
m=()

N=()
M=()
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
    		echo "../LZKP -p0 -i -q${q[$p]} -M${MM[$j]} -t${TT[$j]} -N${N[$p]} -n${n[$p]} -m${m[$p]} -a --port ${PORT} -x"
    		S=$(echo "../LZKP -p0 -i -q${q[$p]} -M${MM[$j]} -t${TT[$j]} -N${N[$p]} -n${n[$p]} -m${m[$p]} -a --port ${PORT} -x") # Multi threaded
    		R="${R}${S} "
    		wait
    	done
    	echo ${R} >> measurements.txt
    done
done


q=("15" "15" "31" "61")
n=("256" "256" "512" "1024")
m=("1024" "4096" "2048" "4096")

N=("2" "4" "8" "16")
M=("40:80" "20:40" "14:27" "10:20")
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
    		echo "../LZKP -p1 -i -q${q[$p]} -M${MM[$j]} -t1 -N${N[$p]} -n${n[$p]} -m${m[$p]} -a --port ${PORT} -x"
    		S=$(echo "../LZKP -p1 -i -q${q[$p]} -M${MM[$j]} -t1 -N${N[$p]} -n${n[$p]} -m${m[$p]} -a --port ${PORT} -x") # Multi threaded
    		R="${R}${S} "
    		wait
    	done
    	echo ${R} >> measurements.txt
    done
done

date +"%Y/%m/%d %H:%M:%S.%N" >> measurements.txt