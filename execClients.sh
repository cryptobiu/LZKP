#!/usr/bin/bash

q=("15" "15" "31" "61")
n=("256" "256" "512" "1024")
m=("1024" "4096" "2048" "4096")

N=("2" "4" "8" "16")
M=("40:80" "20:40" "14:21" "10:20")


PORT=8000
IP=127.0.0.1

NUM_TRIALS=10

touch results.txt
for((p=0;p<${#q[@]};p++)); do
	for((i=0;i<${#N[@]};i++)); do
	    #echo "$i: ${N[$i]}"

	    MM=(${M[$i]//:/ })
	    for((j=0;j<${#MM[@]};j++)); do
	    	#echo "$j: ${MM[$j]}"
	    	for((nt=0;nt<${NUM_TRIALS};nt++)); do
	    		S=$(echo "./LZKP -p1 -q${q[$p]} --ip ${IP} --port ${PORT}") # Single threaded
	    		echo "echo ${S} >> results.txt"
	    		sleep 3
	    	done

	    	for((nt=0;nt<${NUM_TRIALS};nt++)); do
	    		S=$(echo "./LZKP -p1 -q${q[$p]} --ip ${IP} --port ${PORT} -x") # Multi threaded
	    		echo "echo ${S} >> results.txt"
	    		sleep 3
	    	done
	    done
	done
done