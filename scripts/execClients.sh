#!/usr/bin/bash

PORT=8000
IP=127.0.0.1

NUM_TRIALS=10

git log -1 | head -1 > results.txt
date +"%Y/%m/%d %H:%M:%S.%N" >> results.txt

q=()
n=()
m=()

N=()
M=()
tau=()
for((p=0;p<${#q[@]};p++)); do
    MM=(${M[$i]//:/ })
    TT=(${tau[$i]//:/ })

    for((j=0;j<${#MM[@]};j++)); do
    	for((nt=0;nt<${NUM_TRIALS};nt++)); do
    		echo "../LZKP -p0 -q${q[$p]} --ip ${IP} --port ${PORT}"
    		S=$(echo "../LZKP -p0 -q${q[$p]} --ip ${IP} --port ${PORT}") # Single threaded
    		echo ${S} | tee -a results.txt
    		wait
    		sleep 2
    	done

    	for((nt=0;nt<${NUM_TRIALS};nt++)); do
    		echo "../LZKP -p0 -q${q[$p]} --ip ${IP} --port ${PORT} -x"
    		S=$(echo "../LZKP -p0 -q${q[$p]} --ip ${IP} --port ${PORT} -x") # Single threaded
    		echo ${S} | tee -a results.txt
    		wait
    		sleep 2
    	done
    done
done


q=("15" "15" "31" "61")
n=("256" "256" "512" "1024")
m=("1024" "4096" "2048" "4096")

N=("2" "4" "8" "16")
M=("40:80" "20:40" "14:27" "10:20")
for((p=0;p<${#q[@]};p++)); do
    MM=(${M[$i]//:/ })

    for((j=0;j<${#MM[@]};j++)); do
    	for((nt=0;nt<${NUM_TRIALS};nt++)); do
    		echo "./LZKP -p1 -q${q[$p]} --ip ${IP} --port ${PORT}"
    		S=$(echo "./LZKP -p1 -q${q[$p]} --ip ${IP} --port ${PORT}") # Single threaded
    		echo ${S} | tee -a results.txt
    		wait
    		sleep 2
    	done

    	for((nt=0;nt<${NUM_TRIALS};nt++)); do
    		echo "./LZKP -p1 -q${q[$p]} --ip ${IP} --port ${PORT} -x"
    		S=$(echo "./LZKP -p1 -q${q[$p]} --ip ${IP} --port ${PORT} -x") # Multi threaded
    		echo ${S} | tee -a results.txt
    		wait
    		sleep 2
    	done
    done
done

date +"%Y/%m/%d %H:%M:%S.%N" >> results.txt