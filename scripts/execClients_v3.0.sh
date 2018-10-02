#!/usr/bin/bash

PORT=8000
IP=127.0.0.1

NUM_TRIALS=20

git log -1 | head -1 > results.txt
date +"%Y/%m/%d %H:%M:%S.%N" >> results.txt

# Best to copy parameters section from execServers config file
q=("15" "15" "31" "59" "61")
n=("256" "256" "512" "1024" "1024")
m=("1024" "4096" "2048" "4096" "4096")

N=("2" "4" "8" "16" "32" "64")
M=("75:145" "55:105" "55:95" "45:95" "45:85" "45:85")
tau=("34:63" "32:57" "38:57" "26:63" "28:47" "28:49")
for((p=0;p<${#q[@]};p++)); do
    for((i=0;i<${#N[@]};i++)); do
        MM=(${M[$i]//:/ })
        TT=(${tau[$i]//:/ })

        for((j=0;j<${#MM[@]};j++)); do
        	for((nt=0;nt<${NUM_TRIALS};nt++)); do
        		echo "../LZKP -p0 -q${q[$p]} --ip ${IP} --port ${PORT}"
        		S=$(../LZKP -p0 -q${q[$p]} --ip ${IP} --port ${PORT}) # Single threaded
        		echo ${S} | tee -a results.txt
        		wait
        		sleep 1
        	done

        	# for((nt=0;nt<${NUM_TRIALS};nt++)); do
        		# echo "../LZKP -p0 -q${q[$p]} --ip ${IP} --port ${PORT} -x8"
        		# S=$(echo "../LZKP -p0 -q${q[$p]} --ip ${IP} --port ${PORT} -x8") # Single threaded
        		# echo ${S} | tee -a results.txt
        		# wait
        		# sleep 1
        	# done
        done
    done
done

X=("2" "4" "8" "16" "32")
for((x=0;x<${#X[@]};x++)); do
    R=""
    for((nt=0;nt<${NUM_TRIALS};nt++)); do
        echo "../LZKP -p0 -q61 --ip ${IP} --port ${PORT} -x${X[$x]}"
        S=$(../LZKP -p0 -q61 --ip ${IP} --port ${PORT} -x${X[$x]}) # Multi threaded
        echo ${S} | tee -a results.txt
        wait
        sleep 1
    done

    R=""
    for((nt=0;nt<${NUM_TRIALS};nt++)); do
        echo "../LZKP -p0 -q61 --ip ${IP} --port ${PORT} -x${X[$x]}"
        S=$(../LZKP -p0 -q61 --ip ${IP} --port ${PORT} -x${X[$x]}) # Multi threaded
        echo ${S} | tee -a results.txt
        wait
        sleep 1
    done
done


q=("15" "15" "31" "59" "61")
n=("256" "256" "512" "1024" "1024")
m=("1024" "4096" "2048" "4096" "4096")

N=("2" "4" "8" "16" "32" "64")
M=("40:80" "20:40" "14:27" "10:20" "8:16" "7:14")
for((p=0;p<${#q[@]};p++)); do
    for((i=0;i<${#N[@]};i++)); do
        MM=(${M[$i]//:/ })

        for((j=0;j<${#MM[@]};j++)); do
        	for((nt=0;nt<${NUM_TRIALS};nt++)); do
        		echo "../LZKP -p1 -q${q[$p]} --ip ${IP} --port ${PORT}"
        		S=$(../LZKP -p1 -q${q[$p]} --ip ${IP} --port ${PORT}) # Single threaded
        		echo ${S} | tee -a results.txt
        		wait
        		sleep 1
        	done

        	# for((nt=0;nt<${NUM_TRIALS};nt++)); do
        		# echo "../LZKP -p1 -q${q[$p]} --ip ${IP} --port ${PORT} -x8"
        		# S=$(echo "../LZKP -p1 -q${q[$p]} --ip ${IP} --port ${PORT} -x8") # Multi threaded
        		# echo ${S} | tee -a results.txt
        		# wait
        		# sleep 1
        	# done
        done
    done
done

X=("2" "4" "8" "16" "32")
for((x=0;x<${#X[@]};x++)); do
    R=""
    for((nt=0;nt<${NUM_TRIALS};nt++)); do
        echo "../LZKP -p1 -q61 --ip ${IP} --port ${PORT} -x${X[$x]}"
        S=$(../LZKP -p1 -q61 --ip ${IP} --port ${PORT} -x${X[$x]}) # Multi threaded
        echo ${S} | tee -a results.txt
        wait
        sleep 1
    done

    R=""
    for((nt=0;nt<${NUM_TRIALS};nt++)); do
        echo "../LZKP -p1 -q61 --ip ${IP} --port ${PORT} -x${X[$x]}"
        S=$(../LZKP -p1 -q61 --ip ${IP} --port ${PORT} -x${X[$x]}) # Multi threaded
        echo ${S} | tee -a results.txt
        wait
        sleep 1
    done
done

date +"%Y/%m/%d %H:%M:%S.%N" >> results.txt