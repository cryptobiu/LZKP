#!/usr/bin/bash

PORT=8000
IP=127.0.0.1

NUM_TRIALS=10

git log -1 | head -1 > results.txt
date +"%Y/%m/%d %H:%M:%S.%N" >> results.txt

# Best to copy parameters section from execServers config file
q=()
n=()
m=()

N=()
M=()
tau=()
for((p=0;p<${#q[@]};p++)); do
    for((i=0;i<${#N[@]};i++)); do
        MM=(${M[$i]//:/ })
        TT=(${tau[$i]//:/ })

        for((j=0;j<${#MM[@]};j++)); do
        	for((nt=0;nt<${NUM_TRIALS};nt++)); do
        		echo "../LZKP -p0 -q${q[$p]} --ip ${IP} --port ${PORT}"
        		# S=$(echo "../LZKP -p0 -q${q[$p]} --ip ${IP} --port ${PORT}") # Single threaded
        		# echo ${S} | tee -a results.txt
        		# wait
        		# sleep 1
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

X=()
for((x=0;x<${#X[@]};x++)); do
    R=""
    for((nt=0;nt<${NUM_TRIALS};nt++)); do
        echo "../LZKP -p0 -q61 --ip ${IP} --port ${PORT} -x${X[$x]}"
        #S=$(../LZKP -p0 -q61 --ip ${IP} --port ${PORT} -x${X[$x]}) # Multi threaded
        #echo ${S} | tee -a results.txt
        #wait
        #sleep 1
    done

    R=""
    for((nt=0;nt<${NUM_TRIALS};nt++)); do
        echo "../LZKP -p0 -q61 --ip ${IP} --port ${PORT} -x${X[$x]}"
        #S=$(../LZKP -p0 -q61 --ip ${IP} --port ${PORT} -x${X[$x]}) # Multi threaded
        #echo ${S} | tee -a results.txt
        #wait
        #sleep 1
    done
done


q=()
n=()
m=()

N=()
M=()
for((p=0;p<${#q[@]};p++)); do
    for((i=0;i<${#N[@]};i++)); do
        MM=(${M[$i]//:/ })

        for((j=0;j<${#MM[@]};j++)); do
        	for((nt=0;nt<${NUM_TRIALS};nt++)); do
        		echo "../LZKP -p1 -q${q[$p]} --ip ${IP} --port ${PORT}"
        		S=$(echo "../LZKP -p1 -q${q[$p]} --ip ${IP} --port ${PORT}") # Single threaded
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

X=()
for((x=0;x<${#X[@]};x++)); do
    R=""
    for((nt=0;nt<${NUM_TRIALS};nt++)); do
        echo "../LZKP -p1 -q61 --ip ${IP} --port ${PORT} -x${X[$x]}"
        #S=$(../LZKP -p1 -q61 --ip ${IP} --port ${PORT} -x${X[$x]}) # Multi threaded
        #echo ${S} | tee -a results.txt
        #wait
        #sleep 1
    done

    R=""
    for((nt=0;nt<${NUM_TRIALS};nt++)); do
        echo "../LZKP -p1 -q61 --ip ${IP} --port ${PORT} -x${X[$x]}"
        #S=$(../LZKP -p1 -q61 --ip ${IP} --port ${PORT} -x${X[$x]}) # Multi threaded
        #echo ${S} | tee -a results.txt
        #wait
        #sleep 1
    done
done

date +"%Y/%m/%d %H:%M:%S.%N" >> results.txt