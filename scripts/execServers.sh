#!/usr/bin/bash

PORT=8000

NUM_TRIALS=10

git log -1 | head -1 > measurements.txt
date +"%Y/%m/%d %H:%M:%S.%N" | tee -a measurements.txt


# Parameters are going in parallel, first entry of each parameter are corresponding
q=()	# e.g q=("31")
n=()	# e.g n=("256")
m=()	# e.g m=("1024")

N=()	# e.g. N=("2")
M=()	# Split with : for multiple values e.g. M=("40:80")
tau=()
echo "protocol 1" | tee -a measurements.txt
for((p=0;p<${#q[@]};p++)); do
    for((i=0;i<${#N[@]};i++)); do
        MM=(${M[$i]//:/ })
        TT=(${tau[$i]//:/ })

        for((j=0;j<${#MM[@]};j++)); do
        	R=""
        	for((nt=0;nt<${NUM_TRIALS};nt++)); do
        		echo "../LZKP -p0 -i -q${q[$p]} -M${MM[$j]} -t${TT[$j]} -N${N[$i]} -n${n[$p]} -m${m[$p]} -a --port ${PORT}"
        		# S=$(echo "../LZKP -p0 -i -q${q[$p]} -M${MM[$j]} -t${TT[$j]} -N${N[$i]} -n${n[$p]} -m${m[$p]} -a --port ${PORT}") # Single threaded
        		#R="${R}${S} "
        		# wait
        	done
        	echo ${R} | tee -a measurements.txt

        	# R=""
        	# for((nt=0;nt<${NUM_TRIALS};nt++)); do
        		# echo "../LZKP -p0 -i -q${q[$p]} -M${MM[$j]} -t${TT[$j]} -N${N[$i]} -n${n[$p]} -m${m[$p]} -a --port ${PORT} -x8"
        		#S=$(echo "../LZKP -p0 -i -q${q[$p]} -M${MM[$j]} -t${TT[$j]} -N${N[$i]} -n${n[$p]} -m${m[$p]} -a --port ${PORT} -x8") # Multi threaded
        		#R="${R}${S} "
        		# wait
        	# done
        	# echo ${R} | tee -a measurements.txt
        done
    done
done

X=()
for((x=0;x<${#X[@]};x++)); do
    R=""
    for((nt=0;nt<${NUM_TRIALS};nt++)); do
        echo "../LZKP -p0 -i -q61 -M45 -t28 -N64 -n1024 -m4096 -a --port ${PORT} -x${X[$x]}" # Multi threaded
        #S=$(../LZKP -p0 -i -q61 -M45 -t28 -N64 -n1024 -m4096 -a --port ${PORT} -x${X[$x]})
        #R="${R}${S} "
        #wait
    done
    echo ${R} | tee -a measurements.txt

    R=""
    for((nt=0;nt<${NUM_TRIALS};nt++)); do
        echo "../LZKP -p0 -i -q61 -M85 -t49 -N64 -n1024 -m4096 -a --port ${PORT} -x${X[$x]}" # Multi threaded
        #S=$(../LZKP -p0 -i -q61 -M85 -t49 -N64 -n1024 -m4096 -a --port ${PORT} -x${X[$x]})
        #R="${R}${S} "
        #wait
    done
    echo ${R} | tee -a measurements.txt
done


q=()
n=()
m=()

N=()
M=()
echo "protocol 2" | tee -a measurements.txt
for((p=0;p<${#q[@]};p++)); do
    for((i=0;i<${#N[@]};i++)); do
        MM=(${M[$i]//:/ })

        for((j=0;j<${#MM[@]};j++)); do
        	R=""
        	for((nt=0;nt<${NUM_TRIALS};nt++)); do
        		echo "../LZKP -p1 -i -q${q[$p]} -M${MM[$j]} -t1 -N${N[$i]} -n${n[$p]} -m${m[$p]} -a --port ${PORT}"
        		# S=$(echo "../LZKP -p1 -i -q${q[$p]} -M${MM[$j]} -t1 -N${N[$i]} -n${n[$p]} -m${m[$p]} -a --port ${PORT}") # Single threaded
        		#R="${R}${S} "
        		wait
        	done
        	echo ${R} | tee -a measurements.txt

        	# R=""
        	# for((nt=0;nt<${NUM_TRIALS};nt++)); do
        		# echo "../LZKP -p1 -i -q${q[$p]} -M${MM[$j]} -t1 -N${N[$i]} -n${n[$p]} -m${m[$p]} -a --port ${PORT} -x8"
        		#S=$(echo "../LZKP -p1 -i -q${q[$p]} -M${MM[$j]} -t1 -N${N[$i]} -n${n[$p]} -m${m[$p]} -a --port ${PORT} -x8") # Multi threaded
        		#R="${R}${S} "
        		# wait
        	# done
        	# echo ${R} >> measurements.txt
        done
    done
done

X=()
for((x=0;x<${#X[@]};x++)); do
    R=""
    for((nt=0;nt<${NUM_TRIALS};nt++)); do
        echo "../LZKP -p1 -i -q61 -M7 -t1 -N64 -n1024 -m4096 -a --port ${PORT} -x${X[$x]}" # Multi threaded
        #S=$(../LZKP -p1 -i -q61 -M7 -t1 -N64 -n1024 -m4096 -a --port ${PORT} -x${X[$x]})
        #R="${R}${S} "
        #wait
    done
    echo ${R} | tee -a measurements.txt

    R=""
    for((nt=0;nt<${NUM_TRIALS};nt++)); do
        echo "../LZKP -p1 -i -q61 -M14 -t1 -N64 -n1024 -m4096 -a --port ${PORT} -x${X[$x]}" # Multi threaded
        #S=$(../LZKP -p1 -i -q61 -M14 -t1 -N64 -n1024 -m4096 -a --port ${PORT} -x${X[$x]})
        #R="${R}${S} "
        #wait
    done
    echo ${R} | tee -a measurements.txt
done

date +"%Y/%m/%d %H:%M:%S.%N" | tee -a measurements.txt