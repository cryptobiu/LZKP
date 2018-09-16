#!/usr/bin/bash

PORT=8000

NUM_TRIALS=10

git log -1 | head -1 > measurements.txt
date +"%Y/%m/%d %H:%M:%S.%N" >> measurements.txt


# Parameters are going in parallel, first entry of each parameter are corresponding
q=("15" "15" "31" "61")
n=("256" "256" "512" "1024")
m=("1024" "4096" "2048" "4096")

N=("2" "4" "8" "16" "32" "64")
M=("75:145" "55:105" "55:95" "45:95" "45:85" "45:85")
tau=("34:63" "32:57" "38:57" "26:63" "28:47" "28:49")
echo "protocol 1" >> measurements.txt
for((p=0;p<${#q[@]};p++)); do
    for((i=0;i<${#N[@]};i++)); do
        MM=(${M[$i]//:/ })
        TT=(${tau[$i]//:/ })

        for((j=0;j<${#MM[@]};j++)); do
        	R=""
        	for((nt=0;nt<${NUM_TRIALS};nt++)); do
        		echo "../LZKP -p0 -i -q${q[$p]} -M${MM[$j]} -t${TT[$j]} -N${N[$i]} -n${n[$p]} -m${m[$p]} -a --port ${PORT}"
        		S=$(../LZKP -p0 -i -q${q[$p]} -M${MM[$j]} -t${TT[$j]} -N${N[$i]} -n${n[$p]} -m${m[$p]} -a --port ${PORT}) # Single threaded
        		R="${R}${S} "
        		wait
        	done
        	echo ${R} >> measurements.txt

        	# R=""
        	# for((nt=0;nt<${NUM_TRIALS};nt++)); do
        		# echo "../LZKP -p0 -i -q${q[$p]} -M${MM[$j]} -t${TT[$j]} -N${N[$i]} -n${n[$p]} -m${m[$p]} -a --port ${PORT} -x8"
        		#S=$(echo "../LZKP -p0 -i -q${q[$p]} -M${MM[$j]} -t${TT[$j]} -N${N[$i]} -n${n[$p]} -m${m[$p]} -a --port ${PORT} -x8") # Multi threaded
        		#R="${R}${S} "
        		# wait
        	# done
        	# echo ${R} >> measurements.txt
        done
    done
done

X=("2" "4" "8" "16" "32")
for((x=0;x<${#X[@]};x++)); do
    R=""
    for((nt=0;nt<${NUM_TRIALS};nt++)); do
        echo "../LZKP -p0 -i -q61 -M45 -t28 -N64 -n1024 -m4096 -a --port ${PORT} -x${X[$x]}" # Multi threaded
        S=$(../LZKP -p0 -i -q61 -M45 -t28 -N64 -n1024 -m4096 -a --port ${PORT} -x${X[$x]})
        R="${R}${S} "
        wait
    done
    echo ${R} >> measurements.txt

    R=""
    for((nt=0;nt<${NUM_TRIALS};nt++)); do
        echo "../LZKP -p0 -i -q61 -M85 -t49 -N64 -n1024 -m4096 -a --port ${PORT} -x${X[$x]}" # Multi threaded
        S=$(../LZKP -p0 -i -q61 -M85 -t49 -N64 -n1024 -m4096 -a --port ${PORT} -x${X[$x]})
        R="${R}${S} "
        wait
    done
    echo ${R} >> measurements.txt
done


q=("15" "15" "31" "61")
n=("256" "256" "512" "1024")
m=("1024" "4096" "2048" "4096")

N=("2" "4" "8" "16" "32" "64")
M=("40:80" "20:40" "14:27" "10:20" "8:16" "7:14")
echo "protocol 2" >> measurements.txt
for((p=0;p<${#q[@]};p++)); do
    for((i=0;i<${#N[@]};i++)); do
        MM=(${M[$i]//:/ })

        for((j=0;j<${#MM[@]};j++)); do
        	R=""
        	for((nt=0;nt<${NUM_TRIALS};nt++)); do
        		echo "../LZKP -p1 -i -q${q[$p]} -M${MM[$j]} -t1 -N${N[$i]} -n${n[$p]} -m${m[$p]} -a --port ${PORT}"
        		S=$(../LZKP -p1 -i -q${q[$p]} -M${MM[$j]} -t1 -N${N[$i]} -n${n[$p]} -m${m[$p]} -a --port ${PORT}) # Single threaded
        		R="${R}${S} "
        		wait
        	done
        	echo ${R} >> measurements.txt

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

X=("2" "4" "8" "16" "32")
for((x=0;x<${#X[@]};x++)); do
    R=""
    for((nt=0;nt<${NUM_TRIALS};nt++)); do
        echo "../LZKP -p1 -i -q61 -M7 -t1 -N64 -n1024 -m4096 -a --port ${PORT} -x${X[$x]}" # Multi threaded
        S=$(../LZKP -p1 -i -q61 -M7 -t1 -N64 -n1024 -m4096 -a --port ${PORT} -x${X[$x]})
        R="${R}${S} "
        wait
    done
    echo ${R} >> measurements.txt

    R=""
    for((nt=0;nt<${NUM_TRIALS};nt++)); do
        echo "../LZKP -p1 -i -q61 -M14 -t1 -N64 -n1024 -m4096 -a --port ${PORT} -x${X[$x]}" # Multi threaded
        S=$(../LZKP -p1 -i -q61 -M14 -t1 -N64 -n1024 -m4096 -a --port ${PORT} -x${X[$x]})
        R="${R}${S} "
        wait
    done
    echo ${R} >> measurements.txt
done

date +"%Y/%m/%d %H:%M:%S.%N" >> measurements.txt