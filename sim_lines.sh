#!/bin/bash
inpa="25.0 2.5 3 1.5833 50750000 20 1 512 1  10 78.125 00 -1700 1 1"
inpb="0.0 00 100 3000000 200000 0.0 10 10 9000000.0"
echo "$inpa 5.0 0 $inpb" | ./sim_mas.x | sed 's/\ /\n/g' | awk 'BEGIN {FS=":"} {print $2}' > cq5.txt
echo "$inpa 1.0 0.9 $inpb" | ./sim_mas.x | sed 's/\ /\n/g' | awk 'BEGIN {FS=":"} {print $2}' > cq1.txt
