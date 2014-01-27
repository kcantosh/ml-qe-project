#!/bin/bash
#SBATCH -J batch_ti
#SBATCH -o batch_ti.o%j
#SBATCH -N 1 -n 8
#SBATCH -t 24:00:00
#SBATCH -p normal

cd $SLURM_SUBMIT_DIR

#a complete workflow
#note the use of bash commands to handle intermediate
#steps eg., parsing stream output, bc calculator etc

#first unit
#anatase; ibrav==7
aa=3.784
ca=9.515
Oxa=0.0
Oya=0.0
Oza=0.208
ibrav=7
a=3.784
c=9.515
Ox=0.0
Oy=0.0
Oz=0.208




#last unit cell, or thereabouts
#rutile; ibrav==6
#a=4.603 	(step += 0.017)
#c=2.966 	(step -= 0.134)
#Ox=0.3045 	(step += 0.006)
#Oy=0.3045	(step += 0.006)
#Oz=0.0		(step -= 0.004)


#nmr magic angle spinning spectral simulation parameters
inpa="25.0 2.5 3 1.5833 50750000 20 1 128 1  40 78.125 00 -1700 1 1"
inpb="0.0 00 100 3000000 200000 0.0 10 10 9000000.0"


#create some files with comment fields at the top
echo "# Cq eta a c Ox Oy Oz" > anatase_rutile_runs500.txt
echo "# training data for Ox" > Ox_svm500.txt
echo "# training data for Oy" > Oy_svm500.txt
echo "# training data for Oz" > Oz_svm500.txt
echo "# training data for a" > a_svm500.txt
echo "# training data for c" > c_svm500.txt

echo " " 
echo "Starting job on `hostname` at `date`" 
echo " " 

#preform a string of processes in batch:
#	write a qe input file
#	perform scf
#	calculate nmr parameters
#	simulate nmr linshape
#	write features/values to text files
#	iterate unit cell parameters for next iteration

for (( i=1; i<=1000; i++ ))
do
	
	#set the symmetry according to unit cell type
	if [ "$(echo "$a > $c" | bc)" -eq "1" ] 
	then
		ibrav=6 
	fi

	#update the input file
	sed 's/@A@/'$a'/g' anatase_batch.in |
	sed 's/@C@/'$c'/g' |
	sed 's/@IBRAV@/'$ibrav'/g'|
	sed 's/@Ox@/'$Ox'/g' |
	sed 's/@Oy@/'$Oy'/g' |
	sed 's/@Oz@/'$Oz'/g' > anatase_tmp.in

	#scf calculation 
	ibrun ../qe_latest/espresso-4.3/bin/pw.x < anatase_tmp.in > anatase_tmp.out

	#nmr calculation
	efg_vals=$(ibrun ../qe_latest/espresso-4.3/bin/gipaw.x < nmr_anatase_tmp.in |awk '/Ti/&&/Cq/{print $8,$11}')
	#ibrun ../qe_latest/espresso-4.3/bin/gipaw.x < nmr_anatase_tmp.in > foo.out

	#update parameter file
	echo $efg_vals $ibrav $a $c $Ox $Oy $Oz >> anatase_rutile_runs500.txt
	
	#simulate magic angle spinning nmr lineshape for this structure
	tdata=$(echo "$inpa $efg_vals $inpb" | ./sim_mas.x)
	
	#save svm training data
	echo $a $tdata >> a_svm500.txt
	echo $c $tdata >> c_svm500.txt
	echo $Ox $tdata >> Ox_svm500.txt
	echo $Oy $tdata >> Oy_svm500.txt
	echo $Oz $tdata >> Oz_svm500.txt

	rn=$(echo "scale=10; ($RANDOM / 32767)" | bc)

	#calculate next lattice and cell parameters
	a=$(echo "scale=10; $aa + $rn*0.82" | bc)
	c=$(echo "scale=10; $ca - $rn*6.55" | bc)
	Ox=$(echo "scale=10; $Oxa - $rn*0.3" | bc)
	Oy=$(echo "scale=10; $Oya - $rn*0.3" | bc)
	Oz=$(echo "scale=10; $Oza + $rn*0.2" | bc)


done

echo " " 
echo "Completing job on `hostname` at `date`" 
echo " " 
