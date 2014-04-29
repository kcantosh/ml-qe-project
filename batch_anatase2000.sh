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

#using a minimum separation distance of 1.26(O2-) + 0.74(Ti4+)  ~ 2 A


min_dist_sq=4

#nmr magic angle spinning spectral simulation parameters
inpa="25.0 2.5 3 1.5833 50750000 20 1 512 1  10 78.125 00 -1700 1 1"
inpb="0.0 00 100 3000000 200000 0.0 10 10 9000000.0"


#create some files with comment fields at the top
echo "# Cq eta a c Ox Oy Oz" > anatase_rutile_runs2000.txt
echo "# training data for Ox" > Ox_svm2000.txt
echo "# training data for Oy" > Oy_svm2000.txt
echo "# training data for Oz" > Oz_svm2000.txt
echo "# training data for a" > a_svm2000.txt
echo "# training data for c" > c_svm2000.txt

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

for (( i=1; i<=2000; i++ ))
do
	

	dist=$(echo "scale=10; $a*$a*$Ox*$Ox + $a*$a*$Oy*$Oy + $c*$c*$Oz*$Oz" | bc)
	
	echo $sep

	#set the symmetry according to unit cell type
	if [ "$(echo "$a > $c" | bc)" -eq "1" ] 
	then
		ibrav=6 
	else
		ibrav=7
	fi


	if [ "$(echo "$dist > $min_dist_sq" | bc)" -eq "1" ]
	then
	


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
	echo $efg_vals $ibrav $a $c $Ox $Oy $Oz >> anatase_rutile_runs2000.txt
	
	#simulate magic angle spinning nmr lineshape for this structure
	tdata=$(echo "$inpa $efg_vals $inpb" | ./sim_mas.x)
	
	#save svm training data
	echo $a $tdata >> a_svm2000.txt
	echo $c $tdata >> c_svm2000.txt
	echo $Ox $tdata >> Ox_svm2000.txt
	echo $Oy $tdata >> Oy_svm2000.txt
	echo $Oz $tdata >> Oz_svm2000.txt

	fi

	rna=$(echo "scale=10; ($RANDOM / 32767)" | bc)
	rnb=$(echo "scale=10; ($RANDOM / 32767)" | bc)
	rnc=$(echo "scale=10; ($RANDOM / 32767)" | bc)
	rnd=$(echo "scale=10; ($RANDOM / 32767)" | bc)
	#rne=$(echo "scale=10; ($RANDOM / 32767)" | bc)

	#calculate next lattice and cell parameters
	a=$(echo "scale=10; $aa + $rna*0.82" | bc)
	c=$(echo "scale=10; $ca - $rnb*6.55" | bc)
	Ox=$(echo "scale=10; $Oxa - $rnc*0.3" | bc)
	Oy=$(echo "scale=10; $Oya - $rnc*0.3" | bc)
	Oz=$(echo "scale=10; $Oza + $rnd*0.2" | bc)


done

echo " " 
echo "Completing job on `hostname` at `date`" 
echo " " 
