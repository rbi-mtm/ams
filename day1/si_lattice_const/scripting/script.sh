a_start=5.2
a_end=5.8

### how many incremental steps to make
n_step=10

### compute delta value to increment each step:
### delta = ( a_end - a_start ) / n_step
delta=`python -c "print( (${a_end}-${a_start}) / ${n_step} )"`


### loop for istep from 0 to n_step
for istep in `seq 0 ${n_step}`
do
	  ### compute value of a0 for this istep:
	  ### a0 = a_start + (istep * delta)
	  a0=`python -c "print( ${a_start} + ${istep}*${delta})"`

	  ### generate diamond cell with this a0, write xyz to file
	  python gen_diamond.py ${a0} > struc.xyz

	  ### convert xyz format to lammps data format
	  ### (make sure lammps.in in the next step reads from struc.lmp)
	  python ../xyz2lmp.py struc.xyz > struc.lmp

	  ### launch lmp with input file lammps.in, save output to file
	  lmp -log none -in lammps.in > out.lmp

	  ### grep for total energy in the lmp output
	  etot=`grep "TotEng" out.lmp -A 1 | tail -n 1`

	  ## output to screen
    #	echo ${a0} ${etot}
	  printf "%9.6f\t" ${a0}
	  printf "%9.6f\n" ${etot}

done

