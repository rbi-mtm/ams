import ase.io
import sys

if len(sys.argv) != 2:
    print( "**************************************")
    print( "  Incorrect number of arguments!" )
    print( "*********** script usage: ************")
    print( " python  xyz2lmp.py  input_filename  ")
    print( "**************************************")
    sys.exit(-1)

atms = ase.io.read( sys.argv[1], format="extxyz" )
## write to stdout
ase.io.write( "-", atms, format="lammps-data", atom_style="atomic" )
