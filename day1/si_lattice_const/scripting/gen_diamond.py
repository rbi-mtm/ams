import numpy as np
import sys

if len(sys.argv) != 2:
    print("********* Specify input value of a0 as:")
    print(" python "+str(sys.argv[0])+" <a0>")
    print("********* Exiting.")
    sys.exit(-1)

### lattice constant
a0 = np.float64( sys.argv[1] )

### lattice vectors in real space: a_i = [x, y, z]
a1 = [ 1.0, 0.0, 0.0 ]
a2 = [ 0.0, 1.0, 0.0 ]
a3 = [ 0.0, 0.0, 1.0 ]

### put vectors into matrix, multiply by a0 (rescale)
lattice = np.asarray([a1, a2, a3]) * a0


### number of atoms in the unit cell (number of atomic basis vectors)
Natoms = 8

### define an array of size[Natoms,3] to store the vectors of atomic bases
basis = np.ndarray([Natoms,3], dtype=float)

### Position of each atom in the unit cell (atomic basis vectors),
### in fractional/crystal coordinates.
### Values between 0.0 and 1.0
basis[0] = [ 0.0, 0.0, 0.0 ]
basis[1] = [ 0.5, 0.5, 0.0 ]
basis[2] = [ 0.5, 0.0, 0.5 ]
basis[3] = [ 0.0, 0.5, 0.5 ]
basis[4] = [ 0.25, 0.25, 0.25 ]
basis[5] = [ 0.75, 0.75, 0.25 ]
basis[6] = [ 0.75, 0.25, 0.75 ]
basis[7] = [ 0.25, 0.75, 0.75 ]


### define an array of size[Natoms] to store the
### atomic type of each atom in the unit cell.
atm_type = np.ndarray([Natoms], dtype=object)
atm_type[:] = "Si"

### write the unit cell in xyz-format
print(Natoms)
print('Lattice="',*lattice.flatten(),'" properties=species:I:1:pos:R:3')

for n in range(Natoms):
    # compute the real-space coordinates:
    #   matrix-vector multiplication of lattice vectors
    #   with the basis vector in fractional coords (change of basis operation)
    r = np.matmul( lattice.T, basis[n] )
    print( atm_type[n], *r )


