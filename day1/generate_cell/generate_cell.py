import numpy as np

### Lattice constant
a0 = 1.82

### Define lattice vectors in real space: a_i = [x, y, z]
### Scaled by a0
a1 = np.array([ 1.0, 0.0, 0.0 ]) * a0
a2 = np.array([ 0.0, 1.0, 0.0 ]) * a0
a3 = np.array([ 0.0, 0.0, 1.0 ]) * a0

### put vectors into matrix columns
lattice = np.array([a1, a2, a3]).T


### number of atoms in the unit cell (number of atomic basis vectors)
Natoms = 2

### define an array of size[Natoms,3] to store the vectors of atomic bases
basis = np.ndarray([Natoms,3], dtype=float)
### define an array of size[Natoms] to store the atomic type of each atom in the unit cell
atm_type = np.ndarray([Natoms], dtype=object)
atm_type[:] = 1


### Position of each atom in the unit cell (atomic basis vectors),
### in fractional/crystal coordinates.
### Values between 0.0 and 1.0
basis[0] = np.array([ 0.0, 0.0, 0.0 ])
basis[1] = np.array([ 0.5, 0.5, 0.5 ])

### Atomic type of atoms in the basis
atm_type[0] = 1
# atm_type[1] = 2




### Write the unit cell in xyz-format:
### ==================================

### The number of atoms
print( Natoms )

### The simulation box: equal to the lattice vectors
print( 'Lattice="',*a1, *a2, *a3, '" properties=species:I:1:pos:R:3' )

### The atomic positions
for i in range(Natoms):
    ### compute the real-space coordinates:
    ###   matrix-vector multiplication of lattice vectors
    ###   with the basis vector in fractional coords (change of basis operation)
    r = np.matmul( lattice, basis[i] )
    ### output atomic type, and the position vector r
    print( atm_type[i], *r )





### ================================
###
### Create a supercell:
### Replicate the unit cell N times along each direction.
###
### Reminder:
###  Lattice repetition (symmetry operation):
###
###    r = n1*a1 + n2*a2 + n3*a3  ; n_i integers
###
### Strategy to create a supercell:
###  * define how many replications to make in each direction (n_i integers);
###  * loop over replications of the box;
###  * create a `shift` vector for each replication of the box, in fractional coordinates;
###  * apply the shift to each atomic basis vector;
###  * obtain the real-space vector by matrix multiplication.


### How many replicas to create in each direction?
n_replicate_x = 0
n_replicate_y = 0
n_replicate_z = 0

if n_replicate_x > 0 and n_replicate_y > 0 and n_replicate_z > 0:

    ### Output the supercell:

    ### Number of atoms is multiplied by number of replications
    print( Natoms * n_replicate_x * n_replicate_y * n_replicate_z )

    ### The simulation box is also multiplied in each direction
    print('Lattice="', \
          *a1 * n_replicate_x, \
          *a2 * n_replicate_y, \
          *a3 * n_replicate_z, \
          '" properties=species:I:1:pos:R:3')

    ### loop over the replications
    for n1 in range(n_replicate_x):
        for n2 in range(n_replicate_y):
            for n3 in range(n_replicate_z):
                #
                ### construct the `shift` vector for the current box (in fractional coordinates)
                shift = np.array([ n1, n2, n3 ])
                #
                ### loop over atoms in the basis:
                for i in range(Natoms):
                    ### shift the atomic basis vector
                    bas = basis[i] + shift
                    ### obtain real space coords, from original lattice and shifted basis
                    r = np.matmul( lattice, bas )
                    ### output
                    print( atm_type[i], *r )
