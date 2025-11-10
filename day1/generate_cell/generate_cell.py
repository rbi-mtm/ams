import numpy as np

### lattice constant
a0 = 1.82

### lattice vectors in real space: a_i = [x, y, z]
a1 = [ 1.0, 0.0, 0.0 ]
a2 = [ 0.0, 1.0, 0.0 ]
a3 = [ 0.0, 0.0, 1.0 ]

### put vectors into matrix, multiply by a0 (rescale)
lattice = np.asarray([a1, a2, a3]) * a0


### number of atoms in the unit cell (number of atomic basis vectors)
Natoms = 2

### define an array of size[Natoms,3] to store the vectors of atomic bases
basis = np.ndarray([Natoms,3], dtype=float)

### Position of each atom in the unit cell (atomic basis vectors), in fractional/crystal coordinates.
### Values between 0.0 and 1.0
basis[0] = [ 0.0, 0.0, 0.0 ]
basis[1] = [ 0.5, 0.5, 0.5 ]



### write the unit cell in xyz-format
print(Natoms)
print('Lattice="',*lattice.T.flatten(),'" properties=species:I:1:pos:R:3')

for n in range(Natoms):
    # compute the real-space coordinates:
    #   matrix-vector multiplication of lattice vectors
    #   with the basis vector in fractional coords (change of basis operation)
    r = np.matmul( lattice, basis[n] )
    print( n, *r )



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


## how many replicas to create in each direction?
n_replicate_x = 2
n_replicate_y = 3
n_replicate_z = 1

# Rescale the number of atoms in the supercell (for output only)
replicated_Natoms = Natoms * n_replicate_x * n_replicate_y * n_replicate_z

# Rescale the lattice (for output only)
replicated_lat=lattice.copy()
replicated_lat[0] *= n_replicate_x
replicated_lat[1] *= n_replicate_y
replicated_lat[2] *= n_replicate_z


# Output the supercell
print(replicated_Natoms)
print('Lattice="',*replicated_lat.T.flatten(),'" properties=species:I:1:pos:R:3')

# loop over the replications
for n1 in range(n_replicate_x):
    for n2 in range(n_replicate_y):
        for n3 in range(n_replicate_z):
            #
            # construct the `shift` vector for the current box (in fractional coordinates)
            shift = [ n1, n2, n3 ]
            #
            # loop over atoms in the basis:
            for n in range(Natoms):
                # shift the atomic basis vector
                bas = basis[n]+shift
                # obtain real space coords, from original lattice and shifted basis
                r = np.matmul( lattice, bas )
                print( n, *r )
