
# Table of Contents

1.  [Exercise 1: Generating a crystal structure](#orgebf38a1)
    1.  [Part A: Generate an FCC structure in the cubic (conventional) cell:](#org57ddc84)
    2.  [Part B: Generate an FCC structure in the primitive cell:](#org3e05ce3)
    3.  [Part C: Create and visualize a 4x4x4 supercell](#org9f7e3c9)
    4.  [EXTRA-1.1: Performing the inverse operation:](#org22696e6)
    5.  [EXTRA-1.2: constraints on the lattice vectors:](#org98cab6d)
    6.  [EXTRA-1.3: Implementing '`cart_to_cryst`' and '`cryst_to_cart`':](#org79fa1e3)
    7.  [EXTRA-1.4: Implementing the Periodic Boundary Conditions:](#org81b106c)
2.  [Exercise 2: Determine the lattice constant of Si diamond crystal](#org436fe41)
    1.  [Part A: Generate Si crystal in the diamond configuration](#org5d58587)
    2.  [Part B: computing the lattice constant of Si](#org932329e)
3.  [md; melt](#orga0fe926)
    1.  [verlet/leapfrog, total energy, thormostat, damping etc.](#org4dcc41a)
    2.  [computing the force: empirical potentials](#orgb4391b7)
    3.  [LJ/melt in lammps examples](#orgc28c32f)
4.  [create surface](#orgaa645a9)
5.  [minimize surface](#orgd4b62da)
6.  [difference between md and minimization](#org5ff6562)
    1.  [try the melt example with adding damp at the end](#org6e2bc46)
    2.  [](#orge6b6cb7)


<a id="orgebf38a1"></a>

# Exercise 1: Generating a crystal structure


<a id="org57ddc84"></a>

## Part A: Generate an FCC structure in the cubic (conventional) cell:

Face-Centered Cubic (FCC) structure is characterised by one atom at the origin of the cubic cell,
and one atom at the center of each face of the cubic cell (thus 4 atoms per cell).
This is also called the conventional cell of FCC.

The lattice vectors are:

    a1 = [ 1.0, 0.0, 0.0 ]
    a2 = [ 0.0, 1.0, 0.0 ]
    a3 = [ 0.0, 0.0, 1.0 ]

and 4 atomic basis vectors; one at the origin, and one at the center of each face:

    Natoms = 4
    
    basis[0] = [ 0.0, 0.0, 0.0 ]
    basis[1] = [ 0.5, 0.5, 0.0 ]
    basis[2] = [ 0.5, 0.0, 0.5 ]
    basis[3] = [ 0.0, 0.5, 0.5 ]

1.  Modify the '`generate_cell.py`' script to generate the conventional cell of FCC structure, as given above.
    Call the script with python, and redirect the output to a file:
    
        python generate_cell.py > my_output_filename.xyz

2.  Use '`ovito`' program to visualize the xyz output, draw bonds between the atoms.


<a id="org3e05ce3"></a>

## Part B: Generate an FCC structure in the primitive cell:

The cubic cell is however not the primitive cell of FCC.
The primitive cell of FCC has rhombohedral lattice vectors with all angles equal at 60 degrees:

    a1 = [ -0.5, 0.0, 0.5 ]
    a2 = [  0.0, 0.5, 0.5 ]
    a3 = [ -0.5, 0.5, 0.0 ]

and one atom at the origin of the cell:

    Natoms = 1
    basis[0] = [ 0.0, 0.0, 0.0 ]

1.  Modify the '`generate_cell.py`' script to generate the primitive cell of FCC structure. Use '`ovito`' program to visualize the output.


<a id="org9f7e3c9"></a>

## Part C: Create and visualize a 4x4x4 supercell

1.  Modify the '`generate_cell.py`' script to generate 4x4x4 supercells of both structures generated above.

2.  Use '`ovito`' to visualize. The script will generate two structures in the output file:
    one for the unit cell, and one for the supercell. Navigate the structures in '`ovito`'.

3.  Test the structure of your supercell in ovito, using the modifier "Polyhedral Template Matching".
    You might need to tick the boxes "cubic diamond" or "hexagonal diamond".


<a id="org22696e6"></a>

## EXTRA-1.1: Performing the inverse operation:

In the script, the atomic positions are written in the crystal (fractional) coordinates,
which are then transformed by matrix multiplication into the cartesian coordinates for output.
How would you perform the inverse operation? That is, how to go from cartesian coordinates
back to the crystal coordinates, given the matrix of lattice vectors?


<a id="org98cab6d"></a>

## EXTRA-1.2: constraints on the lattice vectors:

In terms of their direction and norm, is there any constraint on the lattice vectors?

Hint: The determinant of the matrix of lattice vectors gives the absolute value of
the volume of the cell; which can also be computed by the scalar triple product: a1.(a2 x a3).

Physically, is it meaningful for the volume of a cell to be zero?
Can the lattice vectors be coplanar, or collinear?
What does it mean if a matrix has determinant equal to zero? What happens to the metric tensor aTa?


<a id="org79fa1e3"></a>

## EXTRA-1.3: Implementing '`cart_to_cryst`' and '`cryst_to_cart`':

Try to implement the transformations from cartesian to crystal, and from crystal to
cartesian coordinates, in your favourite programming language.


<a id="org81b106c"></a>

## EXTRA-1.4: Implementing the Periodic Boundary Conditions:

With the tools from previous question, it is straight-forward to implement the
Periodic Boundary Conditions (PBC).

The atomic positions expressed in crystal coordinates are limited to the range [0:1).
Based on this, we can write a loop to perform PBC on each vector of atomic positions:

    for i in range(Natoms):
       # transform from cartesian to crystal coordinates
       vec[i] = cart_to_cryst( vec[i], lattice )
       # check condition for position: in crystal coordinates should be on the range [0:1)
       vec[i] = periodic( vec[i] )
       # transform from crystal back to cartesian coordinates
       vec[i] = cryst_to_cart( vec[i], lattice )

Try to implement the loop sketched above in your favourite programming language.

PS: It is possible to compute distances between atoms in PBC using such routines.
However, then it's better to shift the range [0:1) to [-0.5:0.5).
Do you have any idea why that is? Try to sketch what is happening.


<a id="org436fe41"></a>

# Exercise 2: Determine the lattice constant of Si diamond crystal


<a id="org5d58587"></a>

## Part A: Generate Si crystal in the diamond configuration

Diamond lattice is composed of two FCC lattices shifted by 0.25.

The convetional cell is cubic, with lattice vectors:

    a1 = [ 1.0, 0.0, 0.0 ]
    a2 = [ 0.0, 1.0, 0.0 ]
    a3 = [ 0.0, 0.0, 1.0 ]

and two unit cells of FCC, shifted by 0.25:

    Natoms = 8
    basis[0] = [ 0.0, 0.0, 0.0 ]
    basis[1] = [ 0.5, 0.5, 0.0 ]
    basis[2] = [ 0.5, 0.0, 0.5 ]
    basis[3] = [ 0.0, 0.5, 0.5 ]
    basis[4] = [ 0.25, 0.25, 0.25 ]
    basis[5] = [ 0.75, 0.75, 0.25 ]
    basis[6] = [ 0.75, 0.25, 0.75 ]
    basis[7] = [ 0.25, 0.75, 0.75 ]

1.  Create the diamond structure in the conventional (cubic) cell, as given above.
    Verify the output structure by visualization.


<a id="org932329e"></a>

## Part B: computing the lattice constant of Si

We can compute the lattice constant of a crystal, by finding the value of \`a0\` which minimizes
the total energy of the cell.
In this exercise, we do the task manually, meaning we explicitly select some values of \`a0\`,
generate the structure with given \`a0\`, compute the corresponding energy, and plot the curve.
The curve can then be fitted with a polynomial to find the exact value of \`a0\` at the minimum.

If you feel confident in scripting, skip to EXTRA-2.1.

1.  Generate the diamond structure with several different values of the lattice constant \`a0\`.
2.  Compute the value of total energy for each \`a0\`.
    Use LAMMPS software (command \`lmp\`) to compute the energy:

To compute the total energy, use the script

compute the total energy for various values of the lattice parameter \`a0\`. Plot the enegy per atom vs value of \`a0\`, and find the minimum.


<a id="orga0fe926"></a>

# md; melt


<a id="org4dcc41a"></a>

## verlet/leapfrog, total energy, thormostat, damping etc.


<a id="orgb4391b7"></a>

## computing the force: empirical potentials

LJ


<a id="orgc28c32f"></a>

## LJ/melt in lammps examples


<a id="orgaa645a9"></a>

# create surface

create (by hand?) Si diamond in conventional cell (8atms per cell)


<a id="orgd4b62da"></a>

# minimize surface

surface as defect, minimize Si diamond with reaxff? maybe sw.


<a id="org5ff6562"></a>

# difference between md and minimization

Taking energy out of the system.
damped dynamics, quickmin/SD, CG, other &#x2026;


<a id="org6e2bc46"></a>

## try the melt example with adding damp at the end


<a id="orge6b6cb7"></a>

## 

