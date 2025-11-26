- [Exercise 1: Generating a crystal structure](#sec-1)
  - [Part A: Generate an FCC structure in the cubic (conventional) cell:](#sec-1-1)
  - [Part B: Generate an FCC structure in the primitive cell:](#sec-1-2)
  - [Part C: Create and visualize a 4x4x4 supercell](#sec-1-3)
  - [EXTRA:](#sec-1-4)
- [Exercise 2: Determine the lattice constant of Si diamond crystal](#sec-2)
  - [Part A: Generate Si crystal in the diamond configuration](#sec-2-1)
  - [Part B: computing the lattice constant of Si](#sec-2-2)
  - [EXTRA:](#sec-2-3)
- [md; melt](#sec-3)
  - [verlet/leapfrog, total energy, thormostat, damping etc.](#sec-3-1)
  - [computing the force: empirical potentials](#sec-3-2)
  - [LJ/melt in lammps examples](#sec-3-3)
- [create surface](#sec-4)
- [minimize surface](#sec-5)
- [difference between md and minimization](#sec-6)
  - [try the melt example with adding damp at the end](#sec-6-1)
  - [](#sec-6-2)


# Exercise 1: Generating a crystal structure<a id="sec-1"></a>

## Part A: Generate an FCC structure in the cubic (conventional) cell:<a id="sec-1-1"></a>

Face-Centered Cubic (FCC) structure is characterised by one atom at the origin of the cubic cell, and one atom at the center of each face of the cubic cell (thus 4 atoms per cell). This is also called the conventional cell of FCC.

The lattice vectors are:

```python
a1 = [ 1.0, 0.0, 0.0 ]
a2 = [ 0.0, 1.0, 0.0 ]
a3 = [ 0.0, 0.0, 1.0 ]
```

and 4 atomic basis vectors; one at the origin, and one at the center of each face:

```python
Natoms = 4

basis[0] = [ 0.0, 0.0, 0.0 ]
basis[1] = [ 0.5, 0.5, 0.0 ]
basis[2] = [ 0.5, 0.0, 0.5 ]
basis[3] = [ 0.0, 0.5, 0.5 ]
```

1.  Modify the `generate_cell.py` script to generate the conventional cell of FCC structure, as given above. Call the script with python, and redirect the output to a file:
    
    ```bash
       python generate_cell.py > my_output_filename.xyz
    ```

2.  Use `ovito` program to visualize the xyz output, draw bonds between the atoms.

## Part B: Generate an FCC structure in the primitive cell:<a id="sec-1-2"></a>

The cubic cell is however not the primitive cell of FCC. The primitive cell of FCC has rhombohedral lattice vectors with all angles equal at 60 degrees:

```python
a1 = [ -0.5, 0.0, 0.5 ]
a2 = [  0.0, 0.5, 0.5 ]
a3 = [ -0.5, 0.5, 0.0 ]
```

and one atom at the origin of the cell:

```python
Natoms = 1
basis[0] = [ 0.0, 0.0, 0.0 ]
```

1.  Modify the `generate_cell.py` script to generate the primitive cell of FCC structure. Use `ovito` program to visualize the output.

## Part C: Create and visualize a 4x4x4 supercell<a id="sec-1-3"></a>

1.  Modify the `generate_cell.py` script to generate 4x4x4 supercells of both structures generated above.

2.  Use `ovito` to visualize. The script will generate two structures in the output file: one for the unit cell, and one for the supercell. Navigate the structures in `ovito`.

3.  Test the structure of your supercell in ovito, using the modifier "Polyhedral Template Matching". You might need to tick the boxes "cubic diamond" or "hexagonal diamond".

## EXTRA:<a id="sec-1-4"></a>

### EXTRA-1.1: Performing the inverse operation:<a id="sec-1-4-1"></a>

In the script, the atomic positions are written in the crystal (fractional) coordinates, which are then transformed by matrix multiplication into the cartesian coordinates for output. How would you perform the inverse operation? That is, how to go from cartesian coordinates back to the crystal coordinates, given the matrix of lattice vectors?

### EXTRA-1.2: constraints on the lattice vectors:<a id="sec-1-4-2"></a>

In terms of their direction and norm, is there any constraint on the lattice vectors?

Hint: The determinant of the matrix of lattice vectors gives the absolute value of the volume of the cell; which can also be computed by the scalar triple product: a1.(a2 x a3).

Physically, is it meaningful for the volume of a cell to be zero? Can the lattice vectors be coplanar, or collinear? What does it mean if a matrix has determinant equal to zero? What happens to the metric tensor aTa?

### EXTRA-1.3: Implementing `cart_to_cryst` and `cryst_to_cart`:<a id="sec-1-4-3"></a>

Try to implement the transformations from cartesian to crystal, and from crystal to cartesian coordinates, in your favourite programming language.

### EXTRA-1.4: Implementing the Periodic Boundary Conditions:<a id="sec-1-4-4"></a>

With the tools from previous question, it is straight-forward to implement the Periodic Boundary Conditions (PBC).

The atomic positions expressed in crystal coordinates are limited to the range [0:1). Based on this, we can write a loop to perform PBC on each vector of atomic positions:

```python
for i in range(Natoms):
    # transform from cartesian to crystal coordinates
    vec[i] = cart_to_cryst( vec[i], lattice )
    # check condition for position: in crystal coordinates should be on the range [0:1)
    vec[i] = periodic( vec[i] )
    # transform from crystal back to cartesian coordinates
    vec[i] = cryst_to_cart( vec[i], lattice )
```

Try to implement the loop sketched above in your favourite programming language.

PS: It is possible to compute distances between atoms in PBC using such routines. However, then it's better to shift the range [0:1) to [-0.5:0.5). Do you have any idea why that is? Try to sketch what is happening.

# Exercise 2: Determine the lattice constant of Si diamond crystal<a id="sec-2"></a>

## Part A: Generate Si crystal in the diamond configuration<a id="sec-2-1"></a>

Diamond lattice is composed of two FCC lattices shifted by 0.25.

The convetional cell is cubic, with lattice vectors:

```python
  a1 = [ 1.0, 0.0, 0.0 ]
  a2 = [ 0.0, 1.0, 0.0 ]
  a3 = [ 0.0, 0.0, 1.0 ]
```

and two unit cells of FCC, shifted by 0.25:

```python
  Natoms = 8

  basis[0] = [ 0.0, 0.0, 0.0 ]
  basis[1] = [ 0.5, 0.5, 0.0 ]
  basis[2] = [ 0.5, 0.0, 0.5 ]
  basis[3] = [ 0.0, 0.5, 0.5 ]
  basis[4] = [ 0.25, 0.25, 0.25 ]
  basis[5] = [ 0.75, 0.75, 0.25 ]
  basis[6] = [ 0.75, 0.25, 0.75 ]
  basis[7] = [ 0.25, 0.75, 0.75 ]
```

1.  Create the diamond structure in the conventional (cubic) cell, as given above. Verify the output structure by visualization.

## Part B: computing the lattice constant of Si<a id="sec-2-2"></a>

We can compute the lattice constant of a crystal, by finding the value of `a0` which minimizes the total energy of the cell. In this exercise, we do the task manually, meaning we explicitly select some values of `a0`, generate the structure with given `a0`, compute the corresponding energy, and plot the curve. The curve can then be fitted with a polynomial to find the exact value of `a0` at the minimum.

If you feel confident in scripting, skip directly to EXTRA-2.1.

1.  Generate the diamond structure with several different values of the lattice constant `a0`. The range of values should be around 5.2 to 5.8.
2.  Compute the value of total energy for each `a0`. Use LAMMPS software to compute the energy.
    
    Launch `lammps` with the command `lmp` as:
    
    ```bash
       lmp -in lammps.in
    ```
    
    this will produce some output on the screen. Search for the line that says: "TotEng". The value printed after it is the total energy of your system.
    
    Save this value into a text file, together with the current value of `a0`. The goal is to produce a file with contents like:
    
    ```txt
       5.20  	-33.795071
       5.25  	-34.151899
       5.30  	-34.414694
       5.35  	-34.588440
       5.40  	-34.677817
       5.45  	-34.687223
       5.50  	-34.620807
       5.55  	-34.482483
       5.60  	-34.275956
       5.65  	-34.004740
       5.70  	-33.672175
       5.75  	-33.281449
       5.80  	-32.835608
    ```
    
    Which can then be plotted by e.g. `gnuplot` to produce something like:
    
    ![img](./figs/a0_etot.png "Plot of total energy versus the lattice constant.")

To compute the total energy, use the script

compute the total energy for various values of the lattice parameter \`a0\`. Plot the enegy per atom vs value of \`a0\`, and find the minimum.

## EXTRA:<a id="sec-2-3"></a>

### EXTRA-2.1: Write a `bash` script<a id="sec-2-3-1"></a>

Create a `bash` script to compute the total energy corresponding to a range of values `a0`. The steps should generally follow the ones outlined in exercise description.

A possible script is written in the file `scripting/script.sh`. Can you make any sense of it?

### EXTRA-2.2: Find a polynomial fit to the data<a id="sec-2-3-2"></a>

Use your favourite programming language to write a functional fit to the data, and extract the value of the minimum `a0`.

### EXTRA-2.3: Compute with the primitive cell of diamond<a id="sec-2-3-3"></a>

Modify the script to create the primitive cell of diamond, and plot the result. What differences do you see (in terms of values of total energy, or `a0`)?

### EXTRA-2.4: Do we need a supercell?<a id="sec-2-3-4"></a>

Repeat the exercise with some supercell instead of the unit cell. Do you get any differences? Why?

# md; melt<a id="sec-3"></a>

## verlet/leapfrog, total energy, thormostat, damping etc.<a id="sec-3-1"></a>

## computing the force: empirical potentials<a id="sec-3-2"></a>

LJ

## LJ/melt in lammps examples<a id="sec-3-3"></a>

# create surface<a id="sec-4"></a>

create (by hand?) Si diamond in conventional cell (8atms per cell)

# minimize surface<a id="sec-5"></a>

surface as defect, minimize Si diamond with reaxff? maybe sw.

# difference between md and minimization<a id="sec-6"></a>

Taking energy out of the system. damped dynamics, quickmin/SD, CG, other &#x2026;

## try the melt example with adding damp at the end<a id="sec-6-1"></a>

## <a id="sec-6-2"></a>
