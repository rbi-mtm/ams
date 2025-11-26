- [Exercise 1: Generating a crystal structure](#sec-1)
  - [Part A: Generate an FCC structure in the cubic (conventional) cell:](#sec-1-1)
  - [Part B: Generate an FCC structure in the primitive cell:](#sec-1-2)
  - [Part C: Create and visualize a 4x4x4 supercell](#sec-1-3)
  - [EXTRA:](#sec-1-4)
- [Exercise 2: Determine the lattice constant of Si diamond crystal](#sec-2)
  - [Part A: Generate Si crystal in the diamond configuration](#sec-2-1)
  - [Part B: computing the lattice constant of Si](#sec-2-2)
  - [EXTRA:](#sec-2-3)
- [Exercise 3: velocity Verlet algorithm, Molecular Dynamics](#sec-3)
  - [Verlet-Stormer](#sec-3-1)
  - [velocity-Verlet](#sec-3-2)
  - [Computing the Force (classical potentials)](#sec-3-3)
  - [Molecular Dynamics](#sec-3-4)
  - [A word on statistical ensembles](#sec-3-5)
  - [Part A: Melting a Lennard-Jones solid](#sec-3-6)
  - [EXTRA:](#sec-3-7)
- [md; melt](#sec-4)
  - [verlet/leapfrog, total energy, thormostat, damping etc.](#sec-4-1)
  - [computing the force: empirical potentials](#sec-4-2)
  - [LJ/melt in lammps examples](#sec-4-3)
- [create surface](#sec-5)
- [minimize surface](#sec-6)
- [difference between md and minimization](#sec-7)
  - [try the melt example with adding damp at the end](#sec-7-1)
  - [](#sec-7-2)


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

1.  Modify the `generate_cell.py` script to generate the conventional cell of FCC structure, as given above. Call the script with `python`, and redirect the output to a file:
    
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

2.  Use `ovito` to visualize. The script will generate two structures in the same output file: one for the unit cell, and one for the supercell. Navigate the structures in `ovito`.

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
    
    Launch `lammps` with the command `lmp`, giving the input file `lammps.in` as:
    
    ```bash
       lmp -in lammps.in
    ```
    
    this will produce some output on the screen. Search for the line that says: "TotEng". The value printed after it is the total energy of your system.
    
    Save the current value of `a0`, together with the corresponding total energy value into a file. The goal is to produce a file with contents (`a0` in first column, and total energy in the second):
    
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

## EXTRA:<a id="sec-2-3"></a>

### EXTRA-2.1: Write a `bash` script<a id="sec-2-3-1"></a>

Create a `bash` script to compute the total energy corresponding to a range of values `a0`. The steps should generally follow the ones outlined in exercise description.

A possible script is written in the file `scripting/script.sh`. Can you make any sense of it?

### EXTRA-2.2: Find a polynomial fit to the data<a id="sec-2-3-2"></a>

Use your favourite programming language to write a functional fit to the data, and extract the value of the minimum `a0`.

### EXTRA-2.3: Compute with the primitive cell of diamond<a id="sec-2-3-3"></a>

Modify the script to work with the primitive cell of diamond, and plot the result. What differences do you see (in terms of values of total energy, or `a0`)?

### EXTRA-2.4: How about a supercell?<a id="sec-2-3-4"></a>

Repeat the exercise with some supercell instead of the unit cell. Do you get any differences? Why?

# Exercise 3: velocity Verlet algorithm, Molecular Dynamics<a id="sec-3"></a>

## Verlet-Stormer<a id="sec-3-1"></a>

In a computer, the first-order derivatives are computed as finite-difference: $$ \frac{\mathrm{d}}{\mathrm{d}t}x(t) = \frac{ x(t+\Delta t) - x(t) }{\Delta t} $$ with a chosen small value of $\Delta t$.

Let's label the positions as function of the timestep:

-   the current positions $x(t)=x_n$;
-   the positions in the next time-step $x(t+\Delta t)=x_{n+1}$;
-   the positions in the previous timestep $x(t-\Delta t)=x_{n-1}$.

Then, the second-order derivatives are computed as difference-of-differences: $$ \frac{\mathrm{d}^2}{\mathrm{d}t^2} x(t) = \frac { \frac{ x_{n+1} - x_n }{\Delta t} - \frac{ x_n - x_{n-1} }{\Delta t} }{\Delta t} = \frac { x_{n+1} - 2x_{n} + x_{n-1} } { \Delta t^2} $$

Then we can write the Newton's equation of motion: $$ F = ma = m \frac{\mathrm{d}^2}{\mathrm{d}t^2} x(t) = m \frac { x_{n+1} - 2x_{n} + x_{n-1} } { \Delta t^2} $$

Reshuffling the terms, to solve for positions in the next time-step: $$ x_{n+1} = 2x_n - x_{n-1} + \frac{F \Delta t^2}{m} $$

Therefore, to compute the future positions of atoms, we need to know the current positions, the positions of the previous time-step, and the force at the current time step.

The instantaneous force is a quantity which depends on the level of theory we want to include, but it can always be computed from the available data. At the simplest level, it's a function of the chemical type, and interatomic distances.

This is the basic Verlet-Stormer algotihrm for the update of atomic positions. Solving it for a number of timesteps $\Delta t$ makes a trajectory, which is just one of the possible solutions of the equation of motion.

## velocity-Verlet<a id="sec-3-2"></a>

The knowledge of positions at the previous time-step of the basic Verlet-Stormer algorithm poses a challenege when initiating the simulation, but can be overcome by the velocity-Verlet algorithm.

The velocity-Verlet algorithm uses particle velocities and forces to calculate the position update. The algorithm follows three steps:

-   calculate the future positions from current positions, current velocity, and current force: $$ x_{n+1} = x_n + v_n \Delta t + \frac{1}{2}F_n\Delta t^2 $$
-   compute the future force $F_{n+1}$ at positions $x_{n+1}$;
-   compute the future velocities: $$ v_{n+1} = v_{n} + \frac{1}{2}(F_n + F_{n+1})\Delta t $$

In order to start the velocity-Verlet algorithm, we need to provide the initial positions, and velocities.

## Computing the Force (classical potentials)<a id="sec-3-3"></a>

The instantaneous force on a configuration of particles $F_n$ can be computed at several different levels of theory. In this section we look at the [Lennard-Jones](https://en.wikipedia.org/wiki/Lennard-Jones_potential) (LJ) potential, which is possibly the simplest pair-potential. It gives the potential energy of two interacting objects, as the function of only the distance $r$ between them. As it does not contain any electronic effects, LJ is often referred to as a "classical" potential. Many other potentials of this type exist.

The LJ potential is defined: $$ V_{LJ}(r) = 4\epsilon \bigg[ \big(\frac{\sigma}{r}\big)^{12} - \big(\frac{\sigma}{r}\big)^{6} \bigg] $$ where $\epsilon$ is the potential depth, and $\sigma$ is the value of distance where $V_{LJ}=0$. It reaches a minimum value at $r=r_{m}=2^{1/6}\sigma$.

Sometimes it is convenient to write it as: $$ V_{LJ}(r) = \frac{A}{r^{12}} - \frac{B}{R^{6}} $$

![img](./figs/LJpot.png "The LJ potential.")

The force can be in general computed as the negative gradient of the potential, but $V_{LJ}$ is spherically symmetric, so the force direction is given by the vector $\hat{r}$, and its magnitude by the simple derivate, with an analytical expression: $$ \vec{F}_{LJ}( \vec{r} ) = -\hat{r}\frac{\mathrm{d}V}{\mathrm{d}r} = \hat{r} 48\epsilon \bigg[ \frac{\sigma^{12}}{r^{13}} - 0.5 \frac{\sigma^{6}}{r^{7}} \bigg] $$ where $\vec{r}$ is the vector connecting two particles, $r=|\vec{r}|$ is its Cartesian norm, and $\hat{r} = \vec{r}/r$ its direction.

## Molecular Dynamics<a id="sec-3-4"></a>

Molecular Dynamics (MD) is a simulation technique which generates a trajectory of configurations, which is a possible solution of the equations of motion. By changing the initial conditions, MD generates a different trajectory which is also a valid solution.

Very often, the concrete method of computing the atomic update step is similar (if not identical) to the one sketched above (the velocity-Verlet). The computation of force can however be much more involved, by including effects at several levels of theory (classical, semi-classical, QM/MM, ab-initio, etc.).

Once we have an equilibrated, and long-enough trajectory, we can compute some thermodynamic properties from the trajectory. Namely, we can compute those that can be expressed as functions of positions and momenta of the particles. For example temperature, expressed as average kinetic energy per particle.

## A word on statistical ensembles<a id="sec-3-5"></a>

Depending on which statistical ensemble we wish to simulate (canonical, micro-canonical, or grand-canonical), different thermodynamic properties are kept constant in the simulation, while others can vary.

If we wish to simulate a system at a certain constant temperature, the particle properties (velocity) have to be scaled to meet the temperature criteria. There exist different schemes for such rescaling, commonly called the "thermostats".

Similarly, if we wish to simulate a system at a certain constant pressure, the schemes are called "barostats".

The different ensembles are often referred to by which properties are constant:

-   micro-canonical NVE: the number of particles/moles $N$, box volume $V$, and the energy $E$. This corresponds to an isolated system, which does not exchange heat or matter with its surrounding. The total energy in NVE is coserved (i.e. the sum of kinetic and potential $E=K+V$).

-   canonical NVT: the number of particles/moles $N$, box volume $V$, and the temperature $T$. The system is allowed to exchange energy (heat) with its surrounding, such that $T$ remains constant. Velocities rescaled with a thermostat, to keep constant the kinetic energy (=temperature).

-   isothermal-isobaric NPT: the number of particles/moles $N$, the pressure $P$, and temperature $T$. The system is allowed to exchange energy (heat), and change the volume (box size) to keep the pressure constant. The pressure is manipulated by rescaling the box, via a barostat.

-   grand-canonical $\mu$VT: the chemical potential $\mu$, volume $V$, and temperature $T$. The system is assumed open, and can exchange heat (energy), and particles (matter, chemical elements) with its surrounding. The exchange of particles is very difficult to properly implement in MD, some Monte-Carlo based methods exist though.

## Part A: Melting a Lennard-Jones solid<a id="sec-3-6"></a>

1.  launch lmp

2.  visulaize output, plot the radial distribution function.

3.  plot the kinetic, potential, and total energy as function of time

## EXTRA:<a id="sec-3-7"></a>

### EXTRA-3.1: significance of `epsilon` and `sigma`<a id="sec-3-7-1"></a>

Imagine you try to model a real molecule with two atoms with the LJ potential. Which experimantal (or computed) values of the molecule would you need, in order to determine `epsilon`? And which to determine `sigma`?

### EXTRA-3.2: implement your own module to compute LJ energy and force:<a id="sec-3-7-2"></a>

Using your favourite programming language, implement functions/routines to compute the total energy, and the force vectors for a given configuration of particles. Use lammps calculations as reference values.

NOTE: keep in mind the distances have to be computed in Periodic Boundary Conditions. NOTE 2: the range of distances where LJ gives a nonzero potential can be quite large, meaning we need to explicitly include a quite large cell in the calculation. In the generic PBC implementation, each image is only made to interact with its nearest-neighbor image, but not the second, third, etc. This is often mitigated by introducing a cutoff range to the LJ interactions, beyond which the interactions are added as analytical expressions based on the density.

### EXTRA-3.3: implement your own Verlet algorithm:<a id="sec-3-7-3"></a>

Using your favourite programming language, implement the Velocity-Verlet algorithm. To compute the energy and forces, use the routines you previously implemented.

Attention to the units, and the $\Delta t$ parameter.

### EXTRA-3.4: compute which Bravais lattice has the lowest LJ energy<a id="sec-3-7-4"></a>

Use the script `generate_cell.py` from the first exercise to generate the Face-Centered Cubic (FCC), Body-Centered Cubic (BCC), and Hexagonal Close-Packed (HCP) structures. Then choose parameters for the `epsilon` and `sigma`, and find the lattice constant `a0` for each crystal lattice. Generate a large supercell with the final parameters of each lattice, and compute the total energy per atom. Which crystal structure has the lowest energy?

The structure with lowest energy is the most stable, and coincidentally appears most often in the nature. It is quite remarkable that a simple LJ potential is already sufficient to reproduce that behaviour of nature. Remember the LJ only has a radial component.

### EXTRA-3.5: LJ substance<a id="sec-3-7-5"></a>

The substance/matter of an LJ simulation is sometimes called the "Lennard-Jonesium". Historically, LJ potential was used quite successfully to model the liquid Argon. Find the values of `epsilon` and `sigma` online to model it, and try to simulate some representative points of its phase diagram (solid, liquid, gas).

# md; melt<a id="sec-4"></a>

## verlet/leapfrog, total energy, thormostat, damping etc.<a id="sec-4-1"></a>

## computing the force: empirical potentials<a id="sec-4-2"></a>

LJ

## LJ/melt in lammps examples<a id="sec-4-3"></a>

# create surface<a id="sec-5"></a>

create (by hand?) Si diamond in conventional cell (8atms per cell)

# minimize surface<a id="sec-6"></a>

surface as defect, minimize Si diamond with reaxff? maybe sw.

# difference between md and minimization<a id="sec-7"></a>

Taking energy out of the system. damped dynamics, quickmin/SD, CG, other &#x2026;

## try the melt example with adding damp at the end<a id="sec-7-1"></a>

## <a id="sec-7-2"></a>
