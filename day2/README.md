# AMS25 / Day 2 / Hands-on Session

Welcome to the second day of the hands-on sessions at **AMS25**!

## Introduction

Last week we learned about general concepts useful for all kinds of atomistic simulations.
Today, we'll be specifically performing **density functional theory (DFT)** calculations.

Compared to last week's simulations that employed classical force fields,
DFT calculations are quantum mechanical.
As a consequence, they are more **accurate**, **informative** and **computationally expensive**.

> [!NOTE]
> On first reading, you might want to skip
> the additional exercises (labeled with **Task**).

## About Quantum ESPRESSO

![QE-logo](https://gitlab.com/QEF/q-e/-/raw/develop/logo.jpg)

We'll be performing all calculations using [**Quantum ESPRESSO (QE)**](https://www.quantum-espresso.org/). QE is [free](https://en.wikipedia.org/wiki/GNU_General_Public_License) and [open-source](https://gitlab.com/QEF/q-e) software.

> [!NOTE]
> QE is only [one of many](https://en.wikipedia.org/wiki/Category:Density_functional_theory_software) software packages implementing DFT.
>
> The differences between the packages are in the intended **use-cases** (some are more optimal for molecules, others for materials),
> **basis sets** (QE uses plane-waves, some use Gaussians or more complicated functions),
> some explicitly take into account **all electrons**, while others (like QE) model the
> electrons near the core using **pseudopotentials**, ...

Some facts about QE:
- created in 2002
- several hundred contributors over the years
- approximately 520k lines of Fortran source code

[What can QE do?](https://www.quantum-espresso.org/what-can-qe-do/)
- **Ground-state calculations**
- **Structural optimizations / molecular dynamics**
- Response properties
- Quantum transport
- ...

## Exercise 1: Basic QE workflow / Silicon

In this exercise, we will be showcasing a basic QE workflow
using the example of silicon in its
naturally occuring (diamond cubic) form.

The workflow consists of:
1. Choose an approximation (exchange-correlation functional)
(we will use [LDA](https://en.wikipedia.org/wiki/Local-density_approximation))
2. Input the crystal structure into QE
3. Get the ground-state electronic density

From the density, we easily get the energy and
[ionic forces](https://en.wikipedia.org/wiki/Hellmann%E2%80%93Feynman_theorem).

However, it turns out that we can learn more about the
details of the ground-state, as well as get a
first approximation to the excited states, by looking
closer into the Kohn-Sham orbitals, which we
obtained on the way, while getting the density.

Therefore, step 4 of the simplified workflow is to investigate
the Kohn-Sham states.

***
### Exercise 1.1: Finding the Crystal Structure

The silicon atomic structure was solved using X-ray diffraction back in the 1920s, i.e. it's very well-known.
Silicon crystallizes in a [diamond crystal structure](https://lampz.tugraz.at/~hadley/ss1/crystalstructure/structures/diamond/diamond.php)
which can be viewed as two identical FCC lattices shifted by $`1/4`$ relative to each other.

> [!NOTE]
> It's [possible](https://www.quantum-espresso.org/faq/faq-input-data/#:~:text=3.1%20DOES%20QE%20USE%20PRIMITIVE%20OR%20CONVENTIONAL%20UNIT%20CELL?,and%20the%20atomic%20positions%20correctly.) to use either the _conventional_ or _primitive_
> unit cells in QE. We will be using the **primitive** cell.

To summarize, these are the parameters we'll be needing:

$$
\begin{align}
&a_\text{lat} = 5.431 \text{ &#197;} = 10.263\ \text{Bohr}, \\
&r_{\text{Si}_1} = (0, 0, 0), \\
&r_{\text{Si}_2} = (0.25, 0.25, 0.25),
\end{align}
$$

where the atomic positions are written in [crystal (fractional) coordinates](https://en.wikipedia.org/wiki/Fractional_coordinates).

Generally, crystal structures can be found at:
- [Cambridge Structural Database (CSD)](https://www.ccdc.cam.ac.uk/)
- [Crystallography Open Database](https://www.crystallography.net/cod/)
- [Materials Project](https://next-gen.materialsproject.org/)
- Publications

***
### Exercise 1.2: Electron Configuration and Pseudopotentials

Before delving in DFT calculations, it is worthwhile to spend time on elementary considerations.
The [atomic number of silicon is 14](https://periodic-table.rsc.org/),
meaning that a neutral atom of silicon has 14 protons and 14 electrons.
Since electrons are fermions, they must obey the Pauli exclusion principle,
i.e., a state with given quantum numbers can be occupied by, at most, one electron.

This gives rise to the following [electron configuration](https://en.wikipedia.org/wiki/Electron_configuration) of a silicon atom:

$$
\begin{align}
1\text{s}^2 2\text{s}^2 2\text{p}^6 3\text{s}^2 3\text{p}^2
\end{align}
$$

or in shorthand notation:

$$
\begin{align}
\[\text{Ne}\]\ 3\text{s}^2 3\text{p}^2.
\end{align}
$$

The 10 electrons (written as [Ne]) are usually defined as the **core** electrons,
while the 4 remaining ones ($`3\text{s}^2 3\text{p}^2`$) are called **valence** electrons.

The _core_ electrons are called as such because they are strongly bound
to the atomic core. In a crystal, core electrons do not participate
in the formation of chemical bonds; rather, they remain bound and **localized**
near their original atom.

On the other hand, the bands near the [Fermi level](https://en.wikipedia.org/wiki/Fermi_level)
(i.e. the conduction and valence bands)
are the ones that primarily determine the interesting reponse properties of a material. These bands are formed by the _valence_ electrons.

Due to the above, QE uses a concept called [_pseudopotentials_](https://en.wikipedia.org/wiki/Pseudopotential),
which is a way to remove the core electrons from the calculation.
In the case of silicon, instead of having 14 electrons
around a core of 14 protons, we can model the system as
4 electrons around an effective core consisting of **14 protons + 10 electrons**.

This effective description of the core manifests itself in
the Kohn-Sham equations by the electron-ion interaction being
replaced by a _pseudopotential_ (i.e. electron-effective ion interaction).

> [!NOTE]
> There is a (fundamentally mathematical) reason why using
> pseudopotentials in QE is, in fact, not a _choice_, but a
> _necessity_. 
>
> Namely, QE employs a [plane-wave basis set](https://en.wikipedia.org/wiki/Basis_set_(chemistry)). Plane-waves are completely delocalized
> in real space. As a consequence, it takes a lot of coefficients
> to represent a spatially localized function in the plane-wave basis.
>
> This is a practical issue since the wave functions
> of core electrons are localized.
> Furthermore, even the valence electron wave functions
> have nodes near the core as they have to be orthogonal
> to the core electron wave functions.
>
> Some DFT codes employ localized basis sets and, as such,
> are able to treat all electrons explicitly.
> One such example is [WIEN2k](https://en.wikipedia.org/wiki/WIEN2k).

Take a look at the pseudopotential file
that we'll be using for silicon calculations:

```bash
less pseudo/Si.pz-vbc.UPF
```

This file is written in the [Unified Pseudopotential Format](https://pseudopotentials.quantum-espresso.org/home/unified-pseudopotential-format).
We can see that `Z_valence = 4.0`.

> [!TIP]
> For simple scrolling with `less`, use the arrows (&uarr;&darr;).
> Exit `less` by typing `q`.

***
### Exercise 1.3: Self-Consistent Field Calculation

After preliminary considerations, let's calculate
the ground state density $$n(\vec{r})$$.
The procedure by which $$n$$ is obtained is called a
**self-consistent field (SCF)** calculation.

Spend some time checking the SCF input file:

```bash
cd ~/ams/day2/exercise_1
less 01_si_scf.in
```

This input is about as simple as it gets.
Many more options are listed in the [QE documentation](https://www.quantum-espresso.org/Doc/INPUT_PW.html).

> [!NOTE]
> Each SCF input file must contain the `&CONTROL`, `&SYSTEM` and `&ELECTRONS` namelists,
> even though `&ELECTRONS` can be left empty (default values will be used).

Let's run the SCF calculation by using `pw.x`, the basic program from the QE suite:

```bash
pw.x -i 01_si_scf.in | tee 01_si_scf.out
```

> [!TIP]
> We'll be using `tee` to simultaneously write the output to a file
> and to the screen to follow calculations in real time.

> [!NOTE]
>
> **FAQ:** How to select `ecutwfc`, the `K_POINTS` grid and such "free" parameters?
>
> These parameters must be selected large enough
> so that the results of the calculation (e.g. energy)
> **do not significantly change**
> if such parameters are increased further.
>
> In principle, one should always perform a _convergence test_ to be sure
> that the selected parameters are large enough.
>
> **Task:** Calculate the dependence of energy for `ecutwfc` values
> of 10, 20, 30, 40, 50, 60 and 70 Ry.

We can open the output file:

```bash
less 01_si_scf.out
```

To understand the output, let's remind ourselves what's an SCF calculation.

We want to obtain the ground state electron density $$n$$.
We can't do that directly for the real many-electron system
(because its Schrödinger equation is computationally intractable).
However, we can construct an auxillary system of noninteracting fermions,
called Kohn-Sham (KS) electrons, that is *solvable* and whose ground state
density is identical (by construction) to the density of the real system!

<p align="center">
  <img src="/day2/figs/scf_diagram.png" width="550">
</p>

> [!NOTE]
> $i$ and $k$ are the quantum numbers of the KS wavefunctions $\phi_{ik}$.
> - $i$ is the KS band; the occupied states are given by $$i \in [1, N_\text{el}/2]$$ (for a **spin-unpolarized** calculation)
>   - In our case, $N_\text{el} = 4 \times 2 = 8$, so we have **4 occupied (degenerate) bands**
> - $k$ is the wave vector; in QE we define the $k$-points for which the calculation will be performed in the `K_POINTS` card

**Question:** We specified a grid of $$12 \times 12 \times 12 = 1728$$ $k$-points. However, how many $k$-points are actually used in the calculation? Can you think of an explanation why?

**Question:** How many iterations did it take to converge the ground state density?
What is the total energy trend over the iterations?<br>
(Hint: use `grep "total energy" si_scf.out | head -n -1`)

**Question:** What are the forces on the two silicon atoms?

> [!TIP]
> It is interesting to check the different contributions to the total energy:
> ```bash
> grep "!" 01_si_scf.out -A7
> ```
> - One-electron contribution: kinetic energy of KS electrons
> - Hartree contribution: classical electrostatic electronic interaction
> - XC contribution: exchange-correlation energy (XC model dependent)
> - Ewald contribution: (pseudo)-ionic interaction computed using the [Ewald method](https://en.wikipedia.org/wiki/Ewald_summation)
>
> The absolute value of energy is pseudopotential dependent; it has no physical meaning.

Let's also take a look at energies for a fixed $k$-point. E.g., for $k=(0,0,0)$, the so called $\Gamma$ point, the KS energies $\epsilon_{ik}$ are:

```bash
k = 0.0000 0.0000 0.0000 (  2109 PWs)   bands (ev):

-5.8889   6.0358   6.0358   6.0358   8.5968   8.5968   8.5968   9.2920
```

There are 8 bands because we specified `nbnds = 8` in the input file.
We also have:

```bash
highest occupied, lowest unoccupied level (ev):     6.0677    6.7311
```

As expected, 4 bands are occupied and 4 are unoccupied.
The occupied band energies are obtained from the SCF cycle,
which fixes the KS Hamiltonian $H_\text{KS}$, allowing us
to also obtain the unoccupied band energies.

***
### Exercise 1.4: Ground State Density

Take a look at the files QE saved after completion of an SCF calculation:

```bash
ls out/silicon.save/
```

We have 72 KS wavefunction files (one for each $k$-point)
and the `charge-density.dat` file, all in binary format.

Let's plot $n(\vec{r})$.
First, we have to convert the `charge-density.dat` file to a more suitable format.
We can do that using the QE postprocessing program, [`pp.x`](https://www.quantum-espresso.org/Doc/INPUT_PP.html):

```bash
pp.x -i 02_si_pp_charge.in | tee 02_si_pp_charge.out
```

Visualize the density using [XCRYSDEN](http://www.xcrysden.org/):

```bash
xcrysden --xsf Si.charge.xsf
```

> [!WARNING]
> If you are on an ARM machine and XCRYSDEN does not work,
> here's a possible fix:
>
> ```bash
> sudo apt-get -y install xcrysden && \
> mkdir -p ~/.xcrysden && \
> echo "set toglOpt(accum) false" > ~/.xcrysden/custom-definitions
> ```

> [!TIP]
> For a clearer image, first click `Modify` &rarr; `Number of Units Drawn`
> and create a supercell. The density itself is visualized using
> `Tools` &rarr; `Data Grid`.

The result should look something like the image below.

<p align="center">
  <img src="/day2/figs/Si_gs_density.png" width="800">
</p>

The charge density is concentrated _between_ atoms,
indicative of covalent bonding.

***
### Exercise 1.5: Non-Self-Consistent Field Calculation

Even though the KS system is basically a fictitious system constructed
as a crutch to obtain $n(\vec{r})$, it turns out that the KS particles
are a good zero-order approximation to real **single-particle excitations (quasiparticles)**
in weakly-correlated systems.

> [!NOTE]
> A popular approach to quasiparticle calculations 
> is the [GW approximation](https://www.cond-mat.de/events/correl11/manuscript/Held.pdf),
> which goes beyond zero-order by approximating the so-called _self-energy_ contribution to
> quasiparticle energy.
>
> A common starting point for a GW calculations is a DFT SCF calculation.

Therefore, a closer look at the KS eigenvalues $\epsilon_{ik}$ and eigenfunctions $\phi_{ik}$
can provide us with useful information about the real system.

First of all, now that we know $n$ from the SCF calculation, we can better sample the
$k$-points by performing a _non-self-consistent field_ (NSCF) calculation (see diagram below).

<p align="center">
  <img src="/day2/figs/nscf_diagram.png" width="350">
</p>

Let's compare the SCF and NSCF input files:

```bash
diff -y 01_si_scf.in 03_si_nscf.in | less
```

and run the NSCF calculation:

```bash
pw.x -i 03_si_nscf.in | tee 03_si_nscf.out
```

**Question:** How many wave-function files were obtained by the NSCF calculation?

***
### Exercise 1.6: Density of States

Now that we have a denser sampling of the KS states in the reciprocal space,
let's look at the [density of states (DOS)](https://en.wikipedia.org/wiki/Density_of_states).

For this purpose, we can use the [`dos.x`](https://www.quantum-espresso.org/Doc/INPUT_DOS.html) postprocessing program:

```bash
dos.x -i 04_si_dos.in | tee 04_si_dos.out
```

We got the output file `Si.dos.dat`. Let's plot it:

```bash
python plot_dos.py Si.dos.dat
```
<p align="center">
  <img src="/day2/figs/Si_dos.png" width="650">
</p>

We can see that there are no KS states around Fermi level.
From that, we can conclude that the single-particle excitation spectrum is also likely gapped,
or in simpler terms: the system is an **insulator**
(in fact, since the gap is small, the system is classified as a **semiconductor**).

**Task:** Plot the integral of the DOS versus energy.

**Question:** What does the DOS integrate to? What's the integral of the DOS up to Fermi energy?

> [!NOTE]
> The band gap is about 0.55 eV, which is **much** smaller than the experimentally
> measured value of about 1.12 eV. Besides the fact that we are using the LDA functional,
> there is a **fundamental** reason why practical DFT calculations underestimate bandgaps.
>
> Namely, it is known that the exact (but unknown) exchange-correlation functional
> has [derivative discontinuities at integer particle numbers](https://physics.stackexchange.com/questions/176419/why-does-density-functional-theory-dft-underestimate-bandgaps),
> while the functionals employed in practical calculations are continuous.
>
> This discontinuity contributes to the fundamental gap,
> leading to DFT band gaps being generally underestimated.

***
### Exercise 1.7: Band Structure & Projected Density of States

We continue with investigating the structure of the KS states.
In this Exercise, we'll plot the band structure, i.e. the dispersion of Kohn-Sham energies
$\epsilon_{i}(\vec{k})$.

These energies are also obtained using an NSCF calculation.
Recall that in Exercise 1.4 we used a **uniform grid** of $k$-points for an NSCF calculation.
This provided us with a uniform sampling of KS states across the [Brillouin zone](https://en.wikipedia.org/wiki/Brillouin_zone),
convenient for calculating the DOS.

On the other hand, to obtain an interesting band structure, we need to select a **path** of
$k$-points through the Brillouin zone for which $\epsilon_{i}(\vec{k})$ will be evaluated.
This can conveniently be done in XCRYSDEN.
We'll be using the following path that connects [high-symmetry points](https://en.wikipedia.org/wiki/Brillouin_zone#Critical_points): $L-\Gamma-X-K-\Gamma$

<p align="center">
  <img src="/day2/figs/Si_band_path.png" width="650">
</p>

The input file is given in `05_si_bands.in`. Let's run the calculation:

```bash
rm out/silicon.save/wfc*
pw.x -i 05_si_bands.in | tee 05_si_bands.out
```

> [!IMPORTANT]
> We first removed the KS wave-functions obtained by the previous NSCF calculation
> so they don't get mixed up with the ones obtained in the `bands.x` calculation.
>
> In a realistic situation, we would first back up the old wave-functions.

**Task:** Reproduce the k-path given in `05_si_bands.in` using XCRYSDEN.

> [!NOTE]
> **FAQ:** How to select a path through the Brillouin zone?
>
> It depends on what we're interested in.
> High-symmetry points are where the bands usually have extremal
> (minimal or maximal) behavior, so we will include them if we want to, e.g., see
> if the system has a [direct or indirect band gap](https://en.wikipedia.org/wiki/Direct_and_indirect_band_gaps).
> The rest of the band structure can usually be deduced from the path
> between the high-symmetry points.
>
> We can also infer other useful quantities, such as
> the [effective mass](https://en.wikipedia.org/wiki/Effective_mass_(solid-state_physics)) of single particle excitations.
> Bands are also a starting point for fitting tight-binding
> models using [maximally-localized Wannier functions](https://wannier.org).
>
> There are tools which generate paths automatically ([SeeK path](https://seekpath.materialscloud.io/))
> but it's best to always check the literature to determine the path.

**Question:** How many wave-function files were obtained by this NSCF (bands) calculation?

We can write the bands in a convenient format using the `bands.x` postprocessing program:

```bash
bands.x -i 06_si_pp_bands.in | tee 06_si_pp_bands.out
```

We'll plot the bands using the `Si.bands.dat.gnu` file.
However, to better understand the nature of the bands, we'll first
project the KS orbitals onto atomic orbitals.
This is done using [`projwfc.x`](https://www.quantum-espresso.org/Doc/INPUT_PROJWFC.html):

```bash
mkdir -p pdos
projwfc.x -i 07_si_projwfc.in | tee 07_si_projwfc.out
```

**Question:** What are the atomic orbitals
that we are projecting the KS states onto?
Hint: Read the ``07_si_projwfc.out`` file and check the
``pdos`` directory.

Using the small helper-program `sumpdos.x`,
let's sum the contributions of the two atoms for
the $s$ and $p$ orbitals, respectively:

```bash
sumpdos.x pdos/*\(s\) > Si.pdos_s.dat
sumpdos.x pdos/*\(p\) > Si.pdos_p.dat
```

We can finally plot the bands and PDOS together:

```bash
python plot_bands_pdos.py
```
<p align="center">
  <img src="/day2/figs/Si_bands_pdos.png" width="1000">
</p>

- The band gap is indirect.
- We can see that the valence band has a mostly $p$-orbital character.
Again, this is indicative of directional (covalent) bonding.

**Task:** Calculate the effective mass of the
[heavy and light holes](https://physics.stackexchange.com/questions/229328/why-do-we-have-heavy-and-light-hole-bands-in-semiconductors)
in the $\Gamma-X$ and $\Gamma-L$ directions by fitting parabolas to the bands.
Compare them with the results on the bottom of
[this page](https://lampz.tugraz.at/~hadley/ss1/semiconductors/silicon_bandstructure.php).

**Task:** Generate a k-path automatically using
[SeeK-path](https://seekpath.materialscloud.io/)
and plot the resulting band structure.
Take care that the generated path is
[discontinous](https://physics.stackexchange.com/questions/545748/what-is-the-meaning-of-vertical-bars-in-paths-of-high-symmetry-points).

***
### Exercise 1.8: Additional Tasks

- Construct a $2\times 2 \times 2$ supercell of silicon
and repeat the calculations above to obtain the band structure.
Since, in the supercell, there are 8 times more electrons,
you get 8 times more occupied bands.
However, no new physics was introduced;
only the periodicity of the system was redefined.
What is the connection between the band
structures of the primitive cell and the supercell?

- Repeat all of the above calculations for [diamond (carbon)](https://en.wikipedia.org/wiki/Diamond). What is the value of the band gap? The necessary pseudopotential is given in the [pseudo](/day2/pseudo) directory.

***
## Exercise 2: DFT in Metallic Systems / Silver

In this Exercise, we'll perform similar elementary calculations
for a metallic system: silver (Ag).

### Exercise 2.1: Fermi Surface

Let's check the SCF input file:

```bash
cd ~/ams/day2/exercise_2
less 01_ag_scf.in
```

There are two main differences compared to silicon calculations:

- We're using the [PBE functional] https://dft.uci.edu/pubs/RCFB08.pdf
- Since Ag is a metal, we have to use [smearing techniques](https://vasp.at/wiki/Smearing_technique)
  to stabilize numeric convergence

Let's run the SCF calculation using 2 CPUs via [MPI](https://en.wikipedia.org/wiki/Message_Passing_Interface):

```bash
mpirun -n 2 pw.x -i 01_ag_scf.in | tee 01_ag_scf.out
```

It's interesting to plot the Fermi surface:

```bash
fs.x -i 02_ag_fs.in | tee 02_ag_fs.out
xcrysden --bxsf silver_fs.bxsf
```
<p align="center">
  <img src="/day2/figs/Ag_fermi_surface.png" width="800">
</p>

The Fermi surface is _almost_ a sphere - we can see
[necks](http://www.jetp.ras.ru/cgi-bin/dn/e_015_01_0049.pdf) around the $L$ points
of the Brillouine zone.

This can be understood in terms of the
[nearly free electron model](https://solidstate.quantumtinkerer.tudelft.nl/test_builds/N_A-equal-0/11_nearly_free_electron_model/)
for systems in which the
[equivalent free electron Fermi sphere extends beyond the Brillouin zone boundary](https://kbose.weebly.com/uploads/1/0/4/9/10492046/fermi_surfaces_from_kittel.pdf).
