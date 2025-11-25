# AMS25 / Day 2 / Hands-on Session

Welcome to the second day of the hands-on sessions at **AMS25**!

## Introduction

Last week we learned about general concepts useful for all kinds of atomistic simulations.
Today, we'll be specifically performing **density functional theory (DFT)** calculations.

Compared to last week's simulations that employed classical force fields,
DFT calculations are quantum mechanical.
As a consequence, they are more **accurate**, **informative** and **computationally expensive**.

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

## Exercise 1: Diamond Silicon

In this exercise, we will be exploring some of the ground-state properties of silicon.

Since we are working within the **density** functional theory framework,
we first need to obtain the **ground state electronic density**.
Let's approach this problem step by step.

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
cd ex1
less si_scf.in
```

This input is about as simple as it gets.
Many more options are listed in the [QE documentation](https://www.quantum-espresso.org/Doc/INPUT_PW.html).

> [!NOTE]
> Each SCF input file must contain the `&CONTROL`, `&SYSTEM` and `&ELECTRONS` namelists,
> even though `&ELECTRONS` can be left empty (default values will be used).

Let's run the SCF calculation by using `pw.x`, the basic program from the QE suite:

```bash
pw.x -i si_scf.in | tee si_scf.out
```

> [!TIP]
> We'll be using `tee` to simultaneously write the output to a file
> and to the screen to follow calculations in real time.

We can open the output file:

```bash
less si_scf.out
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

**Q1:** We specified a grid of $$6 \times 6 \times 6 = 216$$ $k$-points. However, how many $k$-points are actually used in the calculation? Can you think of an explanation why?

**Q2:** How many iterations did it take to converge the ground state density?
What is the total energy trend over the iterations?<br>
(Hint: use `grep "total energy" si_scf.out | head -n -1`)

**Q3:** What are the forces on the two silicon atoms?

> [!TIP]
> It is interesting to check the different contributions to the total energy:
> ```bash
> grep "!" si_scf.out -A7
> ```
> - One-electron contribution: kinetic energy of KS electrons
> - Hartree contribution: classical electrostatic electronic interaction
> - XC contribution: exchange-correlation energy (XC model dependent)
> - Ewald contribution: (pseudo)-ionic interaction computed using the [Ewald method](https://en.wikipedia.org/wiki/Ewald_summation)
>
> The absolute value of energy is pseudopotential dependent; it has no physical meaning.

Let's also take a look at energies for a fixed $k$-point. E.g., for $k=(0,0,0)$, the so called $\Gamma$ point, the KS energies $\epsilon_{ik}$ are:

```bash
k = 0.0000 0.0000 0.0000 (   339 PWs)   bands (ev):

-5.8718   6.0677   6.0677   6.0677   8.6242   8.6242   8.6242   9.3285
```

There are 8 bands because we specified `nbnds=8` in the input file.
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

We have 16 KS wavefunction files (one for each $k$-point)
and the `charge-density.dat` file, all in binary format.

Let's plot $n(\vec{r})$.
First, we have to convert the `charge-density.dat` file to a more suitable format.
We can do that using the QE postprocessing program, [`pp.x`](https://www.quantum-espresso.org/Doc/INPUT_PP.html):

```bash
pp.x -i si_pp_charge.in | tee si_pp_charge.out
```

Visualize the density using [XCRYSDEN](http://www.xcrysden.org/):

```bash
xcrysden --xsf Si.charge.xsf
```

> [!TIP]
> For a clearer image, first click `Modify` &rarr; `Number of Units Drawn`
> and create a supercell. The density itself is visualized using
> `Tools` &rarr; `Data Grid`.

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
diff -y si_scf.in si_nscf.in | less
```

and run the NSCF calculation:

```bash
pw.x -i si_nscf.in | tee si_nscf.out
```

**Q:** How many wave-functions were obtained by the NSCF calculation?
