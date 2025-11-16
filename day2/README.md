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
around a core of 14 protons, we model the system as
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
