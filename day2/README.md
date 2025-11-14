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
> electrons near the core using [**_pseudopotentials_**](https://en.wikipedia.org/wiki/Pseudopotential), ...

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
