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
we first need to obtain the ground state electronic density.
Let's approach this problem step by step.

### Exercise 1.1: Finding the Atomic Structure

The silicon atomic structure was solved using X-ray diffraction back in the 1920s, i.e. it's very well-known.
[Silicon](https://www.chemicalbook.com/article/silicon-diamond-cubic-crystal-structure.htm) crystallizes in a [diamond cubic structure](https://en.wikipedia.org/wiki/Diamond_cubic),
which can be seen as two identical FCC lattices shifted by $1/4$ relative to each other.
