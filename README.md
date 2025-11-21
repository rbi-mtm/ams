# Atomistic Modeling School 2025

Welcome to the **Atomistic Modeling School 2025 ([AMS25](https://rbi-mtm.github.io/teaching/ams25/))!**

The school is held at the [**Institute of Physics (Zagreb, Croatia)**](https://maps.app.goo.gl/SigjffQmq9kXdW2n9) across three days:

- December 3 2025, 09:00 - 17:00
- December 10 2025, 09:00 - 17:00
- December 17 2025, 09:00 - 17:00

This repository contains materials related to the hands-on sessions held at **AMS25**.

## Table of Contents
 1. [Introduction](#1-introduction)
 2. [Schedule](#2-schedule)
 3. [Setting Up](#3-setting-up)
    
    3.1. [Installing the VM](#31-installing-the-vm)

    3.2. [About the VM](#32-about-the-vm)

    3.3. [Local Installation](#33-local-installation)

 4. [Running the Exercises](#4-running-the-exercises)
 5. [Common Issues](#5-common-issues)

    5.1. [Issues Related to the Virtual Machine](#51-issues-related-to-the-vm)

    5.2. [Other Issues](#52-other-issues)

 6. [Software Installation](#6-software-installation)
    
## 1. Introduction

During the hands-on sessions at **AMS25**, you will be running all exercises using a [**virtual machine**](https://en.wikipedia.org/wiki/Virtual_machine) (VM).
The advantages of using a VM for a school like this are:

- a VM is isolated from your host computer; i.e., it is a sandbox in which you can experiment with software without any risk of messing up your main installation;

- we have preinstalled all the software necessary to run the exercises, meaning that we can focus on learning concepts, instead of wasting time on technicalities.

This document should be used in the following way: first, read the chapter on [Setting Up](#3-setting-up) and install the VM. When done, learn a bit [about the VM](#32-about-the-vm). Finally, check how to [run the exercises](#4-running-the-exercises).

***
> [!TIP]
> Some participants might not be able to
> [install the VM](#31-installing-the-vm)
> or run the tutorials.
> Therefore, we encourage you to work in groups.
> Be collegial and share both knowledge and technical resources!

## 2. Schedule

| Start | End | 3.12. | 10.12. | 17.12. |
|:---:|:---:|---|---|---|
| 8:30 | 9:00 | Registration | - | - |
| 9:00 | 9:55 | **Ivor Lončarić**<br>Opening / General Introduction to Atomistic Simulations | **Dino Novko**<br>Phonons, Electrons and Electron-Phonon Coupling | **Dino Novko**<br>Comparing _Ab Initio_ Simulations With Experiments |
| 9:55 | 10:25 | Coffee Break | Coffee Break | Coffee Break |
| 10:25 | 11:20 | **Luka Benić**<br>Fundamentals of Density Functional Theory | **Ryan Requist**<br>_Ab Initio_ Molecular Dynamics | **Ivor Lončarić**<br>Machine learning in Materials Science |
| 11:20 | 11:25 | Break | Break | Break |
| 11:25 | 12:20 | **Miha Gunde**<br>Introduction to Crystallography & Atomistic Simulations | **Bernhard Kretz**<br>Theory and Modelling of Electronic Transport | **Bernhard Kretz**<br>Hands-on Session #1 |
| 12:20 | 13:50 | Lunch | Lunch | Lunch |
| 13:50 | 14:40 | **Miha Gunde**<br>Hands-on Session #1 | **Juraj Ovčar**<br>Hands-on Session #1 | **Bernhard Kretz**<br>Hands-on Session #2 |
| 14:40 | 15:10 | Coffee Break | Coffee Break | Coffee Break |
| 15:10 | 16:00 | **Miha Gunde**<br>Hands-on Session #2 | **Juraj Ovčar**<br>Hands-on Session #2 | **Bernhard Kretz**<br>Hands-on Session #3 |
| 16:00 | 16:05 | Break | Break | Break |
| 16:05 | 17:00 | **Miha Gunde**<br>Hands-on Session #3 | **Juraj Ovčar**<br>Hands-on Session #3 | Panel discussion |

## 3. Setting Up

> [!IMPORTANT]
> If you encounter any issues while setting up, check the [Issues related to the VM](#51-issues-related-to-the-vm).
>
> If you still cannot resolve your error, send a message through the school Slack channel
> or send an email to [Juraj.Ovcar@irb.hr](mailto:Juraj.Ovcar@irb.hr).

***
### 3.1. Installing the VM

 1. Download and install [VirtualBox](https://www.virtualbox.org/wiki/Downloads).
 2. Download the [**ams.ova**](https://drive.google.com/file/d/1Y1tqZAUKk5ZAjnA5kplyvTa9Xvbs9jd_/view?usp=drive_link) file.
 3. Run VirtualBox and click on ``File`` &rarr; ``Import Appliance...``
 4. Under ``File:``, enter the path to the downloaded **ams.ova** file.
 5. Click on ``Settings`` and make sure to specify the ``Machine Base Folder:`` path in a drive that has at least **20GB of free space**.
 6. Click on ``Finish``.

> [!WARNING]
> If you are using a Mac computer with an M-series chip (M1-M4), you should install
> the [macOS / Apple Silicon hosts](https://download.virtualbox.org/virtualbox/7.2.4/VirtualBox-7.2.4-170995-macOSArm64.dmg) platform package in step 3.
>
> You might still experience issues with the VM.
> In that case, send a message through the school Slack channel.

If everything goes well, read the following subsection.

***
### 3.2. About the VM

- The VM is started by selecting **ams** and clicking on ``Start``.

- The user credentials are:

    - username: ``atom``
    - password: ``ams``

- Croatian and English (US) keyboard layouts are installed in the VM, with the Croatian being the default. Switch between them by clicking the flag in the bottom right corner.

- If you wish to do so, when the VM is powered off, right-click on **ams** and select ``Settings`` to modify hardware resource allocation to the VM

***

> [!IMPORTANT]
> This is an [Arch Linux](https://en.wikipedia.org/wiki/Arch_Linux) VM.
> If you don't have previous experience with Linux,
> we encourage you to start the VM, log in, press ``Ctrl+Alt+T`` to activate the Linux terminal
> and spend some time using the [command line](https://linuxcommand.org/).
>
> There are [many great beginner guides](https://labex.io/linuxjourney) available online. Do not be afraid to experiment, since this is a **virtual** machine and you cannot do any harm to your computer. On the other hand, if you mess up something in the VM, you can always just reinstall it.

### 3.3. Local Installation

If you wish to install the software and run the tutorials on your own machine, read the chapter on [Software Installation](#6-software-installation).

## 4. Running the Exercises

TODO

## 5. Common Issues

### 5.1. Issues Related to the VM

- I'm getting an error while installing VirtualBox.

    - Try to resolve your error using the detailed [installation guide for VirtualBox](https://www.virtualbox.org/manual/ch02.html).

- While importing the **ams.ova** file, I'm getting a BIOS-related error, e.g. ``VT-x is disabled in the BIOS for all CPU modes``.

    - You need to enable virtualization in BIOS. Check [Enable Virtualization on Windows](https://support.microsoft.com/en-us/windows/enable-virtualization-on-windows-c5578302-6e43-4b4b-a449-8ced115f58e1) or try to search for a similar guide depending on your hardware.

- While importing the **ams.ova** file, I'm getting ``E_INVALIDARG (0x80070057)``.

    - Make sure you have enough space on the drive on which you are trying to import the VM. Revisit [step 5 of the VM installation](#31-installing-the-vm).

- I'm trying to start the VM on a Linux host, but I'm getting a USB-related error or something like ``ERROR [COM]: aRC=E_ACCESSDENIED (0x80070005)``.

    - You need to add your user account to the [vboxusers group](https://askubuntu.com/questions/377778/how-to-add-users-to-vboxusers-to-enable-usb-usage).

- After starting the VM, I'm seeing a blank screen.

    - Try to [increase video memory](https://askubuntu.com/questions/1134892/ubuntu-18-04-20-04-lts-on-virtualbox-boots-up-but-black-login-screen).

### 5.2. Other Issues

This is a placeholder for issues that may be encountered during testing.

## 6. Software Installation

> [!Warning]
> This chapter is relevant only if you want to **install the software on your own machine.**
>
> Most likely, you **don't** need to do this and instead should be running the exercises
> using the [VM](#31-installing-the-vm).
***

In this section, we give only a concise list of the software installed in the VM since
detailed installation instructions for the majority of the listed software
is easily found online.

- [Arch Linux](https://archlinux.org/)
- [yay](https://github.com/Jguer/yay)
- [Brave](https://brave.com/)
- [Xfce](https://www.xfce.org/)
- [miniconda](https://www.anaconda.com/docs/getting-started/miniconda/main)
- [openmpi](https://www.open-mpi.org/)
- [FFTW](https://www.fftw.org/)
- [LAMMPS](https://www.lammps.org/) with some basic additional packages ``(KSPACE MANYBODY MOLECULE RIGID)``
- [OVITO](https://www.ovito.org/)
- [Quantum ESPRESSO](https://www.quantum-espresso.org/) (only ``pw.x``, ``pp.x`` and ``ph.x``)
- [Xcrysden](http://www.xcrysden.org/)
- CPU version of [MACE](https://github.com/ACEsuit/mace)

Installing the CPU version of MACE is probably the only non-standard installation:

```bash
pip install torch torchvision --index-url https://download.pytorch.org/whl/cpu
pip install mace-torch
```
