import os
import sys
import argparse

import numpy as np
import matplotlib.pyplot as plt

def define_parser():
    """
    Defines a simple parser

    Returns
    -------
    parser : argparse.ArgumentParser
        A simple parser

    """

    parser = argparse.ArgumentParser()

    parser.add_argument('dos_filename', help="path to file containing the DOS", type=str)
    parser.add_argument('-o', '--output_filename',
                        help="if given, saves plot to a file",
                        type=str)

    return parser

def read_columns(filename, usecols=None, skiprows=0):
    """
    Reads columns from a text file and stores them in numpy arrays

    Parameters
    ----------
    filename : str
        Path to the file to be read
    columns : None | list of int
        Which columns to read with 0 being the first; e.g. (1, 2).
        If not given, all columns are read.
    skiprows : int
        How many rows to skip, defaults to 0

    Returns
    -------
    data : numpy.ndarray
        Data read from the text file to numpy array

    """

    if not os.path.exists(filename):
        raise FileNotFoundError(f'Cannot read file {filename}; the file does not exist.')

    data = np.loadtxt(fname=filename, usecols=usecols, skiprows=skiprows)

    return data

def get_fermi_energy_dos(filename, metal=True):
    """
Reads the Fermi energy

    There are two cases:

        1) metal
            the Fermi energy is the highest occupied state energy;
            can be read from the first line of DOS file

        2) insulator
            the Fermi energy is in the middle of the gap;
            the energy needs to be calculated

    Parameters
    ----------
    filename : str
        Path to the DOS file
    metal : bool
        Whether the system is a metal or not

    Returns
    -------
    ef : float
        Fermi energy

    """

    if not os.path.exists(filename):
        raise FileNotFoundError(f'Cannot read file {filename}; the file does not exist.')

    with open(filename, 'r') as f:
        first_line = f.readline()

    if 'EFermi' not in first_line:
        raise ValueError(f'Cannot read Fermi energy from {filename}.')

    # Nice, safe line
    ef = float(first_line.split('EFermi')[-1].strip().split(' ')[-2])

    if not metal:

        # Inefficient to read again, but fine for the school...
        dos_data = read_columns(filename, skiprows=1)
        energy, dos = dos_data[:, 0], dos_data[:, 1]

        # Separate DOS and energy at what QE thinks is the Fermi level
        dos_c = dos[np.where(energy > ef)]
        dos_v = dos[np.where(energy < ef)][::-1]

        energy_c = energy[np.where(energy > ef)]
        energy_v = energy[np.where(energy < ef)][::-1]
        
        # Actual conduction band index
        i_c = np.where(dos_c > 1e-8)
        ec = energy_c[i_c][0]

        # Actual valence band index
        i_v = np.where(dos_v > 1e-8)
        ev = energy_v[i_v][0]

        ef = (ev + ec) / 2.

    return ef


if __name__ == '__main__':
    """
    Plots energy vs. DOS.

    Usage examples:

        python plot_dos.py Si_dos.dat
        python plot_dos.py -o dos_figure.png Si_dos.dat

    """

    # Initialize parser
    parser = define_parser()
    args = parser.parse_args()
    
    # Read the DOS data
    dos_data = read_columns(args.dos_filename, skiprows=1)

    energy = dos_data[:, 0]
    dos = dos_data[:, 1]

    # Get Fermi energy
    ef = get_fermi_energy_dos(args.dos_filename, metal=False)

    # Set Fermi energy as the 0 of energy 
    energy -= ef

    # Plot
    fig, ax = plt.subplots(figsize=(8.0, 6.0))

    # Fermi energy
    ax.axvline(0, ls='--', color='black')

    # Plot DOS
    ax.fill_between(x=energy, y1=dos, where=energy<=0,
                    edgecolor='black', linewidth=1.1,
                    facecolor=(0, 0, 1, 0.5))

    ax.fill_between(x=energy, y1=dos, where=energy>0,
                    edgecolor='black', linewidth=1.1,
                    facecolor=(0, 0, 1, 0.1))

    # Plot limits
    ax.set_ylim(0)

    # Axis labels
    ax.set_xlabel(r'$E-E_\text{F}$ (eV)', fontsize=18)
    ax.set_ylabel(r'$g$ ($\text{eV}^{-1}$)', fontsize=18)

    # Ticks
    ax.minorticks_on()
    ax.tick_params(labelsize=14)

    if args.output_filename:
        plt.savefig(args.output_filename, format='png')

    else:
        plt.show()
