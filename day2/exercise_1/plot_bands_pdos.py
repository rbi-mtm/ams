import numpy as np

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as tck

def get_HSP(HSP_labels, bands_output_filename):
    """
    Get bands high-symmetry points (HSPs)
        labels and x-coordinates

    The HSP labels are inputted by the user.

    The x-coordinates are read from the output of the
        `bands.x` program

    An error is raised if the number of HSP labels is
        not equal to the number of x-coordinates.

    Parameters
    ----------
    hsp_labels : list of str
        HSP labels
    bands_output_filename : str
        Path to the output of the `bands.x` program

    Returns
    -------
    HSP : list of dict
        List of {label: x-coordinate} pairs

    Raises
    ------
    ValueError
        If the number of HSP labels is different
        than the number of found x-coordinates

    """

    HSP = []

    with open(bands_output_filename, 'r') as f:

        for line in f:

            if ('high-symmetry point:' in line):
                HSP_label = HSP_labels[len(HSP)]
                HSP_x = float(line.split()[-1])

                HSP.append({HSP_label: HSP_x})

            if len(HSP) == len(HSP_labels):
                break

        else:
            raise ValueError('Invalid number of x-coordinates of HSPs found.')

    return HSP

def get_number_of_kpts(filename):
    """
    Get number of k-points from a *.dat.gnu file

    Parameters
    ----------
    filename : str
        Path to the *.dat.gnu file

    Returns
    -------
    nk : int
        Number of k-points

    """

    # get number of k points

    with open(filename, mode='r') as fd:

        nk = 0

        lines = fd.readlines()

        for line in lines:

            if len(line.strip()) == 0:
                break

            nk += 1

    return nk

def load_bands(filename, e_fermi=0.0):
    """
    Load k-points and bands from a *.dat.gnu file

    Parameters
    ----------
    filename : str
        Path to the *.dat.gnu file
    e_fermi : float
        This value will be subtracted from each band

    Returns
    -------
    k : numpy.ndarray
        Array of k-point x-coordinates
    bands : numpy.ndarray
        Arrays of bands in eV

    """

    # Get number of k-points
    nk = get_number_of_kpts(filename)

    # Load data
    data = np.loadtxt(filename)

    k = data[:, 0][0:nk]
    bands = np.reshape(data[:, 1], (-1, len(k)))

    bands -= e_fermi

    return k, bands

def load_pdos(pdos_files, e_fermi=0.0):
    """
    Load PDOS(E)

    `pdos_files` should be a dictionary,
        whose (key, value) pairs correspond to
        (label, filename)

    The labels are arbitrary and will be used for plotting
    
    Parameters
    ----------
    pdos_files : dict of str
        Pairs of (arbitrary labels,
        filenames obtained from `sumpdos.x`)
    e_fermi : float
        This value will be subtracted from the PDOS energies

    Returns
    -------
    pdos : dict of numpy.ndarray
        Dictionary of labels and [energy, PDOS] arrays
    
    """
    
    pdos = {}

    for label, filename in pdos_files.items():

        pdos[label] = np.loadtxt(filename, skiprows=1)
        pdos[label][:, 0] -= e_fermi

    return pdos

def load_dos(filename, e_fermi=0.0):
    """
    Loads the total DOS

    Parameters
    ----------
    filename : str
        Path to file containing the total DOS;
        (obtained via projwfc.x)
    e_fermi : float
        This value will be subtracted from the DOS energies

    Returns
    -------
    dos : np.ndarray
        Array of (energy, DOS) values

    """
    
    dos = np.loadtxt(filename, skiprows=1, usecols=[0,1])
    dos[:, 0] -= e_fermi

    return dos

def plot_bands_pdos(HSP, k, bands, pdos, total_dos, savefig=False, ylim=None):
    """
    Plot bands, PDOS and total DOS

    Parameters
    ----------
    HSP : list of dict
        List of {label: x_coordinate} pairs
    k : numpy.ndarray
        Array of k-point x-coordinates
    bands : numpy.ndarray
        Arrays of bands in eV
    pdos : dict of numpy.ndarray
        Dictionary of labels and [energy, PDOS] arrays
    total_dos : np.ndarray
        Array of (energy, DOS) values
    savefig : bool
        If True, saves the figure to Si_bands_pdos.png
    ylim : None | tuple
        If given, sets the y-limits

    """

    gs = gridspec.GridSpec(1, 2,
                           width_ratios=[2.5, 1],
                           wspace=0.05)

    fig = plt.figure(figsize=(12.0, 6.0))

    ax_bands = fig.add_subplot(gs[0])
    ax_dos = fig.add_subplot(gs[1])

    # Plot HSPs

    hsp_labels, hsp_coords = [], []

    for hsp in HSP:

        label, coord = list(hsp.items())[0]

        hsp_labels.append(label)
        hsp_coords.append(float(coord))

        ax_bands.axvline(coord, linewidth=0.5, color='k', alpha=0.5)

    ax_bands.set_xticks(ticks=hsp_coords)
    ax_bands.set_xticklabels(hsp_labels, fontsize=14)
    ax_bands.xaxis.set_tick_params(bottom=False)

    ax_bands.set_xlim(min(k), max(k))

    # Plot bands
    import matplotlib.patheffects as pe

    for i, band in enumerate(bands):
        ax_bands.plot(k, band,
                      linewidth=3.0,
                      alpha=1.0,
                      color='cornflowerblue',
                      path_effects=[pe.Stroke(linewidth=7, foreground='black'), pe.Normal()])


    ax_bands.yaxis.set_minor_locator(tck.AutoMinorLocator())
    ax_bands.yaxis.set_tick_params(labelsize=12, direction='in', which='both')
    ax_bands.set_ylabel(r'$E-E_{F}$ (eV)', fontsize=14)

    # Set ylim

    if ylim:
        ax_bands.set_ylim(ylim)

    # Plot total DOS

    ax_dos.fill_betweenx(y=total_dos[:, 0],
                         x1=total_dos[:, 1],
                         facecolor='gray',
                         alpha=0.6,
                         label='Total DOS')

    # Plot PDOS
    # Dumb code but works...

    labels = pdos.keys()
    pdos_projs = list(pdos.values())

    for i, label in enumerate(labels):

        if '(s)' in label:
            label_s = label
            pdos_proj_s = pdos_projs[i]

        elif '(p)' in label:
            label_p = label
            pdos_proj_p = pdos_projs[i]

    ax_dos.fill_betweenx(y=pdos_proj_s[:, 0],
                         x1=0.0,
                         x2=pdos_proj_s[:, 1],
                         facecolor='mediumseagreen',
                         label=label_s,
                         edgecolor='black')

    ax_dos.fill_betweenx(y=pdos_proj_p[:, 0],
                         x1=pdos_proj_s[:, 1],
                         x2=pdos_proj_s[:, 1] + pdos_proj_p[:, 1],
                         facecolor='salmon',
                         label=label_p,
                         edgecolor='black')

    # Make the DOS ax look nicer

    ax_dos.tick_params(labelleft=False,
                       left=False)

    ax_dos.xaxis.set_tick_params(labelsize=12, direction='in', which='both')
    ax_dos.set_xlabel(r'PDOS (eV$^{-1}$)', fontsize=14)

    ax_dos.set_ylim(ax_bands.get_ylim())
    ax_dos.set_xlim(0.0, np.max(total_dos[:, 1]))

    ax_dos.xaxis.set_minor_locator(tck.AutoMinorLocator())

    ax_dos.legend(frameon=False, fontsize=14)

    # Plot Fermi energy

    ax_bands.axhline(0.0, ls='--', color='black', alpha=0.5)
    ax_dos.axhline(0.0, ls='--', color='black', alpha=0.5)
    
    # Show or save figure
    if savefig:
        plt.savefig('Si_bands_pdos.png', format='png', dpi=300)
    else:
        plt.show()

if __name__ == '__main__':
    """
    Plot bands, PDOS and total DOS

    """

    # Define Fermi energy in eV
    e_fermi = 6.328

    # Get bands high-symmetry points (HSPs) x-coordinates
    HSP = get_HSP(HSP_labels=[r'$L$',
                              r'$\Gamma$',
                              r'$X$',
                              r'$K$',
                              r'$\Gamma$'],
                  bands_output_filename='06_si_pp_bands.out')

    # Read the k-points and bands from the Si.bands.dat.gnu file
    k, bands = load_bands(filename='Si.bands.dat.gnu', e_fermi=e_fermi)

    # Read PDOS
    pdos = load_pdos(pdos_files={'Si (s)': 'Si.pdos_s.dat',
                                 'Si (p)': 'Si.pdos_p.dat'},
                     e_fermi=e_fermi)

    # Read total DOS
    total_dos = load_dos(filename='./pdos/Si.pdos.dat.pdos_tot',
                         e_fermi=e_fermi)

    # Plot bands and PDOS
    plot_bands_pdos(HSP=HSP,
                    k=k,
                    bands=bands,
                    pdos=pdos,
                    total_dos=total_dos,
                    savefig=False,
                    ylim=None)
