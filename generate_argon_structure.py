"""Generate the structure for NEMD simulation"""

import numpy as np
import matplotlib.pyplot as plt

Unitcell = np.array([[0, 0, 0], [0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]])

unified_atomic_mass = 1.66e-27  #kg

sigma = 3.4

Lattice_Constant = sigma * 1.56  # 10K

mass1 = 39.948
mass2 = 39.948 *1

print("Lattice constant: {:.3f}".format(Lattice_Constant))


def duplicate(Nx, Ny, Nz, Unitcell, LC):
    """duplicate the system by Nx, Ny, Nz"""
    natom_in_cell = Unitcell.shape[0]
    ntotal_cell = Nx * Ny * Nz
    print("Total no. of cell:", ntotal_cell)
    print("Natom in cell:", natom_in_cell)
    out_data = np.zeros((natom_in_cell * ntotal_cell, 3))

    icell = -1
    cell_pos = np.zeros((natom_in_cell, 3))
    for ix in np.arange(Nx):
        for iy in np.arange(Ny):
            for iz in np.arange(Nz):
                icell += 1
                x1 = icell * natom_in_cell
                x2 = x1 + natom_in_cell

                cell_pos[:, 0] = ix
                cell_pos[:, 1] = iy
                cell_pos[:, 2] = iz
                atom_pos = (Unitcell + cell_pos) * LC
                out_data[x1:x2, :] = atom_pos
    return out_data


def write_data_lmp(data, Nx, Ny, Nz, LC, name="Ar"):
    """
    write the data in the xyz format
    """
    natom = data.shape[0]
    outfile = name + '_structure.lmp'
    with open(outfile, 'w') as f:
        # write header
        f.write("\n{:d} atoms\n".format(natom))
        f.write("1 atom types\n")
        f.write("{:.5f}\t{:.5f}\txlo xhi\n".format(0, Nx * LC))
        f.write("{:.5f}\t{:.5f}\tylo yhi\n".format(0, Ny * LC))
        f.write("{:.5f}\t{:.5f}\tzlo zhi\n".format(0, Nz * LC))
        f.write("\n")
        f.write("Masses\n\n")
        f.write("1 {:.2f} \n".format(mass1))
        #f.write("2 {:.2f} \n".format(mass1))
        #f.write("3 {:.2f} \n".format(mass1))
        #f.write("4 {:.2f} \n".format(mass2))
        #f.write("5 {:.2f} \n".format(mass2))
        #f.write("6 {:.2f} \n".format(mass2))
        f.write("Atoms\n\n")
        for ii in range(natom):
            fmt = "{:d}\t{:d}\t{:.5f}\t{:.5f}\t{:.5f}\n"
            atom_type = 1
            # if data[ii, 0] < 1 * LC:
            #     atom_type = 1
            # elif 1*LC <= data[ii, 0] < (11) * LC:
            #     atom_type = 2
            # elif 11*LC <= data[ii, 0] < (30) * LC:
            #     atom_type = 3
            # elif 30*LC <= data[ii, 0] < (49) * LC:
            #     atom_type = 4
            # elif 49*LC <= data[ii, 0] < (59) * LC:
            #     atom_type = 5
            # else:
            #     atom_type = 6
            f.write(
                fmt.format(ii + 1, atom_type, data[ii, 0], data[ii, 1],
                           data[ii, 2]))


def write_data_xyz(data, name="Ar"):
    """
    write the data in the xyz format
    """
    natom = data.shape[0]
    outfile = name + '_{}.xyz'.format(natom)
    with open(outfile, 'w') as f:
        f.write("{}\n".format(natom))
        f.write("\n")
        for ii in np.arange(natom):
            fmt = "{}\t" + "{:.5f}\t" * 3 + "\n"
            f.write(fmt.format(name, data[ii, 0], data[ii, 1], data[ii, 2]))


def main():

    nx = 8
    ny = 8
    nz = 8
    LC = Lattice_Constant

    data = duplicate(nx, ny, nz, Unitcell, LC)
    write_data_xyz(data)
    write_data_lmp(data, nx, ny, nz, LC)


if __name__ == "__main__":
    main()
