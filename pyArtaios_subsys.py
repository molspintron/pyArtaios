import os
import numpy as np
import subprocess as sb

def atom_number2name(numbers):
    '''converts list of atom numbers to list of atom names'''
    names = []
    convert_dict = {1:'H', 2:"He", 3: "Li", 4:"Be", 5:"B", 6:"C", 7: "N", 8:"O", 9:"F", 10:"Ne", 11:"Na", 12:"Mg", 13:"Al", 14:"Si", 15:"P", 16:"S", 17:"Cl", 18:"Ar", 26: 'Fe', 27:'Co', 29:'Cu', 79:'Au'}
    for n in numbers:
        names.append(convert_dict[n])

    return names


def atom_name2number(names):
    '''converts list of atom names to list of atom numbers'''
    numbers = []
    convert_dict = {'H':1, "He":2, "Li":3, "Be":4, "B":5, "C":6, "N":7, "O":8, "F":9, "Ne":10, "Na":11, "Mg":12, "Al":13, "Si":14, "P":15, "S":16, "Cl":17, "Ar":18, 'Fe':26, 'Co':27, 'Cu':29, 'Au':79}
    for n in names:
        numbers.append(convert_dict[n])

    return numbers


def write_molden(natoms, atoms, coords, filename = 'test.molden', basis = [], mo = [], homo = [], mocoeff = []):
    '''Routine to write a molden file from gaussian output.

    The input data should be in the format supplied by cclib. If the
    basis is provided, basis set information will be written to the
    output file. If information about the mo energies, the homo/lumos
    and the mo coefficients are given, those will be written to the
    output as well.

    TODO: - add support for open shell systems
          - check [5D10F]
    '''

    # check whether i got atom numbers or atom names
    if isinstance(atoms[0], str):
        atom_names = atoms
        atom_numbers = atom_name2number(atoms)
    else:
        atom_names = atom_number2name(atoms)
        atom_numbers = atoms

    with open(filename, 'w') as file:
        # write header
        file.write('[Molden Format]\n')
        file.write('[Atoms] Angs\n')

        # write coordinates
        for i in range(natoms):
            file.write('{0:^2s} {1: 4.0f} {2: 3.0f} {3:>12.6f} {4:>12.6f} {5:>12.6f}\n'.format(atom_names[i], i+1, atom_numbers[i], coords[i, 0], coords[i, 1], coords[i, 2]))

        if basis != []:
            # write basis set
            if len(basis) > 0:
                # the following keyword is important for the ordering
                # of D and F functions - but i'm not sure whether this
                # is correct in all cases
                # http://cheminf.cmbi.ru.nl/molden/molden_format.html
                file.write('[5D10F]\n')
                file.write('[GTO]\n')
                for i in range(natoms):
                    # new atom
                    file.write('  '+str(i+1)+' 0\n')
                    for b in basis[i]:
                        # ao
                        file.write(' {0:^2s} {1: 3d} {2:1.2f}\n'.format(b[0].lower(), len(b[1]), 1.00))
                        # coefficients
                        for l in b[1]:
                            file.write(' {0: 2.10e} {1: 2.10e}\n'.format(*l))

                    file.write('\n')

        if mo != [] and len(homo) > 0 and mocoeff != []:
            # so far, only alpha electrons!
            convert = 0.036749322176
            if len(mo) > 0:
                file.write('\n[MO]\n')
                # write mo orbitals
                for i in range(len(mo[0])):
                    file.write(' {0:} {1: 4.4f}\n'.format('Ene=', mo[0][i]*convert))
                    file.write(' Spin= Alpha\n')
                    if i <= homo[0]:
                        file.write(' Occup=   2.000000\n')
                    else:
                        file.write(' Occup=   0.000000\n')

                    for j in range(len(mocoeff[0][i])):
                        file.write('{0: 4d}  {1: 4.6f}\n'.format(j+1, mocoeff[0][i][j]))
