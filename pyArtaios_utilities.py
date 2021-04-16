import numpy as np
import cclib
import subprocess as sb

def write_coord(natoms, atoms, coords, freeze = [], filename='coord'):
    '''writes a coord file

    Parameters:
    -----------
    natoms: number of atoms, integer
    atoms: atom symbols, list
    coords: coordinates, numpy array
    freeze: list containing the indices of the atoms to be freezed
    filename: the output filename (default = coord)
    '''

    convert = 1.88972612 # angstrom/5.29177211×10^-9 centimeters
    coordfile = open(filename, 'w')
    coordfile.write('$coord\n')

    for i in range(natoms):
        if i in freeze:
            string = '  {0: 2.6f} {1: 2.6f} {2: 2.6f} {3:^2s} f\n'.format(*coords[i, :]*convert, atoms[i])
        else:
            string = '  {0: 2.6f} {1: 2.6f} {2: 2.6f} {3:^2s}\n'.format(*coords[i, :]*convert, atoms[i])

        coordfile.write(string)

    coordfile.write('$user-defined bonds\n')
    coordfile.write('$end\n')


def write_transportin(settings, left, central, right, filename = 'transport', subsys = False):
    '''Writes transport.in file according to settings and partitioning

    Parameters:
    -----------
    settings : dictionary contained the necessary options, look e.g.
        in artaios class definition or the readme for a list of all
        settings
    left : list of the atom numbers of the left electrode
    central : list of the atom numbers of the central region
    right : list of the atom numbers of the right electrode
    filename : filename for the gaussian log file, with the ".log"
    '''

    # Since python counts from 0 in lists, the partitioning lists have to increased by 1
    part_l = [x+1 for x in left]
    part_c = [x+1 for x in central]
    part_r = [x+1 for x in right]

    target = open('transport.in', 'w')

    # Block: general
    target.write('$general \n')
    target.write(' rbas\n')
    if settings['qcprog'].lower() == 'gaussian':
        target.write(' qcprog g09 \n')
        target.write(' mosfile '+filename+'.log \n')
    elif settings['qcprog'].lower() == 'turbomole':
        target.write(' qcprog tm\n')
        if settings['ri']:
            target.write(' mosfile ridft.out\n')
        else:
            target.write(' mosfile dscf.out\n')

    ## for transport
    if settings['job type'] == 'transport':
        target.write(' do_transport \n')
        target.write(' ham_conv 27.21 \n')
        if settings['print green']:
            target.write(' print_green\n')
        if settings['read green']:
            target.write(' read_green\n')
        target.write(' loewdin_central\n')

    ## for jgreen
    elif settings['job type'] == 'jgreen':
        target.write(' do_jgreen\n')
        target.write(' nelalpha '+str(settings['nelalpha'])+'\n')
        target.write(' nelbeta '+str(settings['nelbeta'])+'\n')
        target.write(' s_a '+str(settings['s_a'])+'\n')
        target.write(' s_b '+str(settings['s_b'])+'\n')

    target.write('$end \n')
    target.write('\n')

    # Block: partitioning
    target.write('$partitioning \n')
    target.write(' leftatoms '   +",".join(map(str, part_l))+'\n')
    if settings['job type'] == 'transport':
        target.write(' centralatoms '+",".join(map(str, part_c))+'\n')
    target.write(' rightatoms '  +",".join(map(str, part_r))+'\n')
    target.write('$end \n')
    target.write('\n')

    # Block: system
    target.write('$system \n')
    target.write(' nspin '+str(settings['spin'])+'\n')
    target.write('$end \n')
    target.write('\n')

    if settings['job type'] == 'transport':
        # Block: energy_range
        target.write('$energy_range \n')
        target.write(' start '+str(settings['elow'])+'\n')
        target.write(' end '  +str(settings['eupp'])+'\n')
        target.write(' steps '+str(settings['esteps'])+'\n')
        target.write('$end \n')
        target.write('\n')

        # Block: electrodes
        target.write('$electrodes \n')
        target.write(' fermi_level '+str(settings['fermi'])+'\n')
        target.write(' self_energy wbl \n')
        target.write(' dos_s 0.036 \n')
        target.write('$end \n')
        target.write('\n')

        if subsys:
            target.write('$subsystem\n')
            target.write(' print_molden\n')
            target.write(' do_diag_central\n')
            target.write(' moldeninfile molden.input\n')
            target.write('$end \n')
            target.write('\n')
        else:
            # Block: local transmissions
            target.write('$local_transmission\n')
            target.write(' bondflux\n')
            target.write('$end\n')
            target.write('\n')

    if settings['qcprog'] == 'turbomole':
        target.write('$debug\n')
        target.write(' xyzfile coord\n')
        target.write('$end\n')

    target.write('\n')
    target.close()

def read_xyz(filename):
    '''Reads an xyz file and returns the structure.

    Parameters:
    -----------
    filename : filename of a xyz file

    Returns:
    --------
    natoms : number of atoms
    atoms : list of atom symbols
    coords: numpy array of cartesian coordinates
    -
    '''

    atoms = []
    coord = []

    with open(filename) as xyz:
        n_atoms = int(xyz.readline())
        title = xyz.readline()
        for i in range(n_atoms):
            line = xyz.readline()
            atom,x,y,z = line.split()
            atoms.append(atom)
            coord.append([float(x), float(y), float(z)])

    xyz.close()
    coords = np.array(coord)

    return n_atoms, atoms, coords


def distance(coord1, coord2):
    '''calculate the distance between two points'''
    dx = coord1[0] - coord2[0]
    dy = coord1[1] - coord2[1]
    dz = coord1[2] - coord2[2]
    dist = np.sqrt(dx**2 + dy**2 + dz**2)
    return dist


def prepare_input_turbomole(natoms, atoms, coords, settings):
    '''sets up turbomole input by running define
    
    Parameters:
    -----------
    natoms : number of atoms
    atoms : list of atom symbols
    coords : numpy array containing the cartesian coordinates
    settings : dictionary with settings for the DFT calculation
    '''

    # write coord file
    # freezing atoms is not used
    write_coord(natoms, atoms, coords, freeze = [], filename='coord')

    # create commands for define
    # proper counting of '\n' is the key
    menu_geometry = ['\n', '\n', 'a coord\n', '*\n', 'no\n', ]
    menu_basis    = ['b\n', 'all '+settings['basis set']+'\n', '*\n']
    menu_guess    = ['eht\n', '\n', str(settings['charge'])+'\n', '\n']
    menu_methods  = ['dft\n', 'on\n', 'func\n', str(settings['functional'])+'\n', 'grid\n', 'm4\n', '\n','scf\n', 'conv\n', '8\n', 'iter\n', '300\n', '\n']
    if settings['ri']:
        menu_methods += ['ri\n', 'on\n', '\n']

    # assemble all commands
    define_commands = menu_geometry + menu_basis + menu_guess + menu_methods + ['*\n']
    # open define with a pipe for input
    # output is piped to DEVNULL/the void
    p = sb.Popen(['define'], stdin=sb.PIPE, stdout = sb.DEVNULL)

    # execute all commands and hope for the best
    for command in define_commands:
        p.stdin.write(command.encode())

    temp = p.communicate()


def prepare_input_gaussian(natoms, atoms, coords, settings, filename = 'transport'):
    '''
    writes a g09 input file fitting to artaios

    Parameters:
    -----------
    natoms : number of atoms
    atoms : list of atom symbols
    coords : numpy array containing the cartesian coordinates
    settings : dictionary with settings for the DFT calculation
    filename : filename for the .com file

    '''
    ad_string = ''
    if settings['nosymmetry']:
        ad_string += ' nosymmetry'

    if settings['ri']:
        ad_string += ' denfit'

    with open(filename+'.com', 'w') as file:
        # first header
        file.write('%chk='+filename+".chk\n")
        if 'gaussian header' in settings.keys():
            file.write(settings['gaussian header']+'\n')
        file.write("%NProcShared="+str(settings['n cores'])+'\n')
        file.write('#P '+settings['functional']+"/"+settings['basis set']+ad_string+' GFPrint\n\n')
        file.write('Transport calculation\n\n')
        file.write(str(settings['charge']) + ' ' + str(settings['multi']))
        file.write('\n')
        # first coords
        for i in range(natoms):
            string = '{0:^2s} {1: 2.6f} {2: 2.6f} {3: 2.6f}\n'.format(atoms[i], *coords[i, :])
            file.write(string)

        # second header
        file.write('\n\n')
        file.write("--Link1--\n")
        file.write("%chk="+filename+".chk\n")
        if 'gaussian header' in settings.keys():
            file.write(settings['gaussian header']+'\n')

        file.write("%NProcShared="+str(settings['n cores'])+'\n')
        file.write('#P '+settings['functional']+"/"+settings['basis set']+ad_string+' GFINPUT guess=read\n')
        file.write('# iop(6/7=3)\n')
        file.write('# iop(5/33=3)\n')
        file.write('# iop(3/33=1)\n\n')
        file.write('Transport calculation\n\n')
        file.write(str(settings['charge']) + ' ' + str(settings['multi']))
        file.write('\n')
        # second coords
        for i in range(natoms):
            string = '{0:^2s} {1: 2.6f} {2: 2.6f} {3: 2.6f}\n'.format(atoms[i], *coords[i, :])
            file.write(string)

        file.write('\n\n')



def partition(natoms, atoms, coords):
    '''
    partition a molecular junction for Artaios
    routine inspired/copied from artaios_define

    Parameters:
    -----------
    natoms : number of atoms
    atoms : list of atom symbols
    coords : numpy array containing the cartesian coordinates

    '''
    '''
    def distance(coord1, coord2):
        # returns distance between two coordinates
        dx = coord1[0] - coord2[0]
        dy = coord1[1] - coord2[1]
        dz = coord1[2] - coord2[2]
        dist = np.sqrt(dx**2 + dy**2 + dz**2)
        return dist
    '''

    # check wether two sulfur atoms are present
    sulfur_pos = []
    part_l = []
    part_r = []
    part_c = []

    # Trying to find two sulfur atoms
    for i in range(natoms):
        if atoms[i].lower() == 's':
            sulfur_pos.append(i)

    # User input to specify terminal sulfur atoms
    if len(sulfur_pos) != 2:
        raise IndexError('Cannot find two terminal sulfur atoms')
    else:
    #     s1 = -1
    #     s2 = -1
    #     while True:
    #
    #         print('')
    #         print(' Cannot find two terminal sulfur')
    #         print(' atoms, specification needed:')
    #         print(' ----------------------------------------------- ')
    #         for i in range(natoms):
    #             if i == s1:
    #                 print(' s1:', i+1, '{0:^2s} {1: 2.6f} {2: 2.6f} {3: 2.6f}\n'.format(atoms[i], *coords[i, :]))
    #             elif i == s2:
    #                 print(' s2:', i+1, '{0:^2s} {1: 2.6f} {2: 2.6f} {3: 2.6f}\n'.format(atoms[i], *coords[i, :]))
    #             else:
    #                 print('    ', i+1, '{0:^2s} {1: 2.6f} {2: 2.6f} {3: 2.6f}\n'.format(atoms[i], *coords[i, :]))
    #         print(' ----------------------------------------------- ')
    #         print(' Enter s1 <number> and s2 <number> to')
    #         print(' specify terminal atoms of the molecule.')
    #         print(' * to exit menu')
    #         command = raw_input(' :')
    #         # Todo: Check wether s1 and s2 is in range
    #         if command[:2]   == 's1':
    #             s1 = int(command[3:])-1
    #         elif command[:2] == 's2':
    #             s2 = int(command[3:])-1
    #         elif command == '*':
    #             if s1 != s2:
    #                 sulfur_pos = [s1,s2]
    #                 break
    #             else:
    #                 print(' ERROR: s1 and s2 cannot be the same atom!')
    #         else :
    #             print(' Unkown command')

    # partitioning junction
        for i in range(natoms):
            if atoms[i].lower() == 'au':
                dist_s1 = distance(coords[i, :], coords[sulfur_pos[0], :])
                dist_s2 = distance(coords[i, :], coords[sulfur_pos[1], :])
                if dist_s1 < dist_s2:
                    part_r.append(i)
                else:
                    part_l.append(i)
            else:
                part_c.append(i)

    return part_l, part_c, part_r
