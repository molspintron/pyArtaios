import numpy as np
from pyArtaios_utilities import read_xyz, distance, partition, prepare_input_gaussian, write_transportin, prepare_input_turbomole
from pyArtaios_subsys import write_molden
from pyArtaios_localtransmission import load_t, plot_local, plot_local_interactive
import subprocess as sb
import os
import cclib

class artaios:

    def __init__(self, natoms, atoms, coords, input_settings):
        '''initializes the artaios class'''

        # initialize switches
        self._ready = False
        self._partitioned = False

        # initialize empty results dictionary
        self.results = {}

        # initialize default settings
        self.settings = {}
        self._default_settings()

        # replace settings with input settings
        for key in input_settings.keys():
            if key in self.settings.keys():
                self.settings[key] = input_settings[key]

        # assign molecule
        self.natoms = natoms
        self.atoms  = atoms
        self.coord  = coords

        # partitioning
        self.left = []
        self.central = []
        self.right = []

        try:
            self.left, self.central, self.right = partition(natoms, atoms, coords)
            self._partitioned = True
        except IndexError:
            pass
            # print('automatic partitioning failed, please define manually')

        # check for consistency
        self._check()


    def calculate(self):
        '''performs all necessary calculations'''

        if not self._ready:
            self._check()

        if self._ready:

            # clear old results
            if len(self.results) > 0:
                self.results = {}

            # create new folder and move old stuff to backup folder
            try:
                os.mkdir('artaios')
            except:
                try:
                    os.rename('artaios', 'artaios_backup')
                    os.mkdir('artaios')
                except:
                    i = 1
                    while True:
                        if not os.path.isdir('artaios_backup_'+str(i)):
                            break

                        i += 1

                    os.rename('artaios', 'artaios_backup_'+str(i))
                    os.mkdir('artaios')

            os.chdir('artaios')

            # perform necessary single point calculation
            print('perform single point calculation...', flush = True)
            self._singlepoint()

            # get mo energies and approx. Fermi level from junction
            try:
                self._get_mo_energies()
            except NotImplementedError:
                print('ERROR: Reading of MO energies with the chosen settings is not supported yet.')

            # run artaios for normal transport calculationi
            print('perform artaios calculation...', flush = True)
            self._artaios_transport()

            # parse local transmissionsi
            print('get local transmissions...', flush = True)
            try:
                self._get_local_transmission()
            except OSError:
                print('ERROR: File for local transmissions not found. Maybe using an old version of ARTAIOS?')
            except NotImplementedError:
                print('ERROR: Multiplicity of 2 for local transmissions not yet supported')

            # run artaios for subsystem MOs
            print('get subsystem MOs...', flush = True)
            try:
                self._get_subsys_energies()
            except NotImplementedError:
                print('ERROR: Calculation of subsystem MOs with the chosen settings is not supported yet.')

            # parse outputs
            os.chdir('..')
            # except:
            #     pass
                #   print('Error: Folder "artaios" already exists')

        else:
            print('class is not ready, please check settings')


    def _artaios_transport(self):
        '''prepares and runs artaios calculation
        to get the transmission
        '''

        # prepare input file
        write_transportin(self.settings, self.left, self.central, self.right)

        # run artaios
        with open('transport.out', 'w') as file:
            sb.call(['artaios', 'transport.in'], stdout = file)

        # get results
        temp = []
        temp.append(np.genfromtxt('transmission.1.dat', usecols = (0)))
        temp.append(np.genfromtxt('transmission.1.dat', usecols = (1)))

        if self.settings['multi'] == 1:
            temp.append(0*np.genfromtxt('transmission.1.dat', usecols = (1)))
           # self.results['transmission'] = np.genfromtxt('transmission.1.dat')
        else:
            temp.append(np.genfromtxt('transmission.2.dat', usecols = (1)))

        self.results['transmission'] = np.array(temp).T


    def _get_subsys_energies(self):
        '''performs subsys calculations with artaios'''

        # check for supported methods
        if self.settings['multi'] == 1 and (self.settings['qcprog'].lower() == 'gaussian' or self.settings['qcprog'].lower() == 'turbomole'):
            # gaussian
            if self.settings['qcprog'].lower() == 'gaussian':
	        # create an extra folder and copy necessary files
                os.mkdir('subsys')
                sb.call(['cp', 'transport.log', 'subsys/transport.log'])
                sb.call(['cp', 'hamiltonian.1', 'subsys/hamiltonian.1'])
                sb.call(['cp', 'overlap', 'subsys/overlap'])
                sb.call(['cp', 'transport.in', 'subsys/transport.in'])

                # go to folder
                os.chdir('subsys')

                # get molden file from gaussian log
                data = cclib.io.ccread('transport.log')
                write_molden(data.natom, data.atomnos, data.atomcoords[-1], 'molden.input', data.gbasis, data.moenergies, data.homos, data.mocoeffs)

            # turbomole
            if self.settings['qcprog'] == 'turbomole':
                p = sb.Popen(['tm2molden'], stdin=sb.PIPE, stdout = sb.DEVNULL)
                p.stdin.write('\n'.encode())
                p.stdin.write('\n'.encode())
                temp = p.communicate()

            # load transport.in
            with open('transport.in', 'r') as file:
                lines = file.readlines()

            # change energy range
            for i in range(len(lines)):
                if 'steps' in lines[i]:
                    lines[i] = ' steps 1\n'

            # add lines for subsystem MO calculation
            lines.append('$subsystem\n')
            lines.append(' print_molden\n')
            lines.append(' do_diag_central\n')
            lines.append(' moldeninfile molden.input\n')
            lines.append('$end\n')

            with open('subsys.in', 'w') as file:
                for l in lines:
                    file.write(l)

            # run artaios
            with open('subsys.out', 'w') as file:
                sb.call(['artaios', 'subsys.in'], stdout = file)

            # get energies
            energies = []
            with open('central.molden_1', 'r') as file:
                lines = file.readlines()

            for l in lines:
                if 'Ene=' in l:
                    energies.append(float(l.split()[1]))

            self.results['subsys energies'] = energies

            if self.settings['qcprog'] == 'gaussian':
                os.chdir('..')

        else:
            raise NotImplementedError


    def _get_local_transmission(self):
        '''get local transmissions'''

        if self.settings['multi'] == 1:
            self.results['local transmissions'] = [*load_t('localtrans.dat')]
        else:
            self.results['local transmissions'] = [*load_t('localtrans.dat')]


    def _get_mo_energies(self):
        '''parses QC output for MO energies'''

        if self.settings['qcprog'].lower() == 'gaussian':
            output = cclib.io.ccread('transport.log')
        elif self.settings['qcprog'].lower() == 'turbomole':
            raise NotImplementedError
            if self.settings['ri']:
                output = cclib.io.ccread('ridft.out')
            else:
                output = cclib.io.ccread('dscf.out')
        else:
            pass

        if output.closed_shell:
            alpha = np.zeros([len(output.moenergies[0]), 2])
            alpha[:, 0] = output.moenergies[0]
            for i in range(output.homos[0]+1):
                alpha[i, 1] = 2

            self.results['mo energies'] = [alpha]

        else:
            alpha = np.zeros([len(output.moenergies[0]), 2])
            alpha[:, 0] = output.moenergies[0]
            for i in range(output.homos[0]+1):
                alpha[i, 1] = 1

            beta = np.zeros([len(output.moenergies[1]), 2])
            beta[:, 0] = output.moenergies[1]
            for i in range(output.homos[1]+1):
                beta[i, 1] = 1

            self.results['mo energies'] = [alpha, beta]


    def plot_local(self, input_energy = -5.0, spin = 1, threshold = 0.1, labels = True):
        '''Displays local transmission for the given energy

        Wrapper for the plot_local routine from
        pyArtaios_localtransmission.py

        Parameters:
        -----------
        input_energy : energy for which to plot the local
            transmission, default is -5.0
        threshold : threshold for which the local transmission
            are drawn. This is given in relative number with
            respect to the highest local transmission
        labels : switch to display the value of
            the local transmission at the arrows, default is
            True

        '''

        # clip spin
        spin = int(np.clip(spin, 1, 2))

        # choose correct local transmission
        local_trans = self.results['local transmissions'][spin-1]

        # check whether chosen energy is available, or choose closest
        if round(float(input_energy), 2) not in local_trans.keys():
            energies = np.fromiter(self.results['local transmission'].keys(), dtype=float)
            ind = np.argmin(abs(energies - round(float(input_energy),2)))
            energy = energies[ind]
        else:
            energy = round(float(input_energy), 2)

        print('Showing local transmissions for '+str(energy)+' eV...')
        # plot_local(self.atoms, self.coord, local_trans[energy], tlabels = labels, threshold = threshold, threshold_method = 'relative')

        plot_local_interactive(self.atoms, self.coord, self.results['local transmissions'])

    def _singlepoint(self):
        '''performs single point calculations'''

        if self.settings['qcprog'].lower() == 'gaussian':
            # gaussian 09 calculation
            # prepate input
            prepare_input_gaussian(self.natoms, self.atoms, self.coord, self.settings)

            # run programs
            sb.call(['nohup', 'g09', 'transport.com'])
            sb.call(['g09_2unform', 'transport.log', str(self.settings['multi'])])

        elif self.settings['qcprog'] == 'turbomole':
            # turbomole calculation
            prepare_input_turbomole(self.natoms, self.atoms, self.coord, self.settings)

            # run turbomole
            if self.settings['ri']:
                with open('ridft.out', 'w') as file:
                    sb.call(['ridft'], stdout = file)
            else:
                with open('dscf.out', 'w') as file:
                    sb.call(['dscf'], stdout = file)

            # extract overlap and hamiltonien
            if self.settings['multi'] == 1:
                sb.call(['tm2unformcl'])
            else:
                sb.call(['tm2unformop'])


    def _default_settings(self):
        '''initialize default settings'''

        self.settings['fermi'] = -5.0		# fermi energy, e.g. for transport or iets calculations
        self.settings['functional'] = 'bp86'	# exc functional
        self.settings['basis set'] = 'lanl2dz'	# basis set
        self.settings['multi'] = 1		# multiplicity
        self.settings['n cores'] = 1		# number of cores
        self.settings['nosymmetry'] = False	# nosymmetry keyword for gaussian
        self.settings['charge'] = 0		# charge of the system
        self.settings['qcprog'] = 'gaussian'    # qc program used for calculation
        self.settings['method'] = 'dft'         # type of qc calculation
        self.settings['ri']     = False         # switch for density fitting

        # settings for transport calculations
        self.settings['eupp']   = 0             # upper energy limit for transmission function
        self.settings['elow']   = -8	        # upper energy limit for transmission function
        self.settings['esteps'] = 400           # number of steps for transmission function
        self.settings['print green'] = False
        self.settings['read green']  = False
        self.settings['gsize']  = 0	        # size of greens matrix
#        self.settings['spin']   = 1     # spin


    def get(self, key):
        '''getter function for results'''
        if key in self.results.keys():
            return self.results[key]
        else:
            print(str(key)+' not available, did you started the calculation?')


    def _check(self, verbose = False):
        '''consistency check

        Performs checks for a bunch of method combinations and
        sets the `ready` switch accordingly
        '''

        self._ready = True

        # list of suppported basis sets for turbomole
        basis_sets_turbomole = ['def2-SVP', 'def2-TZVP', 'def2-SVPP', 'def2-TZVPP']

        # is the program supported:
        if self.settings['qcprog'].lower() == 'gaussian' or self.settings['qcprog'].lower() == 'turbomole':
            pass
        else:
            print('QC program "'+self.settings['qcprog']+'" not supported, please use "gaussian" or "turbomole"')
            self._ready = False

        # check for method
        if self.settings['method'].lower() != 'dft':
            print('computational method not supported, please use DFT')
            self._ready = False

        if self.settings['qcprog'] == 'turbomole':
            if not self.settings['basis set'] in basis_sets_turbomole:
                print('basis set not supported, please use one of those: '+ ', '.join(basis_sets_turbomole))
                self._ready = False

        # check partitioning
        if len(self.left) + len(self.right) + len(self.central) != self.natoms:
            print('partitioning looks wrong, please check')
            self._ready = False

        # check for double entries in the partitioning
        if any(i in self.left for i in self.central) or any(i in self.left for i in self.right) or any(i in self.right for i in self.central):
            print('partitioning looks wrong, entries appear multiple times')
            self._ready = False
