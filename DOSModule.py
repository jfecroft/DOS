import numpy as np
import scipy.constants
from math import pi
from collections import OrderedDict


# define physical constants
EH2K = scipy.constants.physical_constants["hartree-kelvin relationship"][0]
K2EH = scipy.constants.physical_constants["kelvin-hartree relationship"][0]
ATOMICUNITOFTIME = scipy.constants.physical_constants["atomic unit of time"][0]


def get_data(filen):
    """reads a file of the form 1st line heading columns
       2nd lines on values for those headings
       returns a nested dictionary
       converts values to float if possible"""
    open_file = open(filen, "r")
    lines = open_file.readlines()
    keys = lines[0].split()

    data_dict = OrderedDict()

    for line in lines[1:]:
        items = line.split()
        for i, item in enumerate(items):
            try:
                items[i] = float(item)
            except:
                pass
        sub_dict = dict(zip(keys, items))
        data_dict.update({sub_dict['sys']: sub_dict})
    return data_dict


def _read_data(dir_loc, inputfile, qnmax):
    """read in the data from precalculated rovibrational energy levels of the
    system files are og the form 'txtjXXXtxt' where XXX is the L quantum number
    starting at 0
    """
    bound_states = []  # bound states for different L quntum numbers
    for i in xrange(qnmax+1):
        file = dir_loc + inputfile.replace('j', 'j'+str(i))
        try:
            data = np.sort(np.loadtxt(file, usecols=(1,)))
        except IOError:  # maximum quantum number reached set qnmax and return
            qnmax = i-1
            return bound_states, qnmax
        if len(np.atleast_1d(data)) == 1:
            data = np.reshape(data, 1)  # accounts for behavior of 0d arrays
        bound_states.append(data)
    return bound_states, qnmax


class AtomMoleculeDOS:
    def __init__(self, filen):
        self.systems = get_data(filen)

    def _compute_dos(self, J, MQN, bound_states_c, bound_states_d, vmax,
                     energy_range=10.0, energy_offset=0.0):
        """
        computes the dos for an atom moleucle collsiions as described in
        Statistical aspects of ultracold resonant scattering Mayle, Ruzic and
        Bohn
        Variables
        J       - total J quantum number.
        MQN     - MJ projection quantum number.
        bound_states_c - boundstates of the complex obtained from _read_data
        bound_states_d - boundstates of the molecule obtained from _read_data
        vmax    - maximum vibrational quantum number.
        energy_range  - dos = number of states counted/energy range
        energy_offset - calculate dos around energy = energy_offset
        """
        if abs(MQN) > abs(J):
            print 'physically impossible (abs(MJ) > abs(J))'
            return
        limit = (energy_range/2.0)*K2EH  # convert energy_range K to hartrees
        energy_offset *= K2EH
        outside = np.ma.masked_outside  # mask enties in list outside range
        count = np.ma.count  # counts the number of unmasked entries in a list
        abs_gs = bound_states_d[0][0]
        num = 0
        # looping over rotational state of dimer a
        for n in xrange(len(bound_states_d)):
            # looping over all l constant with J
            for l in xrange(len(bound_states_c)):
                # only include pairs which couple to form J
                if abs(n-l) <= J and n+l >= J:
                    # looping over all vibrational levels of a
                    for v in xrange(min(len(bound_states_d[n]), vmax+1)):
                        # degeneracy
                        d = len(xrange(max(-l, MQN-n), min(l, MQN+n)))+1
                        threshold_energy = bound_states_d[n][v]-abs_gs
                        num += count(outside(bound_states_c[l],
                                     -limit-threshold_energy+energy_offset,
                                     limit-threshold_energy+energy_offset))*d

        dos = (float(num)/energy_range)*1.0E-3
        lifetime = dos*1.0E3*EH2K*2.0*pi*ATOMICUNITOFTIME*1.0e9
        # return dos in mK-1 and lifetime in ns
        return dos, lifetime

    def _compute_num_open(self, J, MQN, bound_states_c, bound_states_d,
                          vmax, energy_offset=0.0):
        """
        computes the number of threshold channels of the complex with energy
        less than energy_offset.
        Statistical aspects of ultracold resonant scattering Mayle, Rizic and
        Bohn
        VAriables
        J       - total J quantum number.
        MQN     - MJ projection quantum number.
        bound_states_c - boundstates of the complex obtained from _read_data
        bound_states_d - boundstates of the molecule obtained from _read_data
        vmax    - maximum vibrational quantum number.
        energy_offset - calculate number of threshold channels of the complex
                        with energy less than energy_offset.
        """
        if abs(MQN) > abs(J):
            print 'physically impossible (abs(MJ) > abs(J))'
            return
        energy_offset *= K2EH
        abs_gs = bound_states_d[0][0]
        num_open = 0
        # looping over rotational state of dimer a
        for n in xrange(len(bound_states_d)):
            # looping over all l constant with J
            for l in xrange(len(bound_states_c)):
                # include pairs of rotational levels which couple to form J
                if abs(n-l) <= J and n+l >= J:
                    # looping over all vibrational levels of a
                    for v in xrange(min(len(bound_states_d[n]), vmax+1)):
                        threshold_energy = bound_states_d[n][v]-abs_gs
                        if threshold_energy < energy_offset:
                            num_open += 1
        return num_open

    def get_dos(self, key, J=0, MQN=0, nmax=100, vmax=9999, energy_offset=0.0):
        """
        simple wrapper around _read_data and AtomMoleculeDOS
        Variables
        J       - total J quantum number.
        MQN     - MJ projection quantum number.
        nmax    - maximum allowed rotaional quantum number of the molecule.
        vmax    - maximum allowed vibrational quantum number of the molecule.
        energy_offset - calculate dos around energy = energy_offset
        """
        lmax = nmax + J
        self.bound_states_d, self.nmax = _read_data(
            self.systems[key]['dimer_dirn'],
            self.systems[key]['dimer_filen'],
            nmax)
        self.bound_states_c, self.lmax = _read_data(
            self.systems[key]['cmplx_dirn'],
            self.systems[key]['cmplx_filen'],
            lmax)
        self.dos,  self.lifetime = AtomMoleculeDOS._compute_dos(
            self,
            J,
            MQN,
            self.bound_states_c,
            self.bound_states_d,
            vmax=vmax,
            energy_offset=energy_offset)
        return self.dos, self.lifetime

    def get_num_open(self, key, J=0, MQN=0, nmax=100, vmax=9999,
                     energy_offset=0.0):
        """
        simple wrapper around _read_data and _compute_num_open
        Variables
        J       - total J quantum number.
        MQN     - MJ projection quantum number.
        nmax    - maximum allowed rotaional quantum number of the molecule.
        vmax    - maximum allowed vibrational quantum number of the molecule.
        energy_offset - calculate dos around energy = energy_offset
        """
        lmax = nmax + J
        self.bound_states_d, self.nmax = _read_data(
            self.systems[key]['dimer_dirn'],
            self.systems[key]['dimer_filen'],
            nmax)
        self.bound_states_c, self.lmax = _read_data(
            self.systems[key]['cmplx_dirn'],
            self.systems[key]['cmplx_filen'],
            lmax)
        self.num_open = AtomMoleculeDOS._compute_num_open(
            self,
            J,
            MQN,
            self.bound_states_c,
            self.bound_states_d,
            vmax=vmax,
            energy_offset=energy_offset)
        return self.num_open


class MoleculeMoleculeDOS:
    def __init__(self, filen):
        self.systems = get_data(filen)

    def _compute_dos_MoleculeMolecule(self, J, MJ, nmax, vmax, lmax,
                                      bound_states_d, bound_states_c,
                                      energy_range=10.0):
        """
        computes the dos for an atom moleucle collsiions as described in
        Scattering of Ultracold Molecules in the Highly Resonant Regime --
        Mayle, Ruzic, Quememer and Bohn
        Variables
        J       - total J quantum number.
        MQN     - MJ projection quantum number.
        bound_states_c - boundstates of the complex obtained from _read_data
        bound_states_d - boundstates of the molecule obtained from _read_data
        vmax    - maximum vibrational quantum number.
        lmax    - maximum end-over-end rotational quantum number of the two
                  molecules
        nmax    - maximum rotational quantum number of a single molecule.
        energy_range  - dos = number of states counted/energy range
        """
        if abs(MJ) > abs(J):
            print 'physically impossible (abs(MJ) > abs(J))'
            return
        num = 0  # variable to hold the number of states between limits
        limit = (energy_range/2.0)*K2EH  # convert energy_range K to hartrees
        abs_gs = bound_states_d[0][0]  # energy of the absolute ground state
        for N in xrange(max(0, J-lmax), min(2*nmax, lmax-J)+1):
            for l in xrange(abs(J-N), min(lmax, J+N)+1):
                for n1 in xrange(max(0, N-nmax), nmax+1):
                    for n2 in xrange(abs(N-n1), min(N+n1, nmax)+1):
                        # looping over all vibrational levels of dimer1
                        for v1 in xrange(min(len(bound_states_d[n1]), vmax+1)):
                            # looping over all vibrational levels of dimer2
                            for v2 in xrange(min(len(bound_states_d[n2]),
                                             vmax+1)):
                                threshold_energy = (bound_states_d[n1][v1] -
                                                    abs_gs +
                                                    bound_states_d[n2][v2] -
                                                    abs_gs)
                                if (bound_states_c[l][0] >
                                        limit-threshold_energy):
                                    # the lowest state is higher than highest
                                    # threshold
                                    break
                                else:
                                    start_end = np.searchsorted(
                                        bound_states_c[l],
                                        [-limit-threshold_energy,
                                         limit-threshold_energy],
                                        'left')
                                    d = abs(max(-l, MJ-N)) + min(l, MJ+N) + 1
                                    num += d*(start_end[1] - start_end[0])
        # return dos in per uK i.e.times 10**-6K
        dos = (float(num)/energy_range)*1.0E-6
        # return lifetime in ms
        lifetime = dos*1.0E6*EH2K*2.0*pi*ATOMICUNITOFTIME * 1.0e3
        return dos, lifetime

    def get_dos_MoleculeMolecule(self, key, J=0, MQN=0, nmax=100, vmax=999):
        """
        simple wrapper around _read_data and MoleculeMoleculeDOS
        Variables
        J       - total J quantum number.
        MQN     - MJ projection quantum number.
        nmax    - maximum allowed rotaional quantum number of the molecule.
        vmax    - maximum allowed vibrational quantum number of the molecule.
        energy_offset - calculate dos around energy = energy_offset
        """
        lmax = 2*nmax+J
        self.bound_states_d, self.nmax = _read_data(
            self.systems[key]['dimer_dirn'],
            self.systems[key]['dimer_filen'],
            nmax)
        self.bound_states_c, self.lmax = _read_data(
            self.systems[key]['cmplx_dirn'],
            self.systems[key]['cmplx_filen'],
            lmax)
        self.dos, self.lifetime = MoleculeMoleculeDOS._compute_dos_MoleculeMolecule(
            self, J, MQN, self.nmax, vmax, self.lmax,
            self.bound_states_d, self.bound_states_c)
        return self.dos, self.lifetime
