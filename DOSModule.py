"""
Module for comutation of densities of states for atom molecule collsions
"""
import numpy as np
import scipy.constants
from math import pi
from collections import OrderedDict

# pylint: disable=E1103
# pylint: disable=R0902
# pylint: disable=R0903
# pylint: disable=R0913,R0914

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
            except IOError:
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
        filen = dir_loc + inputfile.replace('j', 'j'+str(i))
        try:
            data = np.sort(np.loadtxt(filen, usecols=(1,)))
        except IOError:  # maximum quantum number reached set qnmax and return
            qnmax = i-1
            return bound_states, qnmax
        if len(np.atleast_1d(data)) == 1:
            data = np.reshape(data, 1)  # accounts for behavior of 0d arrays
        bound_states.append(data)
    return bound_states, qnmax


class AMDOS(object):
    """
    Atom molecule density of states class
    """
    def __init__(self, filen):
        self.systems = get_data(filen)
        self.bound_states_c = None
        self.bound_states_d = None
        self.lmax = None
        self.nmax = None
        self.lifetime = None
        self.dos = None
        self.num_open = None

    @staticmethod
    def _dos(jqn, mqn, bound_states_c, bound_states_d, vmax,
             energy_range=10.0, energy_offset=0.0):
        """
        computes the dos for an atom moleucle collsiions as described in
        Statistical aspects of ultracold resonant scattering Mayle, Ruzic and
        Bohn
        Variables
        jqn       - total jqn quantum number.
        mqn     - mj_qn projection quantum number.
        bound_states_c - boundstates of the complex obtained from _read_data
        bound_states_d - boundstates of the molecule obtained from _read_data
        vmax    - maximum vibrational quantum number.
        energy_range  - dos = number of states counted/energy range
        energy_offset - calculate dos around energy = energy_offset
        """
        if abs(mqn) > abs(jqn):
            raise Exception('physically impossible (abs(mj_qn) > abs(jqn))')
        limit = (energy_range/2.0)*K2EH  # convert energy_range K to hartrees
        energy_offset *= K2EH
        outside = np.ma.masked_outside  # mask enties in list outside range
        count = np.ma.count  # counts the number of unmasked entries in a list
        abs_gs = bound_states_d[0][0]
        num = 0
        # looping over rotational state of dimer a
        for nqn in xrange(len(bound_states_d)):
            # looping over all l constant with jqn
            for lqn in xrange(len(bound_states_c)):
                # only include pairs which couple to form jqn
                if abs(nqn-lqn) <= jqn and nqn+lqn >= jqn:
                    # looping over all vibrational levels of a
                    for vqn in xrange(min(len(bound_states_d[nqn]), vmax+1)):
                        # degeneracy
                        deg = len(xrange(max(-lqn, mqn-nqn),
                                         min(lqn, mqn+nqn)))+1
                        threshold_energy = bound_states_d[nqn][vqn]-abs_gs
                        num += count(
                            outside(bound_states_c[lqn],
                                    -limit-threshold_energy+energy_offset,
                                    limit-threshold_energy+energy_offset))*deg

        dos = (float(num)/energy_range)*1.0E-3
        lifetime = dos*1.0E3*EH2K*2.0*pi*ATOMICUNITOFTIME*1.0e9
        # return dos in mK-1 and lifetime in ns
        return dos, lifetime

    @staticmethod
    def _num_open(jqn, mqn, bound_states_c, bound_states_d,
                  vmax, energy_offset=0.0):
        """
        computes the number of threshold channels of the complex with energy
        less than energy_offset.
        Statistical aspects of ultracold resonant scattering Mayle, Rizic and
        Bohn
        VAriables
        jqn       - total J quantum number.
        mqn     - mj_qn projection quantum number.
        bound_states_c - boundstates of the complex obtained from _read_data
        bound_states_d - boundstates of the molecule obtained from _read_data
        vmax    - maximum vibrational quantum number.
        energy_offset - calculate number of threshold channels of the complex
                        with energy less than energy_offset.
        """
        if abs(mqn) > abs(jqn):
            print 'physically impossible (abs(mj_qn) > abs(jqn))'
            return
        energy_offset *= K2EH
        abs_gs = bound_states_d[0][0]
        num_open = 0
        # looping over rotational state of dimer a
        for nqn in xrange(len(bound_states_d)):
            # looping over all l constant with jqn
            for lqn in xrange(len(bound_states_c)):
                # include pairs of rotational levels which couple to form jqn
                if abs(nqn-lqn) <= jqn and nqn+lqn >= jqn:
                    # looping over all vibrational levels of a
                    for vqn in xrange(min(len(bound_states_d[nqn]), vmax+1)):
                        threshold_energy = bound_states_d[nqn][vqn]-abs_gs
                        if threshold_energy < energy_offset:
                            num_open += 1
        return num_open

    def get_dos(self, key, jqn=0, mqn=0, nmax=100, vmax=9999,
                energy_offset=0.0):
        """
        simple wrapper around _read_data and AMDOS
        Variables
        jqn       - total J quantum number.
        mqn     - mj_qn projection quantum number.
        nmax    - maximum allowed rotaional quantum number of the molecule.
        vmax    - maximum allowed vibrational quantum number of the molecule.
        energy_offset - calculate dos around energy = energy_offset
        """
        lmax = nmax + jqn
        self.bound_states_d, self.nmax = _read_data(
            self.systems[key]['dimer_dirn'],
            self.systems[key]['dimer_filen'],
            nmax)
        self.bound_states_c, self.lmax = _read_data(
            self.systems[key]['cmplx_dirn'],
            self.systems[key]['cmplx_filen'],
            lmax)
        self.dos, self.lifetime = AMDOS._dos(
            self,
            jqn,
            mqn,
            self.bound_states_c,
            self.bound_states_d,
            vmax,
            energy_offset=energy_offset)
        return self.dos, self.lifetime

    def get_num_open(self, key, jqn=0, mqn=0, nmax=100, vmax=9999,
                     energy_offset=0.0):
        """
        simple wrapper around _read_data and _num_open
        Variables
        jqn       - total jqn quantum number.
        mqn     - mj_qn projection quantum number.
        nmax    - maximum allowed rotaional quantum number of the molecule.
        vmax    - maximum allowed vibrational quantum number of the molecule.
        energy_offset - calculate dos around energy = energy_offset
        """
        lmax = nmax + jqn
        self.bound_states_d, self.nmax = _read_data(
            self.systems[key]['dimer_dirn'],
            self.systems[key]['dimer_filen'],
            nmax)
        self.bound_states_c, self.lmax = _read_data(
            self.systems[key]['cmplx_dirn'],
            self.systems[key]['cmplx_filen'],
            lmax)
        self.num_open = AMDOS._num_open(
            jqn,
            mqn,
            self.bound_states_c,
            self.bound_states_d,
            vmax,
            energy_offset)
        return self.num_open


class MMDOS(object):
    """
    class for molecule molecule density of states calculations
    """
    def __init__(self, filen):
        self.systems = get_data(filen)
        self.nmax = None
        self.lmax = None
        self.bound_states_c = None
        self.bound_states_d = None
        self.lifetime = None
        self.dos = None

    @staticmethod
    def _mm_dos(jqn, mj_qn, nmax, vmax, lmax,
                bound_states_d, bound_states_c,
                energy_range=10.0):
        """
        computes the dos for an atom moleucle collsiions as described in
        Scattering of Ultracold Molecules in the Highly Resonant Regime --
        Mayle, Ruzic, Quememer and Bohn
        Variables
        jqn       - total J quantum number.
        mqn     - mj_qn projection quantum number.
        bound_states_c - boundstates of the complex obtained from _read_data
        bound_states_d - boundstates of the molecule obtained from _read_data
        vmax    - maximum vibrational quantum number.
        lmax    - maximum end-over-end rotational quantum number of the two
                  molecules
        nmax    - maximum rotational quantum number of a single molecule.
        energy_range  - dos = number of states counted/energy range
        """
        if abs(mj_qn) > abs(jqn):
            raise Exception('physically impossible (abs(mj_qn) > abs(jqn))')
        num = 0  # variable to hold the number of states between limits
        limit = (energy_range/2.0)*K2EH  # convert energy_range K to hartrees
        abs_gs = bound_states_d[0][0]  # energy of the absolute ground state
        for nqn in xrange(max(0, jqn-lmax), min(2*nmax, lmax-jqn)+1):
            for lqn in xrange(abs(jqn-nqn), min(lmax, jqn+nqn)+1):
                for n1qn in xrange(max(0, nqn-nmax), nmax+1):
                    for n2qn in xrange(abs(nqn-n1qn), min(nqn+n1qn, nmax)+1):
                        # looping over all vibrational levels of dimer1
                        for v1qn in xrange(min(len(bound_states_d[n1qn]),
                                               vmax+1)):
                            # looping over all vibrational levels of dimer2
                            for v2qn in xrange(min(len(bound_states_d[n2qn]),
                                                   vmax+1)):
                                threshold_energy = (
                                    bound_states_d[n1qn][v1qn] - abs_gs +
                                    bound_states_d[n2qn][v2qn] - abs_gs)
                                if (bound_states_c[lqn][0] >
                                        limit-threshold_energy):
                                    # the lowest state is higher than highest
                                    # threshold
                                    break
                                else:
                                    start_end = np.searchsorted(
                                        bound_states_c[lqn],
                                        [-limit-threshold_energy,
                                         limit-threshold_energy],
                                        'left')
                                    deg = (abs(max(-lqn, mj_qn-nqn)) +
                                           min(lqn, mj_qn+nqn) + 1)
                                    num += deg*(start_end[1] - start_end[0])
        # return dos in per uK i.e.times 10**-6K
        dos = (float(num)/energy_range)*1.0E-6
        # return lifetime in ms
        lifetime = dos*1.0E6*EH2K*2.0*pi*ATOMICUNITOFTIME * 1.0e3
        return dos, lifetime

    def mm_dos(self, key, jqn=0, mqn=0, nmax=100, vmax=999):
        """
        simple wrapper around _read_data and MMDOS
        Variables
        jqn       - total J quantum number.
        mqn     - mj_qn projection quantum number.
        nmax    - maximum allowed rotaional quantum number of the molecule.
        vmax    - maximum allowed vibrational quantum number of the molecule.
        energy_offset - calculate dos around energy = energy_offset
        """
        lmax = 2*nmax+jqn
        self.bound_states_d, self.nmax = _read_data(
            self.systems[key]['dimer_dirn'],
            self.systems[key]['dimer_filen'],
            nmax)
        self.bound_states_c, self.lmax = _read_data(
            self.systems[key]['cmplx_dirn'],
            self.systems[key]['cmplx_filen'],
            lmax)
        self.dos, self.lifetime = MMDOS._mm_dos(
            self, jqn, mqn, self.nmax, vmax, self.lmax,
            self.bound_states_d, self.bound_states_c)
        return self.dos, self.lifetime
