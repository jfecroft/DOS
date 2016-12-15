"""
Reproduce published results
"""
from DOSModule import get_data, get_dos
from tabulate import tabulate

FILEN = 'Jfec_AlkaliAtomPlusDimer.yml'
JFEC_ALKALIATOMPLUSDIMER = get_data(FILEN)
FILEN = 'Mayle_AlkaliAtomPlusDimer.yml'
MAYLE_ALKALIATOMPLUSDIMER = get_data(FILEN)
FILEN = 'HeliumPlusHydrocarbon.yml'
HELIUMPLUSHYDROCARBON = get_data(FILEN)


# JQN, MQN = 0, 0  # total J quantum number
# NMAX, VMAX = 100, 9999  # maximum rotational and vibrational quantum numbers
# lmax is determined by JQN and NMAX
# COLLISION_ENERGY = 0.0

header = ['system', 'dos(mK-1)','lt(ns)']
output = []
# printing dos and lt for alkali atom + dimer collisions as published in
# Long-Lived Complexes and Chaos in Ultracold Molecular Collisions
#            J. F. E. Croft and J. L. Bohn, Phys. Rev. A 89, 102714 (2014).
for system, input_dict in JFEC_ALKALIATOMPLUSDIMER.iteritems():
    dos, lt = get_dos(**input_dict)
    output.append([system, dos, lt])


# printing dos and lt for alkali atom + dimer collisions as published in
# statistical Aspects of Ultracold Resonant Scattering --
# M. Mayle, B. P. Ruzic, and J. L. Bohn, Phys. Rev. A 85, 062712 (2012).
for system, input_dict in MAYLE_ALKALIATOMPLUSDIMER.iteritems():
    dos, lt = get_dos(**input_dict)
    output.append([system, dos, lt])

# printing dos and lt for helium + hydrocarbon collisions
# as yet unpublished
for system, input_dict in HELIUMPLUSHYDROCARBON.iteritems():
    dos, lt = get_dos(**input_dict)
    output.append([system, dos, lt])

print tabulate(output, headers=header)
