"""
Reproduce publish molecule moleucle lifetimes
"""
from DOSModule import get_data, get_mm_dos
from tabulate import tabulate

# these can take awhile (days) to run for full converged results

FILEN = 'Jfec_AlkaliDimerPlusDimer.yml'
JFEC_ALKALIDIMERPLUSDIMER = get_data(FILEN)
FILEN = 'Mayle_AlkaliDimerPlusDimer.yml'
MAYLE_ALKALIDIMERPLUSDIMER = get_data(FILEN)
# J, MQN = 0, 0  # total J quantum number
# NMAX, VMAX = 5, 5  # maximum rotational and vibrational quantum numbers
# lmax is determined by J and NMAX


header = ['system', 'dos(mK-1)','lt(ns)']
output = []
# reproducing data published in
# Scattering of Ultracold Molecules in the Highly Resonant Regime --
# M. Mayle, G. Quemener, B. P. Ruzic, and J. L. Bohn, Phys. Rev. A 87, 012709
# (2013).
for system, input_dict in MAYLE_ALKALIDIMERPLUSDIMER.iteritems():
    dos, lt = get_mm_dos(**input_dict)
    output.append([system, dos, lt])

for system, input_dict in JFEC_ALKALIDIMERPLUSDIMER.iteritems():
    dos, lt = get_mm_dos(**input_dict)
    output.append([system, dos, lt])


output.sort()
print tabulate(output, headers=header)
