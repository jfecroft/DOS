"""
Reproduce publish molecule moleucle lifetimes
"""
import DOSModule

# these can take awhile (days) to run for full converged results

FILEN = 'Jfec_AlkaliDimerPlusDimer.yml'
JFEC_ALKALIDIMERPLUSDIMER = DOSModule.MMDOS(FILEN)
FILEN = 'Mayle_AlkaliDimerPlusDimer.yml'
MAYLE_ALKALIDIMERPLUSDIMER = DOSModule.MMDOS(FILEN)
J, MQN = 0, 0  # total J quantum number
NMAX, VMAX = 5, 5  # maximum rotational and vibrational quantum numbers
# lmax is determined by J and NMAX

print 'system  dos(uK-1) lt(ms)'

# reproducing data published in
# Scattering of Ultracold Molecules in the Highly Resonant Regime --
# M. Mayle, G. Quemener, B. P. Ruzic, and J. L. Bohn, Phys. Rev. A 87, 012709
# (2013).
for system in MAYLE_ALKALIDIMERPLUSDIMER.systems.keys():
    dos, lt = MAYLE_ALKALIDIMERPLUSDIMER.mm_dos(system,
                                                jqn=J,
                                                mqn=MQN,
                                                nmax=NMAX,
                                                vmax=VMAX)
    print system, dos, lt

for system in JFEC_ALKALIDIMERPLUSDIMER.systems.keys():
    dos, lt = JFEC_ALKALIDIMERPLUSDIMER.mm_dos(system,
                                               jqn=J,
                                               mqn=MQN,
                                               nmax=NMAX,
                                               vmax=VMAX)
    print system, dos, lt
