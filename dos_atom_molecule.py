"""
Reproduce published results
"""
import DOSModule

FILEN = 'Jfec_AlkaliAtomPlusDimer.txt'
FILEN = 'Jfec_AlkaliAtomPlusDimer.yml'
JFEC_ALKALIATOMPLUSDIMER = DOSModule.AMDOS(FILEN)
FILEN = 'Mayle_AlkaliAtomPlusDimer.txt'
MAYLE_ALKALIATOMPLUSDIMER = DOSModule.AMDOS(FILEN)
FILEN = 'HeliumPlusHydrocarbon.txt'
HELIUMPLUSHYDROCARBON = DOSModule.AMDOS(FILEN)


JQN, MQN = 0, 0  # total J quantum number
NMAX, VMAX = 100, 9999  # maximum rotational and vibrational quantum numbers
# lmax is determined by JQN and NMAX
COLLISION_ENERGY = 0.0

print 'system  dos(mK-1) lt(ns)'
# printing dos and lt for alkali atom + dimer collisions as published in
# Long-Lived Complexes and Chaos in Ultracold Molecular Collisions
#            J. F. E. Croft and J. L. Bohn, Phys. Rev. A 89, 102714 (2014).
for system in JFEC_ALKALIATOMPLUSDIMER.systems.keys():
    dos, lt = JFEC_ALKALIATOMPLUSDIMER.get_dos(system,
                                               jqn=JQN,
                                               mqn=MQN,
                                               nmax=NMAX,
                                               vmax=VMAX,
                                               energy_offset=COLLISION_ENERGY)
    print system, dos, lt


# printing dos and lt for alkali atom + dimer collisions as published in
# statistical Aspects of Ultracold Resonant Scattering --
# M. Mayle, B. P. Ruzic, and J. L. Bohn, Phys. Rev. A 85, 062712 (2012).
for system in MAYLE_ALKALIATOMPLUSDIMER.systems.keys():
    dos, lt = MAYLE_ALKALIATOMPLUSDIMER.get_dos(system,
                                                jqn=JQN,
                                                mqn=MQN,
                                                nmax=NMAX,
                                                vmax=VMAX,
                                                energy_offset=COLLISION_ENERGY)
    print system, dos, lt

# printing dos and lt for helium + hydrocarbon collisions
# as yet unpublished
for system in HELIUMPLUSHYDROCARBON.systems.keys():
    dos, lt = HELIUMPLUSHYDROCARBON.get_dos(system,
                                            jqn=JQN,
                                            mqn=MQN,
                                            nmax=NMAX,
                                            vmax=VMAX,
                                            energy_offset=COLLISION_ENERGY)
    print system, dos, lt
