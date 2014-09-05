import DOSModule
import numpy as np

filen ='Jfec_AlkaliAtomPlusDimer.txt'
Jfec_AlkaliAtomPlusDimer = DOSModule.AtomMoleculeDOS(filen)
filen ='Mayle_AlkaliAtomPlusDimer.txt'
Mayle_AlkaliAtomPlusDimer = DOSModule.AtomMoleculeDOS(filen)
filen ='HeliumPlusHydrocarbon.txt'
HeliumPlusHydrocarbon = DOSModule.AtomMoleculeDOS(filen)


J,MQN = 0,0 #total J quantum number
nmax,vmax = 100, 9999#maximum rotational and vibrational quantum numbers
#lmax is determined by J and nmax
collision_energies = np.logspace(-2,1,11)
collision_energy = 1.0e-9

print 'system  dos(mK-1) lt(ns)'
#printing dos and lt for alkali atom + dimer collisions as published in
#Long-Lived Complexes and Chaos in Ultracold Molecular Collisions
#            J. F. E. Croft and J. L. Bohn, Phys. Rev. A 89, 102714 (2014). 
for system in Jfec_AlkaliAtomPlusDimer.systems.keys():
 Jfec_AlkaliAtomPlusDimer.get_dos(system,J=J,MQN=MQN,nmax=nmax,vmax=vmax,energy_offset=collision_energy)
 NumOpen_Iter = (Jfec_AlkaliAtomPlusDimer.get_num_open(system,J=J,MQN=MQN,nmax=nmax,vmax=vmax,energy_offset=collision_energy) for collision_energy in collision_energies)
 for i in range(8):
  print NumOpen_Iter.next()
# NumOpen = np.fromiter(NumOpen_Iter,np.int_)
# print NumOpen
# print system,Jfec_AlkaliAtomPlusDimer.dos,Jfec_AlkaliAtomPlusDimer.lt, Jfec_AlkaliAtomPlusDimer.num_open
