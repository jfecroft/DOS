import DOSModule

filen ='Jfec_AlkaliAtomPlusDimer.txt'
Jfec_AlkaliAtomPlusDimer = DOSModule.AtomMoleculeDOS(filen)
filen ='Mayle_AlkaliAtomPlusDimer.txt'
Mayle_AlkaliAtomPlusDimer = DOSModule.AtomMoleculeDOS(filen)
filen ='HeliumPlusHydrocarbon.txt'
HeliumPlusHydrocarbon = DOSModule.AtomMoleculeDOS(filen)


J,MQN = 0,0 #total J quantum number
nmax,vmax = 100, 9999#maximum rotational and vibrational quantum numbers
#lmax is determined by J and nmax
collision_energy = 0.0

print 'system  dos(mK-1) lt(ns)'
#printing dos and lt for alkali atom + dimer collisions as published in
#Long-Lived Complexes and Chaos in Ultracold Molecular Collisions
#            J. F. E. Croft and J. L. Bohn, Phys. Rev. A 89, 102714 (2014). 
for system in Jfec_AlkaliAtomPlusDimer.systems.keys():
 dos, lt =Jfec_AlkaliAtomPlusDimer.get_dos(system,J=J,MQN=MQN,nmax=nmax,vmax=vmax,energy_offset=collision_energy)
 print system,dos,lt


#printing dos and lt for alkali atom + dimer collisions as published in
#statistical Aspects of Ultracold Resonant Scattering-- M. Mayle, B. P. Ruzic, and J. L. Bohn, Phys. Rev. A 85, 062712 (2012).
for system in Mayle_AlkaliAtomPlusDimer.systems.keys():
 dos, lt = Mayle_AlkaliAtomPlusDimer.get_dos(system,J=J,MQN=MQN,nmax=nmax,vmax=vmax,energy_offset=collision_energy)
 print system,dos,lt

#printing dos and lt for helium + hydrocarbon collisions 
#as yet unpublished
for system in HeliumPlusHydrocarbon.systems.keys():
 dos, lt = HeliumPlusHydrocarbon.get_dos(system,J=J,MQN=MQN,nmax=nmax,vmax=vmax,energy_offset=collision_energy)
 print system,dos,lt

