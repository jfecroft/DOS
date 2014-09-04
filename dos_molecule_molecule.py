import DOSModule
import sys
#these can take awhile (days) to run for full converged results

filen ='Jfec_AlkaliDimerPlusDimer.txt'
Jfec_AlkaliDimerPlusDimer = DOSModule.MoleculeMoleculeDOS(filen)
filen ='Mayle_AlkaliDimerPlusDimer.txt'
Mayle_AlkaliDimerPlusDimer = DOSModule.MoleculeMoleculeDOS(filen)
J,MQN = 0,0 #total J quantum number
nmax,vmax = 5, 5#maximum rotational and vibrational quantum numbers
#lmax is determined by J and nmax



print 'system  dos(uK-1) lt(ms)'

#reproducing data published in 
#Scattering of Ultracold Molecules in the Highly Resonant Regime -- M. Mayle, G. Quemener, B. P. Ruzic, and J. L. Bohn, Phys. Rev. A 87, 012709 (2013).
for system in Mayle_AlkaliDimerPlusDimer.systems.keys():
 Mayle_AlkaliDimerPlusDimer.get_dos_MoleculeMolecule(system,J=J,MQN=MQN,nmax=nmax,vmax=vmax)
 print system, Mayle_AlkaliDimerPlusDimer.dos, Mayle_AlkaliDimerPlusDimer.lt
 sys.stdout.flush()

for system in Jfec_AlkaliDimerPlusDimer.systems.keys():
 Jfec_AlkaliDimerPlusDimer.get_dos_MoleculeMolecule(system,J=J,MQN=MQN,nmax=nmax,vmax=vmax)
 print system, Jfec_AlkaliDimerPlusDimer.dos, Jfec_AlkaliDimerPlusDimer.lt

