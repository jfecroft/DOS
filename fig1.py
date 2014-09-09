import DOSModule
import numpy as np
import pylab as plt
from jfec_colour import jfec_colour_scheme


filen ='Jfec_AlkaliAtomPlusDimer.txt'
Jfec_AlkaliAtomPlusDimer = DOSModule.AtomMoleculeDOS(filen)
filen ='Mayle_AlkaliAtomPlusDimer.txt'
Mayle_AlkaliAtomPlusDimer = DOSModule.AtomMoleculeDOS(filen)
filen ='HeliumPlusHydrocarbon.txt'
HeliumPlusHydrocarbon = DOSModule.AtomMoleculeDOS(filen)
#choose which set of systems to use
systems = HeliumPlusHydrocarbon
#set up color scheme
colors = jfec_colour_scheme(len(systems.systems.keys()),0,0.8)


J,MQN = 0,0 #total J quantum number
nmax,vmax = 100, 9999#maximum rotational and vibrational quantum numbers
#lmax is determined by J and nmax
collision_energies = np.logspace(-2,3,11)
#collision_energy = 1.0e-9


for system in systems.systems.keys():
 print system
 color = colors.next()
 NumOpen_Iter = (systems.get_num_open(system,J=J,MQN=MQN,nmax=nmax,vmax=vmax,energy_offset=collision_energy) for collision_energy in collision_energies)
 NumOpen_J3_Iter = (systems.get_num_open(system,J=3,MQN=MQN,nmax=nmax,vmax=vmax,energy_offset=collision_energy) for collision_energy in collision_energies)
 Lifetime_Iter = (systems.get_dos(system,J=J,MQN=MQN,nmax=nmax,vmax=vmax)[1] for collision_energy in collision_energies)
 Lifetime_J3_Iter = (systems.get_dos(system,J=3,MQN=MQN,nmax=nmax,vmax=vmax)[1] for collision_energy in collision_energies)
 NumOpen  = np.fromiter(NumOpen_Iter,np.int_)
 NumOpen_J3  = np.fromiter(NumOpen_J3_Iter,np.int_)
 Lifetime = np.fromiter(Lifetime_Iter,np.float)
 Lifetime_J3 = np.fromiter(Lifetime_J3_Iter,np.float)
 plt.plot(collision_energies,Lifetime/NumOpen,label='Helium+'+system.split('_')[1].capitalize(),color=color,linewidth=2.0)
 plt.plot(collision_energies,Lifetime_J3/NumOpen_J3,'--',color=color,linewidth=2.0)


ErangeLine = np.logspace(-2,3,99)
plt.plot(ErangeLine,5.0/ErangeLine,'-.',linewidth=2.0,color='k')
plt.rcParams.update({'font.size': 22})
plt.legend(loc=0,shadow=True)
plt.axvline(x=80.0,linewidth=2.0,color='r')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Collision Energy (K)')
plt.ylabel('Lifetime (ns)')
plt.show()
