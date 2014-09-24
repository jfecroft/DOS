import numpy as np
import scipy.constants
from math import pi
from collections import OrderedDict


#defines some physical constants
EhtoK = scipy.constants.physical_constants["hartree-kelvin relationship"][0]
KtoEh = scipy.constants.physical_constants["kelvin-hartree relationship"][0]
AtomicUnitOfTime = scipy.constants.physical_constants["atomic unit of time"][0]

def get_data(filen):
 """reads a file of the form 1st line heading columns
    2nd lines on values for those headings
    returns a nested dictionary 
    converts values to float if possible"""
 fo = open(filen, "r")
 lines = fo.readlines()
 keys = lines[0].split()

 data_dict = OrderedDict()

 for line in lines[1:]:
  items = line.split()
  for i, item in enumerate(items):
   try: items[i] = float(item)
   except: pass
  sub_dict =  dict(zip(keys,items))
  data_dict.update({sub_dict['sys']:sub_dict})
 return data_dict


def _read_data(dir,inputfile,qnmax):
 """read in the data from precalculated rovibrational energy levels of the system
 files are og the form 'txtjXXXtxt' where XXX is the L quantum number
 starting at 0
 """
 bs = [] #list to store bound states for different L quntum numbers
 for i in xrange(qnmax+1):
  file = dir + inputfile.replace('j','j'+str(i))
  try: 
   data = np.sort(np.loadtxt(file,usecols=(1,)))
  except IOError: #maximum quantum number reached set qnmax and return
   qnmax = i-1
   return bs, qnmax
  if len(np.atleast_1d(data)) == 1:
   data = np.reshape(data,1) #accounts for behavior of 0d arrays
  bs.append(data)
 return bs, qnmax




class AtomMoleculeDOS:
 def __init__(self, filen):
  self.systems = get_data(filen)

 def _compute_dos(self,J,MQN,bs_c,bs_d,vmax,energy_range=10.0,energy_offset=0.0):
 """
 computes the dos for an atom moleucle collsiions as described in
 Statistical aspects of ultracold resonant scattering Mayle, Rizic and Bohn
 VAriables
 J       - total J quantum number.
 MQN     - MJ projection quantum number. 
 bs_c    - list of boundstates of the complex obtained from _read_data
 bs_d    - list of boundstates of the molecule obtained from _read_data
 vmax    - maximum vibrational quantum number.
 energy_range  - dos = number of states counted/energy range
 energy_offset - calculate dos around energy = energy_offset
 """
  if abs(MQN) > abs(J):
   print 'physically impossible (abs(MJ) > abs(J))'
   return
  limit = (energy_range/2.0)*KtoEh                    #convert number between +-energy_range K to hartrees
  energy_offset *= KtoEh
  outside = np.ma.masked_outside                      #masked all enties in in list outside range
  count   = np.ma.count                               #counts the number of unmasked entries in a list
  abs_gs = bs_d[0][0]
  num= 0
  for n in xrange(len(bs_d)):                          #looping over rotational state of dimer a
   for l in xrange(len(bs_c)):                         #looping over all l constant with J
    if abs(n-l) <= J and n+l >= J:                     #only include pairs of rotational levels which can couple to form J
     for v in xrange(min(len(bs_d[n]),vmax+1)):        #looping over all vibrational levels of a
      d = len(xrange(max(-l,MQN-n),min(l,MQN+n)))+1    #degeneracy
      threshold_energy = bs_d[n][v]-abs_gs
      num += count(outside(bs_c[l],-limit-threshold_energy+energy_offset,limit-threshold_energy+energy_offset))*d
 
  dos = (float(num)/energy_range)*1.0E-3
  lt = dos*1.0E3*EhtoK*2.0*pi*AtomicUnitOfTime*1.0e9
 #return dos in mK-1 and lifetime in ns
  return dos,lt

 def _compute_num_open(self,J,MQN,bs_c,bs_d,vmax,energy_range=10.0,energy_offset=0.0):
  if abs(MQN) > abs(J):
   print 'physically impossible (abs(MJ) > abs(J))'
   return
  limit = (energy_range/2.0)*KtoEh                    #convert number between +-energy_range K to hartrees
  energy_offset *= KtoEh
  outside = np.ma.masked_outside                      #masked all enties in in list outside range
  count   = np.ma.count                               #counts the number of unmasked entries in a list
  abs_gs  = bs_d[0][0]
  num_open= 0
  for n in xrange(len(bs_d)):                          #looping over rotational state of dimer a
   for l in xrange(len(bs_c)):                         #looping over all l constant with J
    if abs(n-l) <= J and n+l >= J:                     #only include pairs of rotational levels which can couple to form J
     for v in xrange(min(len(bs_d[n]),vmax+1)):        #looping over all vibrational levels of a
      threshold_energy = bs_d[n][v]-abs_gs
      if threshold_energy < energy_offset: num_open += 1
  return num_open

 def get_dos(self, key,J=0,MQN=0,nmax=100,vmax=9999,energy_offset=0.0):
  lmax  = nmax + J
  self.bs_d, self.nmax =  _read_data(self.systems[key]['dimer_dirn'],self.systems[key]['dimer_filen'],nmax)
  self.bs_c, self.lmax =  _read_data(self.systems[key]['cmplx_dirn'],self.systems[key]['cmplx_filen'],lmax)
  self.dos,  self.lt   =  AtomMoleculeDOS._compute_dos(self,J,MQN,self.bs_c,self.bs_d,vmax=vmax,energy_offset=energy_offset)
  return self.dos, self.lt

 def get_num_open(self, key,J=0,MQN=0,nmax=100,vmax=9999,energy_offset=0.0):
  lmax  = nmax + J
  self.bs_d, self.nmax =  _read_data(self.systems[key]['dimer_dirn'],self.systems[key]['dimer_filen'],nmax)
  self.bs_c, self.lmax =  _read_data(self.systems[key]['cmplx_dirn'],self.systems[key]['cmplx_filen'],lmax)
  self.num_open   =  AtomMoleculeDOS._compute_num_open(self,J,MQN,self.bs_c,self.bs_d,vmax=vmax,energy_offset=energy_offset)
  return self.num_open



class MoleculeMoleculeDOS:
 def __init__(self, filen):
  self.systems = get_data(filen)

 def _compute_dos_MoleculeMolecule(self,J,MJ,nmax,vmax,lmax,bs_d,bs_c,energy_range=10.0):
  """computes the density of states for molecule molecule collisions
     by adding up energy levels in energy_range of threshold"""
  if abs(MJ) > abs(J):
   print 'physically impossible (abs(MJ) > abs(J))'
   return
  num = 0                                             #variable to hold the number of states between limits
  limit = (energy_range/2.0)*KtoEh                    #convert number between +- energy_range K to hartrees
  abs_gs = bs_d[0][0]                                 #the energy of the absolubte ground state
  for N in xrange(max(0,J-lmax),min(2*nmax,lmax-J)+1):
   for l in xrange(abs(J-N),min(lmax,J+N)+1):
    for n1 in xrange(max(0,N-nmax),nmax+1):
     for n2 in xrange(abs(N-n1),min(N+n1,nmax)+1):
      for v1 in xrange(min(len(bs_d[n1]),vmax+1)):         #looping over all vibrational levels of dimer1
       for v2 in xrange(min(len(bs_d[n2]),vmax+1)):        #looping over all vibrational levels of dimer2
        threshold_energy = bs_d[n1][v1]-abs_gs + bs_d[n2][v2]-abs_gs
        if bs_c[l][0] > limit-threshold_energy:
         #the lowest state is higher than highest threshold
         #print 'breaking out of loop', v2
         break #other vibrational levels will only be higher
        else:
         start_end = np.searchsorted(bs_c[l],[-limit-threshold_energy,limit-threshold_energy],'left')
         d = abs(max(-l,MJ-N)) + min(l,MJ+N) + 1
         num += d*(start_end[1] - start_end[0])
 #return dos in per uK i.e.times 10**-6K
  dos = (float(num)/energy_range)*1.0E-6
 #return lifetime in ms 
  lt = dos*1.0E6*EhtoK*2.0*pi*AtomicUnitOfTime * 1.0e3
  return dos, lt


 def get_dos_MoleculeMolecule(self, key,J=0,MQN=0,nmax=100,vmax=999):
  lmax       =  2*nmax+J
  self.bs_d, self.nmax =  _read_data(self.systems[key]['dimer_dirn'],self.systems[key]['dimer_filen'],nmax)
  self.bs_c, self.lmax =  _read_data(self.systems[key]['cmplx_dirn'],self.systems[key]['cmplx_filen'],lmax)
  self.dos,  self.lt   =  MoleculeMoleculeDOS._compute_dos_MoleculeMolecule(self,J,MQN,self.nmax,vmax,self.lmax,self.bs_d,self.bs_c)
  return self.dos,  self.lt



