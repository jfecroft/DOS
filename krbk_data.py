# port scipy
from scipy.constants import physical_constants, electron_mass
amu2au = physical_constants["atomic mass unit-kilogram relationship"][0]/electron_mass
wavenumbers2hartree = physical_constants["inverse meter-hartree relationship"][0]*100

# data from
# K2: http://link.springer.com/article/10.1140/epjd/e2007-00307-2
# KRb : http://journals.aps.org/pra/abstract/10.1103/PhysRevA.76.022511

C6_KRb = 4299.51# a.u.
C6_K2 = 3920.99# a.u.

De_KRb = 4217.81510# cm-1
De_K2 = 4450.90650# cm-1

massK40 = 39.963999*amu2au
mass87Rb = 86.909187*amu2au

data = {'KRb': {'De': De_KRb, 'C6': C6_KRb, 'mass1': mass87Rb, 'mass2': massK40},
        'KK': {'De': De_K2, 'C6': C6_K2, 'mass1': massK40, 'mass2': massK40},
        'KK_Rb': {'De': 2.0*De_KRb, 'C6': 2.0*C6_KRb, 'mass1': 2*massK40, 'mass2': mass87Rb},
        'KRb_K': {'De': De_K2+De_KRb, 'C6': C6_K2+C6_KRb, 'mass1':mass87Rb+massK40, 'mass2': massK40}}


for sys, sys_data in data.iteritems():
        sys_data['De'] = sys_data['De']*wavenumbers2hartree
        sys_data['C12'] = sys_data['C6']**2/(4*sys_data['De'])
        sys_data['sigma'] = (sys_data['C12']/sys_data['C6'])**(1.0/6.0)
        print sys, sys_data
