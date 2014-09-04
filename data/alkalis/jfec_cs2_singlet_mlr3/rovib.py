import subprocess 
import numpy as np
import os
import numpy.ma as ma
import re
from   tempfile import mkstemp, mkdtemp
import shutil
import filecmp
import scipy.constants

#########################################
#replace will search through a file for a specific word an then replace that line
def replace(file, pattern, subst):
    p = re.compile(pattern)
    #Create temp file
    fh, abs_path = mkstemp()
    new_file = open(abs_path,'w')
    old_file = open(file)
    for line in old_file:
        if p.match(line): #using match because for input files only want to replace the currently used variable
         line = pattern + ' = ' + str(subst) + ', \n'
        new_file.write(line)
    #close temp file
    new_file.close()
    os.close(fh)
    old_file.close()
    os.remove(file)
    shutil.move(abs_path, file)

#routine which call 1d_schrodinger eqn solver and returns 
#all the levels below zero outputted
def run_1d_schrodinger(inputfile_name,outputfile_name,L):
 home = os.getcwd()
 replace(inputfile_name, ' L', L) #editing inputfile such that L=L is called
 subprocess.call(home+"/1d_schrodinger.x < " + inputfile_name, stdout=open(os.devnull, 'w'), shell=True)
 return()

############################################
lmin = 499
lmax = 700
inputfile         = 'input_CsCs.txt'
outputfile        = 'fort.10'
sys               = 'cscs'

C6 = 6870.100 #in au  from JCP 132 094105 2010 Coxon
De = 3649.847 # in cm-1  from JCP 132 094105 2010 Coxon
#converting cm-1 to Hartree
De = De*100.0/(scipy.constants.Rydberg*2.0)
C12 = C6**2/(4.0*De)

replace(inputfile, ' C6', C6) #editing inputfile such that L=L is called
replace(inputfile, ' C12', C12) #editing inputfile such that L=L is called
#generate to states of the dimer for different n upto nmax
for i in range(lmin,lmax+1):
 run_1d_schrodinger(inputfile,outputfile,i)
 if filecmp.cmp(outputfile,sys+'_results_j'+str(i-1)+'.dat') != True:
  shutil.copyfile(outputfile,sys+'_results_j'+str(i)+'.dat')
 else:
  break

