DOS
===

Density of states lifetime calculations. Calculates density of states (dos) and lifetime as described in the papers
-Statistical Aspects of Ultracold Resonant Scattering-- M. Mayle, B. P. Ruzic, and J. L. Bohn, Phys. Rev. A 85, 062712 (2012).
-Scattering of Ultracold Molecules in the Highly Resonant Regime -- M. Mayle, G. Quéner, B. P. Ruzic, and J. L. Bohn, Phys. Rev. A 87, 012709 (2013).
Lifetimes are computed from rovibrational energy level data contained in the folder data.
The data files used in the papers and 
-Long-Lived Complexes and Chaos in Ultracold Molecular Collisions  J. F. E. Croft and J. L. Bohn, Phys. Rev. A 89, 102714 (2014). 
are provided.

the text file .txt contain directories and file names for each system.

1d_schrodinger.F can be used to compute rovibrational energy levels for a morse potential and can simply be modified for any other 1d potential.

running the following commands will run a simple example calculation.
$python dos_molecule_molecule.py should provide the same output as the text file output_dos_molecule_molecule.txt 
$python dos_atom_molecule.py should provide the same output as the text file output_dos_atom_molecule.txt 

