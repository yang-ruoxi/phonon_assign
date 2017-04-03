# M mode

* M2+ mode, M_3+ mode are both in-phase tilting modes, the difference is the origin shift. (M2+: (1/2,1/2,0); M3+: (0,0,0))
  They both lead the parent cubic structure to P4/mbm phase, space group number 127. 
  
* M2- mode is collective motions of B cation and X anion, two elements can have different amplitude so that creates a Jah-Teller type of distortion: five B-X bonding lengthens and one shortens. 
  This leads the cubic phase to 129 P4/nmm. 
  
* M3- mode is A cation moving out-of phase. 
  This also leads to 129 P4/nmm space group.
  
* M5-: This is collective motions of all atoms, including A cation anti-phase displacement and octahedral distortion. 
  This leads to 10 P2/m space group. 
  

The above modes are found in ABX3 compounds. 

First, a modulation structure needs to be created: `phonopy --dim='x x x' modulation.conf`, e.g. where 

`MODULATION = 2 2 2, 0.5 0.5 0 1 20`

shoud be included in the .conf file. For details check the phonopy website.  

Once you get the MPOSCAR, rename it to POSCAR, and check its symmetry by doing `phonopy --symmetry` (because phonopy takes POSCAR by default).

Now you get the symmetry, go to Stokes' ISODISTORT (http://stokes.byu.edu/iso/isotropy.php), check what modes lead to such space group, and you can "view distortion" on Stokes' website and compare the distortion to your MPOSCAR distortion. 

By doing so you can decide which mode is responsible for the certain branch at certain q-points. 

