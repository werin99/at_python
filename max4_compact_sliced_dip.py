# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 09:07:07 2019

@author: svewer

Adopted from a MATLAB input file 191024

Using sliced dipole magnets the fractional tunes (not integers) and chromaticity 
matches the AT MATLAB results.


%+++++++++++++++++++++
% Reduced input file by Sverker from Magnus original
% 160107
%+++++++++++++++++++++

%  Loads an Accelerator Toolbox lattice model of the MAX IV 3 GeV storage
%  ring into the cell array THERING.
%
%  Author: Magnus Sjöström, MAX IV laboratory (2014-11-11). Messed up by SW 151221.
%


"""

import numpy as np
import at
import matplotlib.pyplot as plt


# Corrector definitions
kick=np.zeros(2)
corrh         =     at.elements.Corrector('corrh',0.0, kick  )
corrv         =     at.elements.Corrector('corrv',0.0, kick  )
corr_d        =     at.elements.Corrector('corrd',0.025, kick  )
Corr          =     [corr_d, corrv, corr_d, corr_d, corrh, corr_d]


# Straight Section definitions

str0500  = at.elements.Drift('str0500', 0.500000)
strx403  = at.elements.Drift('strx403', 0.403110)
str0321  = at.elements.Drift('str0321', 0.321000)
str0302  = at.elements.Drift('str0302', 0.302000)
strx203  = at.elements.Drift('strx203', 0.203110)
str0100  = at.elements.Drift('str0100', 0.100000)
str0075  = at.elements.Drift('str0075', 0.075000)
str0010  = at.elements.Drift('str0010', 0.010000)
strx006  = at.elements.Drift('strx006', 0.006080)

str_sw1 = at.elements.Drift('str_sw1', 0.354)
str_sw2 = at.elements.Drift('str_sw2', 0.225)
str_sw3 = at.elements.Drift('str_sw3', 0.36268)

LongStr         = at.elements.Drift('LongStr', 2.32100)
ShortStr        = at.elements.Drift('ShortStr', 1.302000)

# Quadrupoles
qf      = at.elements.Quadrupole('qf'    , 0.150000, 4.030076)
qfm   	= at.elements.Quadrupole('qfm'   , 0.150000,  3.773995 )   
qfend 	= at.elements.Quadrupole('qfend' , 0.250000,3.653849) 
qdend 	= at.elements.Quadrupole('qdend' , 0.250000, -2.503663)   

# Sextupoles
sd    = at.elements.Sextupole('sd'    , 0.10000, -116.625229 )
sdend	= at.elements.Sextupole('sdend' , 0.10000,  -170.000000 )   
sfm  	= at.elements.Sextupole('sfm'   , 0.10000,  170.000000)   
sfo 	= at.elements.Sextupole('sfo'   , 0.10000,  174.000000 )  
sfi   = at.elements.Sextupole('sfi'   , 0.10000, 207.412038)

# Bending magnet definitions
dftot = at.elements.Bend('dip', 1.22378, 0.052357678, -0.6970889)
dmctot = at.elements.Bend('dmctot', 0.75424, 0.026191039, -0.5578045)

# SLICED Bending magnet definitions
pi=np.pi
d0 = at.elements.Bend('dip', 0.36189000000,1.094181*pi/180,-0.864858)
df1 = at.elements.Bend('dip', 0.050000,0.151199*pi/180,-0.864908)
df2 = at.elements.Bend('dip',0.050000, 0.151101*pi/180,-0.866059)
df3 = at.elements.Bend('dip',0.050000, 0.101861*pi/180,-0.551829)
df4 = at.elements.Bend('dip',0.050000,0.001569*pi/180,+0.011759)
df5 = at.elements.Bend('dip',0.050000,0.000089*pi/180,-0.000128)

dm5 = at.elements.Bend('dipm',0.050000,0.000090*pi/180,-0.000129)
dm4 = at.elements.Bend('dipm',0.050000,0.001579*pi/180,+0.011839)
dm3 = at.elements.Bend('dipm',0.050000,0.102549*pi/180,-0.555557)
dm2 = at.elements.Bend('dipm',0.050000,0.152122*pi/180,-0.871910)
dm1 = at.elements.Bend('dipm',0.050000, 0.152220*pi/180,-0.870751)

ds6 = at.elements.Bend('dipm',0.050000,0.001070*pi/180,+0.006608)
ds5 = at.elements.Bend('dipm',0.050000,0.050729*pi/180,-0.271428)
ds4 = at.elements.Bend('dipm',0.050000,0.074672*pi/180,-0.425119)
ds3 = at.elements.Bend('dipm',0.050000,0.076248*pi/180,-0.426048)
ds2 = at.elements.Bend('dipm',0.050000,0.114983*pi/180,-0.584884)
ds1 = at.elements.Bend('dipm',0.050000,0.152049*pi/180,-0.870351)
ds0 = at.elements.Bend('dipm',0.204240,0.621695*pi/180,-0.870701)

# Gather dipole slices
Dip  = [df5, df4, df3, df2, df1, d0];
Dipm = [ds6, ds5, ds4, ds3, ds2, ds1, ds0, dm1, dm2, dm3, dm4, dm5];

# Dipole magnet units
Dip     =   Dip + Dip[::-1]


# RING STRUCTURE 
#--------------------------------------------------------------------------

# Virtual units containing split quadrupoles with central sextupole
Sqfm =      [ qfm, str0075, sfm, str0075, qfm ]
Sqfo =      [ qf,  str0075, sfo, str0075, qf ]
Sqfi =      [ qf,  str0075, sfi, str0075, qf ]

# Construct magnet sections
Cell = [str0100, str0100, strx203, sd, str0010]+ Dip+[ str0010, sd, strx403 ]
Cell_corr = [str0100]+ Corr+[ strx203, sd, str0010]+Dip+[str0010, sd, strx403 ]

Uc1     = Sqfm + Cell
Uc1_corr     = Sqfm + Cell_corr
Uc2     = Sqfo + Cell
Uc3     = Sqfi + Cell+ Sqfi 
Uc4     = Uc2[::-1]
Uc5     = Uc1[::-1]

Mc     = [ str_sw1, qfend, str_sw2, qdend, strx006]+Dipm+[str_sw3, sdend ]

Achr  = [ LongStr]+ Mc+[ShortStr]+ Uc1 + Uc2 + Uc3+ Uc4+ Uc5+ [ShortStr]+ Mc[::-1]+[ LongStr ]

Achr_corr  = [ LongStr]+ Mc+[ ShortStr,]+Uc1_corr + Uc2 + Uc3 + Uc4+ Uc5 +[ShortStr]+ Mc[::-1] +[LongStr ]

Ring    =  Achr_corr + 19*Achr
#Ring= Achr
#Ring = [ Achr_corr ];

#Build lattice
ELIST = Ring;
# ELIST = achr_corr;

#uildlat(ELIST);
THERING=at.lattice.Lattice(ELIST, energy=3,periodicity=1)
length=np.size(Ring)
refpts=np.r_[0:length+1]

optics=at.physics.linopt(ELIST, refpts=refpts, get_chrom=True)
beta=optics[3]['beta']
disp=optics[3]['dispersion']

_, tunes, chroms, twiss = THERING.linopt(0.0, 0, get_chrom=True, coupled=False)

print('Tunes:',tunes)
print('Beta: \n',twiss['beta'])

print('Circumference =',THERING.circumference)
print('Periodicity =',THERING.periodicity)

print('Chromaticity =',chroms)

s=at.lattice.get_s_pos(ELIST)

plt.plot(s,beta[:,0],color='green')
plt.plot(s,beta[:,1])
axes=plt.gca()
axes.set_xlabel('s (m)')
axes.set_ylabel('beta (m)')
axes.set_title('Beta functions')
plt.show()




plt.plot(s,disp[:,0],'r')
axes=plt.gca()
axes.set_xlabel('s (m)')
axes.set_ylabel('dispersion (m)')
axes.set_title('Dispersion function')
plt.show()

