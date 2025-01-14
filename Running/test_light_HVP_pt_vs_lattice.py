# compute the light HVP PI(Q^2) - PI(Q^2/4) in perturbation theory 
# to be  compared with non-perturbative determination
from adlerpy.adler_routines import Particle
from adlerpy.adler_sm import adler_charm_pert # The charm quark connected contribution to the Adler function
from adlerpy.adler_sm import adler_bottom_pert # The bottom quark connected contribution to the Adler function
from adlerpy.adler_sm import adler_light_pert # The light quark connected contribution to the Adler function
from adlerpy.adler_sm import adler_OZI_pert # The quark disconnected contribution to the Adler function
from adlerpy.adler_sm import adler_massless_connected
from adlerpy.adler_sm import alphas
from matplotlib import pyplot as plt  # To make some plots. 
import pandas as pd  # to create tables. 
from scipy.integrate import quad # to do the integrals. 
import numpy as np
import time

import bootstrap
import statistics

## create configs 
SAMPLES = 1000
NBIN = 100

MZ_config = bootstrap.bootstrap(91.1876, 0.0021, SAMPLES, NBIN ) # Z boson mass PDG
aZ_config = bootstrap.bootstrap(0.1183, 0.0007, SAMPLES, NBIN ) # alpha_z at Z pole from FLAG24
# aZ_config = bootstrap.bootstrap(0.1185, 0.0016, SAMPLES, NBIN ) # alpha_z at Z pole from adlerpy paper
# aZ_config = bootstrap.bootstrap(0.1185, 0.0030, SAMPLES, NBIN ) # alpha_z at Z pole from Mainz22

# MSbar quark masses
mUp_config = bootstrap.bootstrap(2.16*0.001, 0.07*0.001, SAMPLES, NBIN)     # PDG
mDown_config = bootstrap.bootstrap(4.7*0.001, 0.07*0.001, SAMPLES, NBIN)    # PDG
mStrange_config = bootstrap.bootstrap(93.5*0.001, 0.8*0.001, SAMPLES, NBIN) # PDG

mC_config = bootstrap.bootstrap(1.278, 0.006, SAMPLES, NBIN) # scale invariant mC (FLAG)
# mC_config = bootstrap.bootstrap(1.27, 0.02, SAMPLES, NBIN) # scale invariant mC (PDG)
mB_config = bootstrap.bootstrap(4.171, 0.02, SAMPLES, NBIN) # scale invariant mB (FLAG)

# heavy quark pole mass
mC_Pole_config = bootstrap.bootstrap(1.67, 0.07, SAMPLES, NBIN) # scale invariant mC (PDG)
mB_Pole_config = bootstrap.bootstrap(4.78, 0.06, SAMPLES, NBIN) # scale invariant mC (PDG)

# Condensates
GG_config = bootstrap.bootstrap(0.012, 0.012, SAMPLES, NBIN)  # gluon condensate https://arxiv.org/pdf/2302.01359
QQ_config = bootstrap.bootstrap(0.0013, 0.007, SAMPLES, NBIN) # quark condensate https://arxiv.org/pdf/2302.01359

alpha0=1/137.035999180
QVAL = np.sqrt(7) # GeV
mmu=0.105658 # muon mas

isov_contrib = np.empty(NBIN)

def integrand_light(Q2,aZ,Mz,particles, GG, qq, nloops, QED):
    Q = np.sqrt(Q2)
    return 1 / (12*np.pi**2) * (adler_light_pert(aZ, Mz, Q=Q, particles=particles, GG=GG,qq=qq, nloops=nloops, QED=QED))/Q2;

def Pi_light(Q2, aZ, Mz, particles, GG, qq, nloops, QED):
    return quad(integrand_light, Q2/4, Q2, args=(aZ,Mz,particles,GG,qq,nloops,QED))[0]


for k in range(0, NBIN):
    Mz = MZ_config[k]
    aZ = aZ_config[k]
    mup = mUp_config[k]
    mdw = mDown_config[k]
    mst = mStrange_config[k]
    mc0 = mC_config[k]
    mb0 = mB_config[k]
    mc0_pole = mC_Pole_config[k]
    mb0_pole = mB_Pole_config[k]
    #Here you define the SM particles 
    up=Particle("up",x=[mup,2,3],mudec=0.001,mpole=None,mpole_on=False)
    down=Particle("down",x=[mdw,2,3],mudec=0.001,mpole=None,mpole_on=False)
    strange=Particle("strange",x=[mst,2,3],mudec=0.001,mpole=None,mpole_on=False)
    charm=Particle("charm",x=[mc0,mc0,4],mudec=2*(mc0),mpole=mc0_pole,mpole_on=False)
    bottom=Particle("bottom",x=[mb0,mb0,5],mudec=2*mb0,mpole=mb0_pole,mpole_on=False)
    top=Particle("top",x=[164,164,6],mudec=164,mpole=164,mpole_on=False) # not necessary in this case.
    particle_list=[up,charm,down,strange,bottom,top]

    GG=GG_config[k];
    qq=QQ_config[k];
    nloops=5
    QED=True


    isov_contrib[k] =  3 / 2 * 0.5 * Pi_light(QVAL**2, aZ, Mz, particle_list, GG, qq, nloops, QED)
    # isov_contrib[k] *= 1e5

mean_tot = np.mean(isov_contrib) 
std_tot = np.std(isov_contrib) 
print("     -At ", nloops, " loops: ", mean_tot, std_tot, "\n")
