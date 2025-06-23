# compute the light HVP PI(Q^2) - PI(Q^2/4) in perturbation theory 
# to be  compared with Harvey's mathematica determination
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

Mz = 91.1876 # Z boson mass PDG
aZ_config = bootstrap.bootstrap(0.1183, 0.0007, SAMPLES, NBIN ) # alpha_z at Z pole from FLAG24

# MSbar quark masses
mup = 2.16*0.001   # [GeV] PDG
mdw = 4.7*0.001    # [GeV] PDG
mst = 93.5*0.001   # [GeV] PDG

mC_config = bootstrap.bootstrap(1.278, 0.006, SAMPLES, NBIN) # scale invariant mC (FLAG)
mB_config = bootstrap.bootstrap(4.171, 0.02, SAMPLES, NBIN) # scale invariant mB (FLAG)

# heavy quark pole mass
mC_Pole_config = bootstrap.bootstrap(1.67, 0.07, SAMPLES, NBIN) # scale invariant mC (PDG)
mB_Pole_config = bootstrap.bootstrap(4.78, 0.06, SAMPLES, NBIN) # scale invariant mC (PDG)

# Condensates
GG_config = bootstrap.bootstrap(0.012, 0.012, SAMPLES, NBIN)  # gluon condensate https://arxiv.org/pdf/2302.01359
QQ_config = bootstrap.bootstrap(-0.0013, 0.007, SAMPLES, NBIN) # quark condensate https://arxiv.org/pdf/2302.01359

alpha0=1/137.035999180
mmu=0.105658 # muon mas
QVAL = np.sqrt(10) # GeV

isov_contrib = np.empty(NBIN)

def integrand_light(Q2,aZ,Mz,particles, GG, qq, nloops, QED):
    Q = np.sqrt(Q2)
    return 1 / (12*np.pi**2) * (adler_light_pert(aZ=aZ, Mz=Mz, Q=Q, particles=particles, GG=GG,qq=qq, nloops=nloops, QED=QED))/Q2;

def Pi_light(Q2, aZ, Mz, particles, GG, qq, nloops, QED):
    return quad(integrand_light, Q2/4, Q2, args=(aZ,Mz,particles,GG,qq,nloops,QED))[0]


def adler_light(aZ, Mz, Q, particles, nloops, GG, qq, QED):
    adlight=adler_light_pert(aZ=aZ,Mz=Mz,Q=Q,particles=particles,nloops=nloops,GG=GG,qq=qq,QED=QED)
    return adlight

emin = 1.0
emax = 100.0
q_list_log = np.linspace(np.log(emin), np.log(emax), num=400) # log(Q) in GeV
q_list = np.exp(q_list_log) 

adler_tot = np.empty((len(q_list), NBIN))

for k in range(0, NBIN):
    aZ = aZ_config[k]
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

    GG= GG_config[k];
    qq= QQ_config[k];
    nloops=4
    QED=False

    for q in range(0,len(q_list)):
        ad_q = adler_light(aZ, Mz, q_list[q], particles=particle_list, nloops=nloops, GG=GG, qq=qq, QED=QED)
        adler_tot[q,k] = ad_q

    # isov_contrib[k] =  3 / 2 * 0.5 * Pi_light(QVAL**2, aZ, Mz, particle_list, GG, qq, nloops, QED) # with charge factor
    # isov_contrib[k] =  3 / 2  * Pi_light(QVAL**2, aZ, Mz, particle_list, GG, qq, nloops, QED) # without charge factor

# mean_tot = np.mean(isov_contrib) 
# std_tot = np.std(isov_contrib) 
# print(" -At Q^2 ", QVAL**2 , " GeV^2, with ", nloops, " loops: ", mean_tot, std_tot, "\n")

# saving data

adler_mean = np.mean(adler_tot, axis=1)
adler_err = np.std(adler_tot, axis=1)

with open("/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/deltalpha-analysis/new_strategy/pert_running/comparison_harvey_code/data/adler_light_4loop.dat", "w") as ff:
    print("#idx    q_gev    adler_mean  adler_err", file=ff)
    count=0
    for qq, admean, aderr in zip(q_list, adler_mean, adler_err):
        print("%d %e  %e  %e" % (count, qq, admean, aderr), file = ff)
        count+=1
