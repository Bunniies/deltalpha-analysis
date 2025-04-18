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

Mz = 91.1876 
# aZ_config = bootstrap.bootstrap(0.1179, 0.0010, SAMPLES, NBIN ) # PDG 2020 (KM)
aZ_config = bootstrap.bootstrap(0.1183, 0.0007, SAMPLES, NBIN ) # FLAG 2024

# MSbar quark masses
mUp = 2.16*0.001      # [GeV] PDG
# mUp = 2.41*0.001      # [GeV] KM
mDown = 4.7*0.001     # [GeV] PDG
# mDown = 4.4*0.001     # [GeV] KM
mStrange = 93.5*0.001 # [GeV] PDG
# mStrange = 85.2*0.001 # [GeV] KM

mC_config = bootstrap.bootstrap(1.278, 0.006, SAMPLES, NBIN) # scale invariant mC (FLAG)
# mC_config = [1.3 for _ in range(0,NBIN)] # scale invariant mC (KM)
mB_config = bootstrap.bootstrap(4.171, 0.02, SAMPLES, NBIN) # scale invariant mB (FLAG)
# mB_config = [4.2 for _ in range(0,NBIN)] # scale invariant mB (KM)

# heavy quark pole mass
mC_Pole_config = bootstrap.bootstrap(1.67, 0.07, SAMPLES, NBIN) # scale invariant mC (PDG)
# mC_Pole_config = bootstrap.bootstrap(1.35, 0.08, SAMPLES, NBIN) # scale invariant mC (KM)
mB_Pole_config = bootstrap.bootstrap(4.78, 0.06, SAMPLES, NBIN) # scale invariant mC (PDG)
# mB_Pole_config = bootstrap.bootstrap(4.40, 0.12, SAMPLES, NBIN) # scale invariant mC (KM)

# Condensates
GG_config = bootstrap.bootstrap(0.012, 0.012, SAMPLES, NBIN)  # gluon condensate https://arxiv.org/pdf/2302.01359
QQ_config = bootstrap.bootstrap(0.0013, 0.007, SAMPLES, NBIN) # quark condensate https://arxiv.org/pdf/2302.01359

alpha0=1/137.035999180
mmu=0.105658 # muon mass


def adlertot(aZ,Mz,Q,particles,nloops,mpole_on, mulow,GG,qq):
    adlight=adler_light_pert(aZ=aZ,Mz=Mz,Q=Q,particles=particles,nloops=nloops,GG=GG,qq=qq,QED=True)
    adc=adler_charm_pert(aZ=aZ,Mz=Mz,Q=Q,particles=particles,mulow=mulow,mpole_on=mpole_on,nloops=nloops,GG=GG,QED=True)
    adb = adler_bottom_pert(aZ=aZ, Mz=Mz,Q=Q, particles=particles, mpole_on=mpole_on, cut_low_as3=1.3, nloops=nloops, GG=GG, QED=True)
    addisc = adler_OZI_pert(aZ=aZ, Mz=Mz,Q=Q, particles=particles,nloops=nloops, QED=True, GG=GG, qq=qq)
    return adc+adlight+adb+addisc

emin = 1.0
emax = 100.0
q_list_log = np.linspace(np.log(emin), np.log(emax), num=400) # log(Q) in GeV
q_list = np.exp(q_list_log) 
# q_list = np.linspace(1.0, 100, num=400) # Q in GeV

adler_tot = np.empty((len(q_list), NBIN))

for k in range(0, NBIN):
    aZ = aZ_config[k]
    mc0 = mC_config[k]
    mb0 = mB_config[k]
    mc0_pole = mC_Pole_config[k]
    mb0_pole = mB_Pole_config[k]
    GG=GG_config[k];
    qq=QQ_config[k];
    nloops=4 # perturbative order alphas^5
    mulow=1.273     
    mpole_on = False


    #Here you define the SM particles 
    up=Particle("up",x=[mUp,2,3],mudec=0.001,mpole=None,mpole_on=False)
    down=Particle("down",x=[mDown,2,3],mudec=0.001,mpole=None,mpole_on=False)
    strange=Particle("strange",x=[mStrange,2,3],mudec=0.001,mpole=None,mpole_on=False)
    charm=Particle("charm",x=[mc0,mc0,4],mudec=2*(mc0),mpole=mc0_pole,mpole_on=False)
    bottom=Particle("bottom",x=[mb0,mb0,5],mudec=2*mb0,mpole=mb0_pole,mpole_on=False)
    top=Particle("top",x=[164,164,6],mudec=164,mpole=164,mpole_on=False) # not necessary in this case.
    particle_list=[up,charm,down,strange,bottom,top]

    for q in range(0,len(q_list)):
        ad_q = adlertot(aZ, Mz, q_list[q], particles=particle_list, nloops=nloops, mpole_on=mpole_on, mulow=mulow, GG=GG, qq=qq)
        adler_tot[q,k] = ad_q


# saving data

adler_mean = np.mean(adler_tot, axis=1)
adler_err = np.std(adler_tot, axis=1)



with open("/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/running_jegerlhener/running_AdlerPy/using_my_input_values/log_spaced_mom/adlerPy_4loops.dat", "w") as ff:
    print("#idx    q_gev    adler_mean  adler_err", file=ff)
    count=0
    for qq, admean, aderr in zip(q_list, adler_mean, adler_err):
        print("%d %e  %e  %e" % (count, qq, admean, aderr), file = ff)
        count+=1


