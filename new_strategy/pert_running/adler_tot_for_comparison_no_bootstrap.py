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
aZ = 0.1183  # FLAG 2024
# daZ = 0.0007 # FLAG 2024
daZ = 0.0016 # Old error
# MSbar quark masses
mUp = 2.16*0.001      # [GeV] PDG
mDown = 4.7*0.001     # [GeV] PDG
mStrange = 93.5*0.001 # [GeV] PDG

mC = 1.278  # scale invariant mC (FLAG)
dmC = 0.006 # scale invariant mC (FLAG)
mB = 4.171  # scale invariant mB (FLAG)
dmB = 0.02  # scale invariant mB (FLAG)

mC_Pole= 1.67   # scale invariant mC (PDG)
dmC_Pole = 0.07 # scale invariant mC (PDG)
mB_Pole = 4.78  # scale invariant mB (PDG)
dmB_Pole = 0.06 # scale invariant mB (PDG)

# Condensates
GG =0.012   # gluon condensate https://arxiv.org/pdf/2302.01359
dGG = 0.012
QQ = -0.0013  # quark condensate https://arxiv.org/pdf/2302.01359
dQQ = 0.007
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

  
nloops=5 # perturbative order alphas^5
mulow=1.273     
mpole_on = False

#Here you define the SM particles 
up=Particle("up",x=[mUp,2,3],mudec=0.001,mpole=None,mpole_on=False)
down=Particle("down",x=[mDown,2,3],mudec=0.001,mpole=None,mpole_on=False)
strange=Particle("strange",x=[mStrange,2,3],mudec=0.001,mpole=None,mpole_on=False)
charm=Particle("charm",x=[mC,mC,4],mudec=2*(mC),mpole=mC_Pole,mpole_on=False)
bottom=Particle("bottom",x=[mB,mB,5],mudec=2*mB,mpole=mB_Pole,mpole_on=False)
top=Particle("top",x=[164,164,6],mudec=164,mpole=164,mpole_on=False) # not necessary in this case.
particle_list=[up,charm,down,strange,bottom,top]

charm_pluserr=Particle("charm",x=[mC+dmC,mC+dmC,4],mudec=2*(mC+dmC),mpole=mC_Pole+dmC_Pole,mpole_on=False)
bottom_pluserr=Particle("bottom",x=[mB+dmB,mB+dmB,5],mudec=2*(mB+dmB),mpole=mB_Pole+dmB_Pole,mpole_on=False)
particle_list_pluserr=[up,charm_pluserr,down,strange,bottom_pluserr,top]


adler_tot = np.empty(len(q_list))
adler_tot_err = np.empty(len(q_list))


for q in range(0,len(q_list)):
    adcentral = adlertot(aZ, Mz, q_list[q], particles=particle_list, nloops=nloops,mpole_on=mpole_on,mulow=mulow,GG=GG,qq=QQ)
    ad_plus_as = adlertot(aZ+daZ, Mz, q_list[q], particles=particle_list, nloops=nloops,mpole_on=mpole_on,mulow=mulow,GG=GG,qq=QQ)
    ad_plus_GG = adlertot(aZ, Mz, q_list[q], particles=particle_list, nloops=nloops,mpole_on=mpole_on,mulow=mulow,GG=GG+dGG,qq=QQ)
    ad_plus_QQ = adlertot(aZ, Mz, q_list[q], particles=particle_list, nloops=nloops,mpole_on=mpole_on,mulow=mulow,GG=GG,qq=QQ+dQQ)
    ad_plus_pert = adlertot(aZ, Mz, q_list[q], particles=particle_list, nloops=nloops-1,mpole_on=mpole_on,mulow=mulow,GG=GG,qq=QQ)
    ad_plus_masses = adlertot(aZ, Mz, q_list[q], particles=particle_list_pluserr, nloops=nloops,mpole_on=mpole_on,mulow=mulow,GG=GG,qq=QQ)

    # error
    error = np.sqrt(
        (ad_plus_as - adcentral)**2 +
        (ad_plus_GG - adcentral)**2 +
        (ad_plus_QQ - adcentral)**2 +
        (ad_plus_pert - adcentral)**2 + 
        (ad_plus_masses - adcentral)**2
    )

    adler_tot[q] = adcentral
    adler_tot_err[q] = error


# saving data


with open("/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/running_jegerlhener/running_AdlerPy/using_my_input_values/log_spaced_mom_no_bootstrap/adlerPy.dat", "w") as ff:
    print("#idx    q_gev    adler_mean  adler_err", file=ff)
    count=0
    for qq, admean, aderr in zip(q_list, adler_tot, adler_tot_err):
        print("%d %e  %e  %e" % (count, qq, admean, aderr), file = ff)
        count+=1


