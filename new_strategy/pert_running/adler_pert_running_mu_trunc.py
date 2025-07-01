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
SAMPLES = 10000
NBIN = 100

## input values 
Mz = 91.1876 

aZ_val = 0.1183
# aZ_err = 0.0007 # new err
aZ_err = 0.0016 # old err
aZ_config = bootstrap.bootstrap(aZ_val, aZ_err, SAMPLES, NBIN ) # FLAG 2024

mUp = 2.16*0.001      # [GeV] PDG
mDown = 4.7*0.001     # [GeV] PDG
mStrange = 93.5*0.001 # [GeV] PDG

mC_mean = 1.278
mC_err = 0.006
mC_config = bootstrap.bootstrap(mC_mean, mC_err, SAMPLES, NBIN)

mB_mean = 4.171
mB_err = 0.02
mB_config = bootstrap.bootstrap(mB_mean, mB_err, SAMPLES, NBIN)

GG_val = 0.012
GG_err = 0.012
GG_config = bootstrap.bootstrap(GG_val, GG_err, SAMPLES, NBIN)

QQ_val = -0.0013
QQ_err = 0.0007
QQ_config = bootstrap.bootstrap(QQ_val, QQ_err, SAMPLES, NBIN)

alpha0=1/137.035999180 # alpha thompson limit
muon_m=0.105658 # muon mass

# particle list for central value of Adler function
up_=Particle("up",x=[mUp,2,3],mudec=0.001,mpole=None,mpole_on=False)
down_=Particle("down",x=[mDown,2,3],mudec=0.001,mpole=None,mpole_on=False)
strange_=Particle("strange",x=[mStrange,2,3],mudec=0.001,mpole=None,mpole_on=False)
charm_=Particle("charm",x=[mC_mean,mC_mean,4],mudec=2*(mC_mean),mpole=None,mpole_on=False)
bottom_=Particle("bottom",x=[mB_mean,mB_mean,5],mudec=2*mB_mean,mpole=None,mpole_on=False)
top_=Particle("top",x=[164,164,6],mudec=164,mpole=None,mpole_on=False) # not necessary in this case.
particle_list_centra_val=[up_,charm_,down_,strange_,bottom_,top_]

def adlertot(aZ,Mz,Q,particles,mu,nloops,mpole_on, mulow,GG,qq):
    adlight=adler_light_pert(aZ=aZ,Mz=Mz,Q=Q,particles=particles,mu=mu,nloops=nloops,GG=GG,qq=qq,QED=True)
    adc=adler_charm_pert(aZ=aZ,Mz=Mz,Q=Q,particles=particles,mu=mu,mulow=mulow,mpole_on=mpole_on,nloops=nloops,GG=GG,QED=True)
    adb = adler_bottom_pert(aZ=aZ, Mz=Mz,Q=Q, particles=particles,mu=mu, mpole_on=mpole_on, cut_low_as3=1.3, nloops=nloops, GG=GG, QED=True)
    addisc = adler_OZI_pert(aZ=aZ, Mz=Mz,Q=Q, particles=particles,mu=mu,nloops=nloops, QED=True, GG=GG, qq=qq)
    return adc, adlight, adb, addisc

emin = 1.0 #1.0
emax = 100.0
q_list_log = np.linspace(np.log(emin), np.log(emax), num=400) # log(Q) in GeV
q_list = np.exp(q_list_log) 
q_list[391] = Mz # setting to mZ the q_list to ease the integration.
# q_list = [1.8224332442626667]

adler_bins_l = np.empty((len(q_list), NBIN))
adler_bins_c = np.empty((len(q_list), NBIN))
adler_bins_b = np.empty((len(q_list), NBIN))
adler_bins_ozi = np.empty((len(q_list), NBIN))

# storing central values
adler_central_l = np.empty(len(q_list)) # light
adler_central_c = np.empty(len(q_list)) # charm
adler_central_b = np.empty(len(q_list)) # bottom
adler_central_ozi = np.empty(len(q_list)) # disc

nloops=5
mulow=1.273
mpole_on=False

for q in range(0,len(q_list)):


    # Adler mean values
    adc, adlight, adb, adozi = adlertot(aZ_val, Mz, q_list[q],particles=particle_list_centra_val,mu=None,nloops=nloops,mpole_on=mpole_on, mulow=mulow, GG=GG_val, qq=QQ_val)
    
    adler_central_c[q] = adc
    adler_central_l[q] = adlight
    adler_central_b[q] = adb
    adler_central_ozi[q] = adozi

    # taking care of decoupling
    if q_list[q]*np.sqrt(2) < 2*mC_mean:
        mu_config = np.random.uniform(q_list[q]/np.sqrt(2), q_list[q]*np.sqrt(2), NBIN)

    elif q_list[q]*np.sqrt(2) > 2*mC_mean and q_list[q]/np.sqrt(2) < 2*mC_mean:
        if q_list[q] < 2*mC_mean:
            mu_config = np.random.uniform(q_list[q]/np.sqrt(2), 2*mC_mean, NBIN)
        else:
            mu_config = np.random.uniform(2*mC_mean, q_list[q]*np.sqrt(2), NBIN)

    elif q_list[q]*np.sqrt(2) < 2*mB_mean and q_list[q]/np.sqrt(2) > 2*mC_mean:
        mu_config = np.random.uniform(q_list[q]/np.sqrt(2), q_list[q]*np.sqrt(2), NBIN)

    elif q_list[q]*np.sqrt(2) > 2*mB_mean and q_list[q]/np.sqrt(2) < 2*mB_mean:
        if q_list[q] < 2*mB_mean:
            mu_config = np.random.uniform(q_list[q]/np.sqrt(2), 2*mB_mean, NBIN)
        else:
            mu_config = np.random.uniform(2*mB_mean, q_list[q]*np.sqrt(2), NBIN)

    elif q_list[q]/np.sqrt(2) > 2*mB_mean:
        mu_config = np.random.uniform(q_list[q]/np.sqrt(2), q_list[q]*np.sqrt(2), NBIN)

    print("Q= ", q_list[q], " Q/sqrt(2)= ", q_list[q]/np.sqrt(2), " Q*sqrt(2)= ", q_list[q]*np.sqrt(2), "\n min(mu)= ", min(mu_config), " max(mu)= ", max(mu_config), "\n mudec_c= ", 2*mC_mean, " mudec_b= ", 2*mB_mean, "\n")
    for k in range(0,NBIN):
        
        mu = mu_config[k]
        aZ = aZ_config[k]
        mc0 = mC_config[k]
        mb0 = mB_config[k]
        GG = GG_config[k]
        qq = QQ_config[k]
         
        up=Particle("up",x=[mUp,2,3],mudec=0.001,mpole=None,mpole_on=False)
        down=Particle("down",x=[mDown,2,3],mudec=0.001,mpole=None,mpole_on=False)
        strange=Particle("strange",x=[mStrange,2,3],mudec=0.001,mpole=None,mpole_on=False)
        charm=Particle("charm",x=[mc0,mc0,4],mudec=2*mC_mean,mpole=None,mpole_on=False)
        bottom=Particle("bottom",x=[mb0,mb0,5],mudec=2*mB_mean,mpole=None,mpole_on=False)
        top=Particle("top",x=[164,164,6],mudec=164,mpole=None,mpole_on=False) # not necessary in this case.
        particle_list=[up,charm,down,strange,bottom,top]
        # try: 
        adc_bin, adl_bin, adb_bin, adozi_bin = adlertot(aZ, Mz, q_list[q], particles=particle_list, mu=mu, nloops=nloops,mpole_on=mpole_on, mulow=mulow,GG=GG, qq=qq)

        adler_bins_c[q,k] = adc_bin
        adler_bins_l[q,k] = adl_bin
        adler_bins_b[q,k] = adb_bin
        adler_bins_ozi[q,k] = adozi_bin
        # except Exception:

#
adler_err_l = np.std(adler_bins_l, axis=1)
adler_err_c = np.std(adler_bins_c, axis=1)
adler_err_b = np.std(adler_bins_b, axis=1)
adler_err_ozi = np.std(adler_bins_ozi, axis=1)

##
print("# Finished.")
print("# Nloops: ", nloops)
print("# mpole_on: ", mpole_on)


# saving data
import os
path_save = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/running_jegerlhener/running_AdlerPy/adler_per_comp/mu_as_trunc_err"
fname = "adlerPy_mutrunc_old_err.dat"

with open(os.path.join(path_save, fname), "w") as ff:
        print("# idx    q_gev    adl_val    adl_err       adc_val      adc_err     adb_val      adb_err     adozi_val      adozi_err ", file=ff)
        count = 0
        for q in range(0,len(q_list)): 
            print("%d %e %e %e %e %e %e %e %e %e " % (count, q_list[q], adler_central_l[q], adler_err_l[q], adler_central_c[q], adler_err_c[q], adler_central_b[q], adler_err_b[q], adler_central_ozi[q], adler_err_ozi[q]), file=ff)
            count +=1
