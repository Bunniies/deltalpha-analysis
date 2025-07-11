# perturbative running computed with bootstrap
# charm and light truncation error included
# results from this code are used in the paper 

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
# aZ_err = 0.0007 # new error
aZ_err = 0.0016 # old error
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
    adlight=adler_light_pert(aZ=aZ,Mz=Mz,Q=Q,particles=particles,mu=mu,nloops=nloops,GG=GG,qq=qq,QED=False)
    adc=adler_charm_pert(aZ=aZ,Mz=Mz,Q=Q,particles=particles,mu=mu, mpole_on=mpole_on,nloops=nloops,GG=GG,QED=False, cut_low_as3=1.0, cut_high_as3=20)
    adb = adler_bottom_pert(aZ=aZ, Mz=Mz,Q=Q, particles=particles,mu=mu, mpole_on=mpole_on, cut_low_as3=1.3, nloops=nloops, GG=GG, QED=False)
    addisc = adler_OZI_pert(aZ=aZ, Mz=Mz,Q=Q, particles=particles,mu=mu,nloops=nloops, QED=False, GG=GG, qq=qq)
    return adc, adlight, adb, addisc

emin = 1.0
emax = 100.0
q_list_log = np.linspace(np.log(emin), np.log(emax), num=400) # log(Q) in GeV
q_list = np.exp(q_list_log) 
q_list[391] = Mz # setting to mZ the q_list to ease the integration.

# q_list = [16, 20, 30, 50, 70, 80]

adler_bins_l = np.empty((len(q_list), NBIN))
adler_bins_c = np.empty((len(q_list), NBIN))
adler_bins_b = np.empty((len(q_list), NBIN))
adler_bins_ozi = np.empty((len(q_list), NBIN))

# storing central values
adler_central_l = np.empty(len(q_list)) # light
adler_central_c = np.empty(len(q_list)) # charm
adler_central_b = np.empty(len(q_list)) # bottom
adler_central_ozi = np.empty(len(q_list)) # disc

nloops=4
mulow=1.273
mpole_on=False

err_trunc_adl = np.empty(len(q_list))
err_trunc_adc = np.empty(len(q_list))
err_trunc_adb = np.empty(len(q_list))
err_trunc_adozi = np.empty(len(q_list))

for q in range(0,len(q_list)):


    # Adler mean values
    adc, adlight, adb, adozi = adlertot(aZ_val, Mz, q_list[q],particles=particle_list_centra_val,mu=None,nloops=nloops,mpole_on=mpole_on, mulow=mulow, GG=GG_val, qq=QQ_val)
    # truncation error light
    adc_trunc, adlight_trunc, adb_trunc, adozi_trunc = adlertot(aZ_val, Mz, q_list[q],particles=particle_list_centra_val,mu=None,nloops=nloops-1,mpole_on=mpole_on, mulow=mulow, GG=GG_val, qq=QQ_val)
    # truncation error charm
    adc_trunc, _, adb_trunc, adozi_trunc = adlertot(aZ_val, Mz, q_list[q],particles=particle_list_centra_val,mu=None,nloops=nloops-3,mpole_on=mpole_on, mulow=mulow, GG=GG_val, qq=QQ_val)

    adler_central_c[q] = adc
    adler_central_l[q] = adlight
    adler_central_b[q] = adb
    adler_central_ozi[q] = adozi

    err_trunc_adl[q] = adlight - adlight_trunc
    err_trunc_adc[q] = adc - adc_trunc
    err_trunc_adb[q] = adb - adb_trunc
    err_trunc_adozi[q] = adozi - adozi_trunc

    mu_config = np.random.uniform(q_list[q]/np.sqrt(2), q_list[q]*np.sqrt(2), NBIN)
    
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
        charm=Particle("charm",x=[mc0,mc0,4],mudec=2*(mc0),mpole=None,mpole_on=False)
        bottom=Particle("bottom",x=[mb0,mb0,5],mudec=2*mb0,mpole=None,mpole_on=False)
        top=Particle("top",x=[164,164,6],mudec=164,mpole=None,mpole_on=False) # not necessary in this case.
        particle_list=[up,charm,down,strange,bottom,top]

        adc_bin, adl_bin, adb_bin, adozi_bin = adlertot(aZ, Mz, q_list[q], particles=particle_list, mu=None, nloops=nloops,mpole_on=mpole_on, mulow=mulow,GG=GG, qq=qq)

        adler_bins_c[q,k] = adc_bin
        adler_bins_l[q,k] = adl_bin
        adler_bins_b[q,k] = adb_bin
        adler_bins_ozi[q,k] = adozi_bin
      


adler_err_l = np.sqrt(np.std(adler_bins_l, axis=1)**2 + err_trunc_adl**2)
adler_err_c = np.sqrt(np.std(adler_bins_c, axis=1)**2 + err_trunc_adc**2)
adler_err_b = np.sqrt(np.std(adler_bins_b, axis=1)**2 + err_trunc_adb**2)
adler_err_ozi = np.sqrt(np.std(adler_bins_ozi, axis=1)**2 + err_trunc_adozi**2)

##
print("# Finished.")
print("# Nloops: ", nloops)
print("# mpole_on: ", mpole_on)
print("# mu set to: ", "None ")
print("# QED set to False in light charm and bottom")


# saving data
import os
path_save = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/running_jegerlhener/running_AdlerPy/adler_per_comp/std_trunc_err"
fname = "adlerPy_trunc_err_included_alpha_old_err_QED_false.dat"

with open(os.path.join(path_save, fname), "w") as ff:
        print("# idx    q_gev    adl_val    adl_err       adc_val      adc_err     adb_val      adb_err     adozi_val      adozi_err ", file=ff)
        count = 0
        for q in range(0,len(q_list)): 
            print("%d %e %e %e %e %e %e %e %e %e " % (count, q_list[q], adler_central_l[q], adler_err_l[q], adler_central_c[q], adler_err_c[q], adler_central_b[q], adler_err_b[q], adler_central_ozi[q], adler_err_ozi[q]), file=ff)
            count +=1


## check Charm truncation error

# def adler_charm(aZ, Mz, Q, particles, mu, mpole_on, nloops, GG, QED):
    # adc=adler_charm_pert(aZ=aZ,Mz=Mz,Q=Q,particles=particles,mu=mu,mpole_on=mpole_on,nloops=nloops,GG=GG,QED=QED, cut_high_as3=20, cut_low_as3=1.0)
    # return adc
# 
# qval_c = 30
# adc = adler_charm(aZ_val, Mz, qval_c, particles=particle_list_centra_val, mu=None, nloops=3, mpole_on=mpole_on, GG=GG_val, QED=True)
# adc_0l = adler_charm(aZ_val, Mz, qval_c, particles=particle_list_centra_val, mu=None, nloops=0, mpole_on=mpole_on,  GG=GG_val, QED=True)
# adc_1l = adler_charm(aZ_val, Mz, qval_c, particles=particle_list_centra_val, mu=None, nloops=1, mpole_on=mpole_on,  GG=GG_val, QED=True)
# adc_2l = adler_charm(aZ_val, Mz, qval_c, particles=particle_list_centra_val, mu=None, nloops=2, mpole_on=mpole_on,  GG=GG_val, QED=True)
# 
# print("3loops: ", adc)
# print("2loops: ", adc_2l)
# print("1loops: ", adc_1l)
# print("0loops: ", adc_0l)