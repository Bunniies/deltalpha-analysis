
from adlerpy.adler_routines import Particle
from adlerpy.adler_sm import adler_charm_pert # The charm quark connected contribution to the Adler function
from adlerpy.adler_sm import adler_bottom_pert # The bottom quark connected contribution to the Adler function
from adlerpy.adler_sm import adler_light_pert # The light quark connected contribution to the Adler function
from adlerpy.adler_sm import adler_OZI_pert # The quark disconnected contribution to the Adler function
from adlerpy.adler_sm import adler_massless_connected
from scipy.optimize import fsolve
from adlerpy.adler_routines import adler_o_masless_i_heavy
from adlerpy.adler_routines import alphas
from adlerpy.adler_routines import adler_massless_connected

from adlerpy.adler_routines import Particle
from scipy.integrate import quad
from matplotlib import pyplot as plt 
import numpy as np 
import pandas as pd



import bootstrap
import statistics

#Here you define the SM particles 
up=Particle("up",x=[2.16*0.001,2,3],mudec=0.001,mpole=None,mpole_on=False)
down=Particle("down",x=[4.7*0.001,2,3],mudec=0.001,mpole=None,mpole_on=False)
strange=Particle("strange",x=[0.0935,2,3],mudec=0.001,mpole=None,mpole_on=False)
charm=Particle("charm",x=[1.278,1.278,4],mudec=2*1.278,mpole=1.67,mpole_on=False)
bottom=Particle("bottom",x=[4.171,4.171,5],mudec=2*4.171,mpole=4.25,mpole_on=False)
top=Particle("top",x=[164,164,6],mudec=164,mpole=164,mpole_on=False)
particle_list=[up,charm,down,strange,bottom,top] 


GG=0.012;  # Gluon condensate
dGG=0.012; # Gluon condensate error 100%
qq=-0.0013 # quark condensate
dqq=0.0007 # quark condensate error
aZ=0.1183; # central value
Mz=91.1876
# daZ=0.0007 # new error
daZ=0.0016 # old error

nloops=4;  # perturbative order alphas^5
mulow=1.273     

def adlertot(aZ,Mz,Q,particles,nloops, mulow,GG,qq):
    adlight=adler_light_pert(aZ=aZ,Mz=Mz,Q=Q,particles=particles,nloops=nloops,GG=GG,qq=qq,QED=True)
    adc=adler_charm_pert(aZ=aZ,Mz=Mz,Q=Q,particles=particles,mulow=mulow,nloops=nloops,GG=GG,QED=True)
    adb = adler_bottom_pert(aZ=aZ, Mz=Mz,Q=Q, particles=particles, cut_low_as3=1.3, nloops=nloops, GG=GG, QED=True)
    addisc = adler_OZI_pert(aZ=aZ, Mz=Mz,Q=Q, particles=particles,nloops=nloops, QED=True, GG=GG, qq=qq)
    return adc, adlight, adb, addisc

emin = 1.0 #1.0
emax = 100.0
q_list_log = np.linspace(np.log(emin), np.log(emax), num=400) # log(Q) in GeV
q_list = np.exp(q_list_log) 
q_list[391] = Mz 

adler_central_l = np.empty(len(q_list)) # light
adler_central_c = np.empty(len(q_list)) # charm
adler_central_b = np.empty(len(q_list)) # bottom
adler_central_ozi = np.empty(len(q_list)) # disc

adler_error_l = np.empty(len(q_list)) # light
adler_error_c = np.empty(len(q_list)) # charm
adler_error_b = np.empty(len(q_list)) # bottom
adler_error_ozi = np.empty(len(q_list)) # disc

for q in range(0, len(q_list)):
    adc, adl, adb, adozi = adlertot(aZ,Mz,q_list[q],particles=particle_list,nloops=nloops,mulow=mulow,GG=GG,qq=qq) # central value
    adcplusas, adlplusas, adbplusas, adoziplusas = adlertot(aZ+daZ,Mz,q_list[q],particles=particle_list,nloops=nloops,mulow=mulow,GG=GG,qq=qq) # shift from alphas
    adcpluspert, adlpluspert, adbpluspert, adozipluspert = adlertot(aZ,Mz,q_list[q],particles=particle_list,nloops=nloops-1,mulow=mulow,GG=GG,qq=qq) # truncation error
    adcplusGG, adlplusGG, adbplusGG, adoziplusGG =adlertot(aZ,Mz,q_list[q],particles=particle_list,nloops=nloops,mulow=mulow,GG=GG+dGG,qq=qq) # parametric shift GG
    adcplusqq, adlplusqq, adbplusqq, adoziplusqq =adlertot(aZ,Mz,q_list[q],particles=particle_list,nloops=nloops,mulow=mulow,GG=GG,qq=qq+dqq) # parametric shift qq

    errorl = np.sqrt((adl - adlplusas)**2 + (adl - adlpluspert)**2 + (adl - adlplusGG)**2 + (adl - adlplusqq)**2)
    errorc = np.sqrt((adc - adcplusas)**2 + (adc - adcpluspert)**2 + (adc - adcplusGG)**2 + (adc - adcplusqq)**2)
    errorb = np.sqrt((adb - adbplusas)**2 + (adb - adbpluspert)**2 + (adb - adbplusGG)**2 + (adb - adbplusqq)**2)
    errorozi = np.sqrt((adozi - adoziplusas)**2 + (adozi - adozipluspert)**2 + (adozi - adoziplusGG)**2 + (adozi - adoziplusqq)**2)

    adler_central_l[q] = adl
    adler_central_c[q] = adc
    adler_central_b[q] = adb
    adler_central_ozi[q] = adozi

    adler_error_l[q] = errorl
    adler_error_c[q] = errorc
    adler_error_b[q] = errorb
    adler_error_ozi[q] = errorozi

# saving data
import os
path_save = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/running_jegerlhener/running_AdlerPy/adler_per_comp/rodolfo_err"
fname = "adlerPy_rodolfo_old_alpha_err.dat"

with open(os.path.join(path_save, fname), "w") as ff:
        print("# idx    q_gev    adl_val    adl_err       adc_val      adc_err     adb_val      adb_err     adozi_val      adozi_err ", file=ff)
        count = 0
        for q in range(0,len(q_list)): 
            print("%d %e %e %e %e %e %e %e %e %e " % (count, q_list[q], adler_central_l[q], adler_error_l[q], adler_central_c[q], adler_error_c[q], adler_central_b[q], adler_error_b[q], adler_central_ozi[q], adler_error_ozi[q]), file=ff)
            count +=1
