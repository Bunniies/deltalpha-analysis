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
QVAL = 8. # GeV
mmu=0.105658 # muon mass

isov_contrib = np.empty(NBIN)

def integrand_light(Q2,aZ,Mz,particles, GG, qq, nloops, QED):
    Q = np.sqrt(Q2)
    return 1 / (12*np.pi**2) * (adler_light_pert(aZ, Mz, Q=Q, particles=particle_list, GG=GG,qq=qq, nloops=nloops, QED=QED))/Q2;

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

    GG=0.0#GG_config[k];
    qq=0.0#QQ_config[k];
    nloops=5
    QED=True
    a=aZ/np.pi
    changelight=1/137.036/np.pi*(6.305580589584867*a**2 + 25.72254648430382*a**3)


    # isov_contrib[k] =   (16. * alpha0**2 * mmu**2 )/(9 * Q2VAL) * quad(integrand_light, Q2VAL/4, Q2VAL, args=(aZ,Mz,particle_list,GG,qq,nloops,QED))[0]
    # isov_contrib[k] =  0.5 / (137.036 * np.pi)**2 *  (16. * np.pi**2 * mmu**2 )/(9 * QVAL**2) * quad(integrand_light, QVAL**2/4, QVAL**2, args=(aZ,Mz,particle_list,GG,qq,nloops,QED))[0]
    isov_contrib[k] =  3 / 2 * 0.5 / (137.036 * np.pi)**2 *  (16. * np.pi**2 * mmu**2 )/(9 * QVAL**2) * Pi_light(QVAL**2, aZ, Mz, particle_list, GG, qq, nloops, QED)
    isov_contrib[k] *= 1e10

print("Final isov results at ", QVAL, " GeV2")
mean_tot = np.mean(isov_contrib)
std_tot = np.std(isov_contrib)
print(" -Total: ", mean_tot, std_tot)




## checks
adlight_up = adler_light_pert(aZ=aZ, Mz=Mz, Q=5, particles=particle_list, nloops=nloops, GG=GG, qq=qq, QED=True)

adlight_down = adler_light_pert(aZ=aZ, Mz=Mz, Q=2.5, particles=particle_list, nloops=nloops, GG=GG, qq=qq, QED=True)

print(adlight_up)
print(adlight_down)

## checks plot

def adlertot(aZ,Mz,Q,particles,nloops,GG,qq):
    adlight=adler_light_pert(aZ=aZ,Mz=Mz,Q=Q,particles=particles,nloops=nloops,GG=GG,qq=qq,QED=True)
    # adc=adler_charm_pert(aZ=aZ,Mz=Mz,Q=Q,particles=particles,mpole_on=False,nloops=nloops,GG=GG,QED=True)
    # adb=adler_bottom_pert(aZ=aZ,Mz=Mz,Q=Q,particles=particles,mpole_on=False,nloops=nloops,GG=GG,QED=True)
    # adozi=adler_OZI_pert(aZ,Mz,Q,particles,nloops=nloops,QED=QED,GG=GG,qq=qq)
    return adlight#+adc+adb+adozi

q_list = np.arange(1.2,90,0.05)
adlist=[adlertot(aZ=aZ,Mz=Mz,Q=q,particles=particle_list,nloops=5,GG=0,qq=0) for q in q_list]

fig, ax = plt.subplots()
#plt.plot(qplot,aslat,color="blue",linewidth=0.8)
plt.plot(q_list,adlist,color="blue",linewidth=0.9)
plt.grid()
plt.xticks(fontsize=14)
ax.spines['bottom'].set_linewidth(1.3)
ax.spines['left'].set_linewidth(1.3)
ax.spines['top'].set_linewidth(1.3)
ax.spines['right'].set_linewidth(1.3)
ax.set(xlabel='$Q^2$ (GeV)', ylabel='$\hat{\\alpha}_s(M_Z)$')
ax.set_xlim([0,20])
ax.set_ylim([2,4])
plt.show()


## conversion between MSbar and m pole 


# def check_conversion(asMq, Nl, mSbar):
    # LO = 1 + 4 / 3 * (asMq / np.pi)
    # NLO = (asMq / np.pi)**2 * (-1.0414 * Nl + 13.4434)
    # NNLO = (asMq / np.pi)**3 * (0.6527*Nl**2 - 26.655*Nl + 190.595)
    # return  mSbar * (LO + NLO + NNLO)
# 
# alpha_mu = alphas(0.1185, 91.1876, 1.278, particle_list)
# 
# mpole_converted = check_conversion(alpha_mu, 4, 1.278 )
# 
# print("alpha: ", alpha_mu)
# print("m pole: ", mpole_converted)