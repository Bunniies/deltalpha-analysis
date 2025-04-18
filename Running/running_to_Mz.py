from adlerpy.adler_routines import Particle
from adlerpy.adler_sm import adler_charm_pert # The charm quark connected contribution to the Adler function
from adlerpy.adler_sm import adler_bottom_pert # The bottom quark connected contribution to the Adler function
from adlerpy.adler_sm import adler_light_pert # The light quark connected contribution to the Adler function
from adlerpy.adler_sm import adler_OZI_pert # The quark disconnected contribution to the Adler function
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

# inputs with no error
Mz = 91.1876 # PDG
aZ_config = bootstrap.bootstrap(0.1183, 0.0007, SAMPLES, NBIN ) # alpha_z at Z pole from FLAG24

# MSbar quark masses
mUp_config = bootstrap.bootstrap(2.16*0.001, 0.07*0.001, SAMPLES, NBIN)     # PDG
mDown_config = bootstrap.bootstrap(4.7*0.001, 0.07*0.001, SAMPLES, NBIN)    # PDG
mStrange_config = bootstrap.bootstrap(93.5*0.001, 0.8*0.001, SAMPLES, NBIN) # PDG
mC_config = bootstrap.bootstrap(1.278, 0.006, SAMPLES, NBIN) #  MS_c (FLAG)
mB_config = bootstrap.bootstrap(4.171, 0.02, SAMPLES, NBIN)  # MS_B (FLAG)

# heavy quark pole mass
mC_Pole_config = bootstrap.bootstrap(1.67, 0.07, SAMPLES, NBIN) # scale invariant mC (PDG)
mB_Pole_config = bootstrap.bootstrap(4.78, 0.06, SAMPLES, NBIN) # scale invariant mC (PDG)

# Condensates
GG_config = bootstrap.bootstrap(0.012, 0.012, SAMPLES, NBIN)  # gluon condensate https://arxiv.org/pdf/2302.01359
QQ_config = bootstrap.bootstrap(0.0013, 0.007, SAMPLES, NBIN) # quark condensate https://arxiv.org/pdf/2302.01359


alpha0=1/137.035999180
nloops=5
 
charm_contrib  = np.empty(NBIN)
bottom_contrib = np.empty(NBIN)
light_contrib  = np.empty(NBIN)
disc_contrib   = np.empty(NBIN)

for k in range(0,NBIN):
    aZ = aZ_config[k]
    mup = mUp_config[k]
    mdw = mDown_config[k]
    mst = mStrange_config[k]
    mc0 = mC_config[k]
    mb0 = mB_config[k]
    mc0_pole = mC_Pole_config[k]
    mb0_pole = mB_Pole_config[k]

    qq=QQ_config[k];
    GG=GG_config[k];
   
    Intervals_charm  = [np.sqrt(4),2*mc0*1.3,8,9,12,15,Mz]
    Intervals_bottom = [np.sqrt(4),2*mb0*0.7,2*mb0*1.3,15,20,30,Mz]
    Intervals_light  = [np.sqrt(4),9,12,15,Mz]
    mpole_on = True

    #Here you define the SM particles 
    up=Particle("up",x=[mup,2,3],mudec=0.001,mpole=None,mpole_on=False)
    down=Particle("down",x=[mdw,2,3],mudec=0.001,mpole=None,mpole_on=False)
    strange=Particle("strange",x=[mst,2,3],mudec=0.001,mpole=None,mpole_on=False)
    charm=Particle("charm",x=[mc0,mc0,4],mudec=2*(mc0),mpole=mc0_pole,mpole_on=mpole_on)
    bottom=Particle("bottom",x=[mb0,mb0,5],mudec=2*mb0,mpole=mb0_pole,mpole_on=mpole_on)
    top=Particle("top",x=[164,164,6],mudec=164,mpole=164,mpole_on=False) # not necessary in this case.
    particle_list=[up,charm,down,strange,bottom,top]

    nloops=5
    QED=True
    a=aZ/np.pi

    # define the function integrals
    def integrand_charm(Q,aZ,Mz,particles,GG,qq,nloops,QED, mpole_on):
        return 2*alpha0/3/np.pi*(adler_charm_pert(aZ=aZ,Mz=Mz,Q=Q,particles=particles,mpole_on=mpole_on,cut_low_as3=1.0,nloops=nloops,GG=GG,QED=QED))/Q;

    def integrand_bottom(Q,aZ,Mz,particles,GG,qq,nloops,QED,mpole_on):
        return 2*alpha0/3/np.pi*(adler_bottom_pert(aZ=aZ,Mz=Mz,Q=Q,particles=particles,mpole_on=mpole_on,cut_low_as3=1.3,nloops=nloops,GG=GG,QED=QED))/Q;

    def integrand_light(Q,aZ,Mz,particles,GG,qq,nloops,QED):
        return 2*alpha0/3/np.pi*(adler_light_pert(aZ=aZ,Mz=Mz,Q=Q,particles=particles,nloops=nloops,GG=GG,qq=qq,QED=QED))/Q;

    def integrand_disconnected(Q,aZ,Mz,particles,GG,qq,nloops,QED):
        return 2*alpha0/3/np.pi*(adler_OZI_pert(aZ=aZ,Mz=Mz,Q=Q,particles=particles,nloops=nloops,QED=QED,GG=GG,qq=qq))/Q;

    # compute charm contribution 
    Integrals=[]
    for i in range(len(Intervals_charm)-1):
        Integrals.append(quad(integrand_charm, 
                              Intervals_charm[i],
                              Intervals_charm[i+1],args=(aZ,Mz,particle_list,GG,qq,nloops,QED,mpole_on))[0])
    changecharm=1/137.036/np.pi*(4.203720393056578*a**2 + 17.148364322869213*a**3)
    dacharm=sum(Integrals)+changecharm
    charm_contrib[k] = dacharm 

    # compute bottom contribution
    Integrals=[]

    mb=bottom.mrun(aZ=aZ,Mz=Mz,mu=Mz,particles=particle_list)
    changebottom=1/137.036/np.pi*(-4*mb**2/3/Mz**2+1.0509300982641445*a**2 + 4.287091080717303*a**3)

    for i in range(len(Intervals_bottom)-1):
        Integrals.append(quad(integrand_bottom, 
                              Intervals_bottom[i],
                              Intervals_bottom[i+1],args=(aZ,Mz,particle_list,GG,qq,nloops,QED,mpole_on))[0])
    dabottom=sum(Integrals)+changebottom
    bottom_contrib[k] = dabottom 

    # compute light contribution
    Integrals=[]
    changelight=1/137.036/np.pi*(6.305580589584867*a**2 + 25.72254648430382*a**3)

    for i in range(len(Intervals_light)-1):
        Integrals.append(quad(integrand_light, 
                              Intervals_light[i],
                              Intervals_light[i+1],args=(aZ,Mz,particle_list,GG,qq,nloops,QED))[0])
    dalight=sum(Integrals)+changelight
    light_contrib[k] = dalight 

    # compute disconnected contribution 
    Integrals=[]

    for i in range(len(Intervals_light)-1):
        Integrals.append(quad(integrand_disconnected, 
                              Intervals_light[i],
                              Intervals_light[i+1],args=(aZ,Mz,particle_list,GG,qq,nloops,QED))[0])
    dasiconnected=sum(Integrals)
    disc_contrib[k] = dasiconnected

## Summing everything up 
print("Contribution per flavour")
mean_charm = np.mean(charm_contrib)
std_charm = np.std(charm_contrib)
print(" -charm: ", mean_charm, std_charm)

mean_bottom = np.mean(bottom_contrib)
std_bottom =  np.std(bottom_contrib)
print(" -bottom: ", mean_bottom, std_bottom)

mean_light = np.mean(light_contrib)
std_light =  np.std(light_contrib)
print(" -light: ", mean_light, std_light)

mean_disc = np.mean(disc_contrib)
std_disc =  np.std(disc_contrib)
print(" -disc: ", mean_disc, std_disc)

print("Total contribution for interval starting at Q0^2=", Intervals_charm[0])
mean_tot = mean_charm + mean_bottom + mean_light + mean_disc
std_tot = np.sqrt(std_charm**2 + std_bottom**2 + std_light**2 + std_disc**2)
print(" -Total: ", mean_tot, std_tot)

print("Total contribution from bootstrap")
total_contrib = charm_contrib + bottom_contrib + light_contrib + disc_contrib
mean_tot = np.mean(total_contrib)
std_tot =  np.std(total_contrib)
print(" -Total: ", mean_tot, std_tot)

