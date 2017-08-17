#!/usr/bin/python



# -*- coding: utf-8 -*-

#Created on Mon Jan 25 13:22:33 2016

#@author: marisa
#Standalone_Taufit_simu.py

import argparse
import datetime
import os, sys
import pypsr_standalone as psr
import matplotlib.pyplot as plt
import numpy as np

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

#"""Read and print header information"""
"""Define options to the script"""
parser = argparse.ArgumentParser()
parser.add_argument('-sx','--spectralindex',type=float,
                    help="Provide the spectral index alpha, where tau propto freq^-alpha")
parser.add_argument('-fr','--freq',nargs = '*', type=float,
                    help="Provide the frequency range in MHz")
parser.add_argument('-nb','--nbins',type=int,
                    help="Provide the number of bins (interger)")
parser.add_argument('-nch','--nchan',type=int,
                    help="Provide the number of frequency channels (interger)")
parser.add_argument('-m','--method',
                    help="Choosing scattering method. Choose between 'iso' and 'onedim'")                   
parser.add_argument('-dt','--datatag',
                    help="Tag written into saved txt files")
parser.add_argument('-snr','--snr',type=float,
                    help="Provide signal to noise ratio")
parser.add_argument('-dc','--dutycycle', type=float,
                    help="Provide the dutycycle of pulse as percentage")
parser.add_argument('-tatf','--tauatfrequency',nargs='*', type=float,
                    help="Provide a tau value (sec) at a given frequency (MHz), e.g. 0.5 150")
parser.add_argument('-p','--pulseperiod', type=float,
                    help="Provide the pulse period (sec)")                
                   

args = parser.parse_args()

"""Allocate variable names to the parsed options"""
alpha= args.spectralindex
freqlow = args.freq[0]
freqhigh = args.freq[1]
pulseperiod = args.pulseperiod
meth = args.method
nbins = args.nbins
nch = args.nchan
datat = args.datatag
snr = args.snr
dutyc = args.dutycycle
tatf_tau = args.tauatfrequency[0]
tatf_freq = args.tauatfrequency[1]
tatf_freqGHz = tatf_freq/1000.

if datat is None:
    datat = 'Simualted'

"""Create folder to save to"""
newpath = r'./SimulatedData'
if not os.path.exists(newpath):
    os.makedirs(newpath)
    
freqsimu = np.linspace(freqlow,freqhigh,nch) #(in MHz)
freqsimuGHz = freqsimu/1000.

"""Calculate tau value at each frequency channel, using
psr.tauatfreq(oldfreq,oldtau,newfreq,specindex)"""  
"""Assumes constant specindex"""

tausecs_simu = []
for i in freqsimuGHz:
    tausec = psr.tauatfreq(tatf_freqGHz,tatf_tau,i,alpha)
    tausecs_simu.append(tausec)
    

pulsar = 'Simulated'
print0 = "Pulsar : %s" %pulsar
print1 = "Number of channels: %d" %nch
print2 = "Number of bins: %d" %nbins
print3 = "Taus (ms):"

for k in range(4):
    print eval('print{0}'.format(k))

for j in range(nch):
    print "Freq %.1f MHz: %.2f" %(freqsimu[j], tausecs_simu[j]*1000) 
   
profilexaxis = np.linspace(0,pulseperiod,nbins)

freqMHz = freqsimu
freqGHz = freqMHz/1000.

fluxspecindx = 1.6
print "Using a flux spectral index of -%.1f to scale the intrinsic Gaussian amplitude" %fluxspecindx

datas = []
profiles =[]

for i in range(nch):
    profile, data = psr.simulate(pulseperiod,tausecs_simu[i],dutyc,fluxspecindx,freqGHz[i],freqlow/1000.,nbins,snr,meth)
    comp_rms = psr.find_rms(data,nbins)
    datas.append(data)
    profiles.append(profile)

   
"""Plotting starts"""

if meth in 'onedim':
    print "METH IS ONEDIM"
    datacol = 'b'
if meth in 'iso':
    print "METH IS ISO"
    datacol = 'r'
else:
    print "NO METH SPECIFIED"
    
##PLOT DATA##

plt.close('all')

figg = plt.figure(1,figsize=(16,10))
figg.subplots_adjust(left = 0.055, right = 0.98, wspace=0.35,hspace=0.35,bottom=0.05)    
for j in range(nch):
    plt.subplot(3,4,j+1)
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.plot(profilexaxis,datas[j],datacol,alpha = 0.25)
#    plt.plot(profilexaxis,profiles[j],'m', alpha = 0.4,lw=2.0)
    plt.title('PSR %s at %.1f MHz' %(pulsar, freqMHz[j]))
    plt.xlim(xmax=pulseperiod)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.xlabel('time (sec)')
    plt.ylabel('intensity (mJy)')
#    plt.ylim(ymax=1.0)
        


plt.show()


RMSs = []
for i in range(nch):
    peak = np.max(datas[i])
    rms = peak/snr
    RMSs.append(rms)
meanRMS = np.mean(RMSs)

"""Save txt files"""
    
timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H:%M:%S")

txtpath = r'./SimulatedDatafiles'
if not os.path.exists(txtpath):
    os.makedirs(txtpath)
txtfile = '%s_%s_data_%s.txt' %(pulsar,meth,timestamp)
FPtxt = os.path.join(txtpath,txtfile) 
headr0 = 'File: %s Src: %s Nsub: %d Nch: %d Npol: %d Nbin: %d RMS: %.2f' %('SimulatedFile',pulsar,1,nch,1,nbins,meanRMS)

#1. Col 0s, Freqch, bins, intensity
col0 = np.zeros(nbins*nch,dtype=int)
colfre = []
colnbs = []
colnb = np.linspace(1,nbins,nbins,dtype=int)
for i in range(nch):
    colf = i*np.ones(nbins,dtype=int)
    colfre.append(colf)
    colnbs.append(colnb)
col1 = np.array(colfre).reshape(nch*nbins,1).flatten()
col2 = np.array(colnbs).reshape(nch*nbins,1).flatten()
col3 = np.array(datas).reshape(nch*nbins,1).flatten()
colarray = np.array([col0,col1,col2,col3]).T
np.savetxt(FPtxt,colarray,fmt=['%d','%d','%d','%f'],header=headr0)

bw = (freqhigh-freqlow)/nch


f = open(FPtxt, "r")
contents = f.readlines()
f.close()
for i in range(nch):
    headr1 = '# MJD(mid): %s P: %f Freq: %f BW: %f \n' %('Simulated',pulseperiod,freqMHz[i],bw)
    contents.insert((i+1) + nbins*i, headr1)
f = open(FPtxt, "w")
contents = "".join(contents)
f.write(contents)
f.close()


print "%s saved in %s" %(txtfile,txtpath)



####FIT MARK DATA
#markt = np.loadtxt("/Users/marisa/python-workingdir/LC6_030/B1937+21/MarkData/data2col.txt")
#markunlog = 10**(markt[:,1]/2)
#from lmfit import Model, conf_interval, printfuncs
#from lmfit import minimize, Parameter, Parameters, fit_report
#from lmfit.models import LinearModel, PowerLawModel, ExponentialModel, QuadraticModel
#markbins = markt[:,0]
#markbins = np.delete(markbins,[2020,2044])
#markunlog = np.delete(markunlog,[2020,2044])
#expmod = ExponentialModel()
#exppars = expmod.guess(markunlog[2000:-1], x=markbins[2000:-1])
#expout = expmod.fit(markunlog[2000:-1], exppars, x=markbins[2000:-1])
#plt.plot(markbins[2000:-1],expout.best_fit,'k-'
