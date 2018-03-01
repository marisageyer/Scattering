#!/usr/bin/python



# -*- coding: utf-8 -*-

#Created on Mon Jan 25 13:22:33 2016

#@author: marisa
#Standalone_Taufit_simu.py

import argparse
import os, sys
import pypsr_standalone as psr
import matplotlib.pyplot as plt
import lmfit
from lmfit import Model, conf_interval, printfuncs
from lmfit import minimize, Parameter, Parameters, fit_report
from lmfit.models import LinearModel, PowerLawModel, ExponentialModel, QuadraticModel
import numpy as np
from scipy import special
from scipy import stats
import DataReadIn as dri

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

#"""Read and print header information"""
"""Define options to the script"""
parser = argparse.ArgumentParser()
parser.add_argument('-f','--filename',
                    help="Provide the pathway to the data file, containing first column of frequencies second column of tau values")

args = parser.parse_args()

"""Allocate variable names to the parsed options"""
filepath = args.filename

print "Expect format as produced by Taufitter.py (3 col: Freq(MHz), Tau(s), TauErr(s))"
pulsar, data = dri.read_Taufitter(filepath)


freqMHz = data[:,0]
freqGHz = data[:,0]/1000.
taus = data[:,1] 
if np.shape(data)[1] >= 3:
    tauserr = data[:,2]

taus_order, taus2_order = np.zeros(len(freqMHz)), np.zeros(len(freqMHz))
taus_order_err, taus2_order_err = np.zeros(len(freqMHz)), np.zeros(len(freqMHz))

if 'aniso' in filepath:
    tau2s = data[:,3] 
    tau2serr = data[:,4]
    
    for i in range(len(freqMHz)):
        if taus[i] > tau2s[i]:
            taus_order[i] = taus[i]
            taus2_order[i] = tau2s[i]
            taus_order_err[i] = tauserr[i]
            taus2_order_err[i] = tau2serr[i]     
        else:
            taus_order[i] = tau2s[i]
            taus2_order[i] = taus[i]
            taus_order_err[i] = tau2serr[i]   
            taus2_order_err[i] =  tauserr[i]

    taus = taus_order
    tau2s = taus2_order
    tauserr = taus_order_err
    tau2serr = taus2_order_err


print "Txt file imported"

"""CALCULATE FITS TO TAUS"""
powmod = PowerLawModel()
powparstau = powmod.guess(taus,x=freqGHz)

if np.shape(data)[1] >= 3:
    powout = powmod.fit(taus,powparstau,x=freqGHz,weights=1/(tauserr))
else:
    powout = powmod.fit(taus,powparstau,x=freqGHz,weights=1)
print(fit_report(powout.params))    
fit = powout.best_fit    
alpha = -powout.best_values['exponent']
alphaerr = powout.params['exponent'].stderr
amp = powout.best_values['amplitude']
amperr = powout.params['amplitude'].stderr

if 'aniso' in filepath:
    powparstau2 = powmod.guess(tau2s,x=freqGHz)
    powout2 = powmod.fit(tau2s,powparstau2,x=freqGHz,weights=1/(tau2serr))

    print(fit_report(powout2.params))    
    fit2 = powout2.best_fit    
    alpha2 = -powout2.best_values['exponent']
    alphaerr2 = powout2.params['exponent'].stderr
    amp2 = powout2.best_values['amplitude']
    amperr2 = powout2.params['amplitude'].stderr
    
    taugeo = np.sqrt(taus*tau2s)
    taugeoerr = np.sqrt(tauserr*tau2serr)
    powparstauG = powmod.guess(taugeo,x=freqGHz)
    powoutG = powmod.fit(taugeo,powparstauG,x=freqGHz,weights=1/(taugeoerr))

    print(fit_report(powoutG.params))    
    fitG = powoutG.best_fit    
    alphaG = -powoutG.best_values['exponent']
    alphaerrG = powoutG.params['exponent'].stderr
    ampG = powoutG.best_values['amplitude']
    amperrG = powoutG.params['amplitude'].stderr
    
  
   
plt.close('all')

fig = plt.figure(figsize=(7,6))

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.title('PSR %s'%pulsar)
if np.shape(data)[1] >= 3:
    plt.errorbar(freqMHz,taus,yerr=tauserr,fmt='k*',markersize=9.0,capthick=2,linewidth=1.5,label=r'$\alpha = %.1f \pm %.1f$' %(alpha,alphaerr))
else:
    plt.plot(freqMHz,taus,'k*',markersize=9.0,linewidth=1.5,label=r'$\alpha = %.1f \pm %.1f$' %(alpha,alphaerr))
plt.plot(freqMHz,fit,'k--',linewidth=1.5)
if 'aniso' in filepath:
    plt.errorbar(freqMHz,tau2s,yerr=tau2serr,fmt='r^',markersize=9.0,capthick=2,linewidth=1.5,label=r'$\alpha_2 = %.1f \pm %.1f$' %(alpha2,alphaerr2))
    plt.plot(freqMHz,fit2,'r:',linewidth=1.5)
    plt.errorbar(freqMHz,taugeo,yerr=taugeoerr,fmt='co',markersize=9.0,capthick=2,linewidth=1.5,label=r'$\alpha_G = %.1f \pm %.1f$' %(alphaG,alphaerrG))
    plt.plot(freqMHz,fitG,'c-.',linewidth=1.5)
    
plt.xscale('log')
plt.yscale('log')
plt.yticks(fontsize=14)
ticksMHz = (freqMHz).astype(np.int)[0:len(freqMHz):2] 
plt.xticks(ticksMHz,ticksMHz,fontsize=14)
plt.legend(fontsize=14)
plt.xlabel(r'$\nu$ (MHz)',fontsize=16, labelpad=15.0)
plt.ylabel(r'$\tau$ (sec)',fontsize=16)
plt.xlim(xmin = 0.95*freqMHz[0],xmax=1.05*freqMHz[-1])
plt.gcf().subplots_adjust(bottom=0.15)




"""Create folder to save to"""
picpath = r'./SpectraPlots'
if not os.path.exists(picpath):
    os.makedirs(picpath)

Spectraplot = '%s_tauspectrum.png'  % pulsar
fileoutputtau = os.path.join(picpath,Spectraplot)
plt.savefig(fileoutputtau, dpi=150)
print 'Saved %s in %s' %(Spectraplot,picpath)


