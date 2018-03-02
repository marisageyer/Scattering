#!/usr/bin/python2.7

#Created on Mon Jan 25 13:22:33 2016

#@author: marisa

import argparse
import os, sys
import math
import pypsr_standalone as psr
import matplotlib.pyplot as plt
import numpy as np
from lmfit.models import LinearModel

import DataReadIn as dri

#from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
#rc('text', usetex=False)

"""Define options to the script"""
parser = argparse.ArgumentParser()
parser.add_argument('-f','--filename',
                    help="Provide the pathway to the data files")
parser.add_argument('-p','--period',type=float,
                    help="Provide the period of the pulsar in seconds")
parser.add_argument('-s','--simulate',
                    help="Choosing this option leads to simulated data. For broadening function: choose between 'onedim', 'iso', 'aniso'.")                    
parser.add_argument('-m','--method',
                    help="Choosing method to fit data or simulation. Choose between 'onedim' and 'iso'")                   
parser.add_argument('-dc','--datacycle',
                    help="The type of data. Used to label filenames.")
parser.add_argument('-t','--template',
                    help="filepath to txt file, containing a high frequency EPN profile. It will use the fixed with of this profile as the Guassian intrinsic template.")
args = parser.parse_args()

"""Allocate variable names to the parsed options"""
filepath = args.filename
pulseperiod = args.period
simu = args.simulate
meth = args.method
datac = args.datacycle
temp = args.template


if simu is None:       
    print "\n Reading in data \n"
    pulsars = []
#    pulsar, nch, nbins, lm_rms = dri.read_header(filepath)
    pulsar, nch, nbins,nsub, lm_rms, tsub = dri.read_headerfull(filepath)
    pulsars = [pulsar for i in range(nch)]
            
    print0 = "Pulsar name: %s" %pulsar
    print1 = "Number of channels: %d" %nch
    print2 = "Number of bins: %d" %nbins
    print3 = "RMS: %f" %lm_rms
    print4 = "Tsub: %f sec" %tsub 
    for k in range(4):
            print eval('print{0}'.format(k))
    
else:
    print "\n Simulating data \n"
    print(" \n 1. The %s broadening function" %simu)
    if simu in ('iso','onedim', 'aniso'):
        while True:
            try:
                taudivP = raw_input("Express max. tau as a fraction of the pulseperiod: ")
                taudivP = float(taudivP)
            except ValueError:
                print("Try again. Enter an int or float.")
                continue
            else:
                break
    else:
        while True:
            try:
                taudivP1 = raw_input("Express max. tau1 as a fraction of the pulseperiod: ")
                taudivP1 = float(taudivP1)
                taudivP2 = raw_input("Express max. tau2 as a fraction of the pulseperiod: ")
                taudivP2 = float(taudivP2)
                taudivP = np.array([taudivP1, taudivP2])
            except ValueError:
                print("Try again. Enter an int or float.")
                continue
            else:
                break
    print(" \n 2. The simulated pulsar")
    if pulseperiod is None:  
        while True:
            try:
                pulseperiod = raw_input("Express pulse period in seconds: ")
                pulseperiod = float(pulseperiod)
            except ValueError:
                print("Try again. Enter an int or float.")
                continue
            else:
                break  
    while True:
        try:
            dutycycle = raw_input("Express duty cycle as a percentage of pulse period, %s sec: " %pulseperiod)
            dutycycle = float(dutycycle)
        except ValueError:
            print("Try again. Enter an int or float.")
            continue
        else:
            break
    print(" \n 3. Data and observing properties")
    while True:
        try:
            freqlow, freqhigh, incr = raw_input("Choose a lowest and highest frequency and increment (MHz), e.g 50 120 5: ").split()
            freqlow, freqhigh, incr = float(freqlow), float(freqhigh), float(incr)
        except ValueError:
            print("Try again. Separate 3 floats with a space.")
            continue
        else:
            break

    while True:
        try:            
            snr = raw_input("Choose peak SNR: ")
            snr = int(snr)
        except ValueError:
            print("Try again. Enter an int or float.")
            continue
        else:
            break
    
    nbins = 512
    freqsimu = np.arange(freqlow,freqhigh+incr,incr)
    nch = len(freqsimu)
    propconst = taudivP*pulseperiod*(freqlow/1000.)**4
    propconst = np.array(propconst)
    
    tausecs = []
    for i in range(len(propconst)):
        tausec = propconst[i]*(freqsimu/1000.)**(-4)
        tausecs.append(tausec)
    tausecs = np.array(tausecs).transpose()
    
    pulsars = []
    for i in range(nch):
        if simu == 'aniso':
            pulsar = r'Simul.: $\tau_1 = %.2f$, $\tau_2 = %.2f$ ' %(tausecs[i][0],tausecs[i][1])
            pulsars.append(pulsar)
        else:
            pulsar = r'Simul: $\tau_1 = %.2f$' %tausecs[i]
            pulsars.append(pulsar)
    pulsar = 'Simulated'
    print0 = "Pulsar : %s" %pulsar
    print1 = "Number of channels: %d" %nch
    print2 = "Number of bins: %d" %nbins
    print3 = ""
    for k in range(4):
        print eval('print{0}'.format(k))


## Define time axis, and time/bins conversions

if pulseperiod == None: 
    print "Using Tsub in header to convert bins to time. Note Tsub here is full phase time, corresponding to nbins."
    pulseperiod = tsub

profilexaxis = np.linspace(0,pulseperiod,nbins)
pbs = pulseperiod/nbins
tbs = tsub/nbins

# Prepare to fit profile using tau_fitter

obtainedtaus = []
lmfittausstds = []
obtainedtaus2 = []
lmfittausstds2 = []

bestparamsall = []
bestparams_stdall = []
correls = []
redchis = []

freqmsMHz =[]
freqcsMHz =[]
noiselessmodels =[]
results = []
comp_rmss = []
comp_fluxes= []
comp_SNRs =[]
datas = []
climbvals =[]



halfway = nbins/2.

for i in range(nch):
    if simu is None:
            print "\n Channel %d" %i
	    # read in data 
            data, freqc, freqm = dri.read_data(filepath,i,nbins)
            freqmsMHz.append(freqm)
            freqcsMHz.append(freqc)
            # roll the data of lowest freq channel to middle of bins 
	    if i ==0:
                peakbin = np.argmax(data)
                shift = int(halfway -int(peakbin))
                print 'peak bin at lowest freq channel:%d' %peakbin
            else:
                peakbin = peakbin
                shift = int(halfway - int(peakbin)) 
            data = np.roll(data,shift)
            print "Rolling data by -%d bins" %shift
    else:
        freqmsMHz = freqsimu
        data = psr.simulate(pulseperiod,tausecs[i],dutycycle,-1.6,freqGHz[i],freqlow/1000.,nbins,snr,simu)
    comp_rms = psr.find_rms(data,nbins)
   
    if temp is None: 
        if meth is None:
            print "No fitting method was chosen. Will default to an isotropic fitting model. \n Use option -m with 'onedim' to change."
            result, noiselessmodel, besttau, taustd, bestparams, bestparams_std, redchi, corsig = psr.tau_fitter(data,nbins)
            climbval = psr.returnclimb(np.linspace(1,nbins,nbins),bestparams[1],bestparams[0],bestparams[2],besttau,bestparams[3],nbins)
        elif meth == 'iso':
            result, noiselessmodel, besttau, taustd, bestparams, bestparams_std, redchi, corsig = psr.tau_fitter(data,nbins)
            climbval = psr.returnclimb(np.linspace(1,nbins,nbins),bestparams[1],bestparams[0],bestparams[2],besttau,bestparams[3],nbins)
        elif meth == 'onedim':
            result, noiselessmodel, besttau, taustd, bestparams, bestparams_std, redchi, corsig = psr.tau_1D_fitter(data,nbins)
            climbval = psr.returnclimb1D(np.linspace(1,nbins,nbins),bestparams[1],bestparams[0],bestparams[2],besttau,bestparams[3],nbins)
        #elif meth in ('fix'):
        #     print "This uses fixed width values, as given by fixval array"
        #     fixvals = np.array([ 735.96503772,  570.80448298,  283.9382534 ,  220.03587608,
        #157.5615805 ,  142.65554165,  106.30809045,  114.00090536])
        #      noiselessmodel, besttau, taustd, bestparams, bestparams_std, redchi = psr.tau_fitter_fix1D(data,nbins,fixvals[i])
        else:
             print "Incorrect fitting method. Choose from iso or onedim"
    else:
        templatefile = np.loadtxt(temp)
        if meth is None:
            tempdata = templatefile[:,3]
            noiselessmodel, besttau, taustd, bestparams, bestparams_std, redchi = psr.tau_tempfitter(data,nbins,tempdata)
        elif meth == 'iso':
            tempdata = templatefile[:,3]
            noiselessmodel, besttau, taustd, bestparams, bestparams_std, redchi = psr.tau_tempfitter(data,nbins,tempdata)
        elif meth == 'onedim':
            tempdata = templatefile[:,3]
            noiselessmodel, besttau, taustd, bestparams, bestparams_std, redchi = psr.tau_1D_tempfitter(data,nbins,tempdata)                
        else:
           print "Incorrect fitting method. Choose from iso or onedim or multi (for analytical double peak)"
    


    """ INPUT SECTION ENDS """

            
    comp_flux = psr.find_modelflux(noiselessmodel,nbins)
    comp_SNR_model = psr.find_peaksnr(noiselessmodel,comp_rms)

    print 'Estimated SNR (from model peak and rms): %.2f' % comp_SNR_model
    comp_SNR =  psr.find_peaksnr_smooth(data,comp_rms)
    print 'Estimated SNR (from data peak and rms): %.2f' % comp_SNR

    print 'Channel Tau (ms): %.2f \pm %.2f ms' %(besttau,taustd)
    """Append all arrays - to have length of number of frequency channels"""
    
    obtainedtaus.append(besttau)    
    lmfittausstds.append(taustd)
    bestparamsall.append(bestparams)
    bestparams_stdall.append(bestparams_std)
    
    redchis.append(redchi)
    correls.append(corsig)          
    noiselessmodels.append(noiselessmodel)
    results.append(result)
    
    comp_rmss.append(comp_rms)
    comp_fluxes.append(comp_flux)
    comp_SNRs.append(comp_SNR)
    datas.append(data)
    climbvals.append(climbval)


"Pick out the correlations of sigma with Amp to use in flux errors"

cor_sigA = np.zeros(len(correls))
for i in range(len(correls)):
    if correls[i] is None:
        cor_sigA[i] = 0
    elif correls[i] is not None:
        cor_sigA[i] = correls[i]['A']
            
## insert a SNR cutoff here later if want to.
## have removed it for now

print "Using no SNR cutoff" 

data_highsnr =np.array(datas)
model_highsnr = np.array(noiselessmodels)

taus_highsnr = np.array(obtainedtaus)
lmfitstds_highsnr = np.array(lmfittausstds)

taussec_highsnr = taus_highsnr*pbs
lmfitstdssec_highsnr = lmfitstds_highsnr*pbs

climb_highsnr = np.array(climbvals)
    
freqMHz_highsnr = np.array(freqmsMHz)
fluxes_highsnr = np.array(comp_fluxes)
corsigA_highsnr = np.array(cor_sigA)
rms_highsnr = np.array(comp_rmss)
        
number_of_plotted_channels = len(data_highsnr)
npch = number_of_plotted_channels

"""Array with all the other fitting parameters: sigma, A, etc."""
bestpT = np.transpose(bestparamsall)
bestpT_std = np.transpose(bestparams_stdall)

print "Number of plotted channels: %d/%d" %(npch, nch)

bestpT_highSNR = bestpT
bestpT_std_highSNR = bestpT_std

"""Calculate fits for parameters sigma and mu"""

"""Fit a DM model to delta mu"""
delnuarray = [-(1/freqMHz_highsnr[-1]**2-1/freqMHz_highsnr[i]**2) for i in range(npch)] ##in MHz
delmuarray = [(bestpT_highSNR[1][-1] - bestpT_highSNR[1][i])*pbs for i in range(npch)] ##in seconds
delmu_stdarray = [(bestpT_std_highSNR[1][-1] - bestpT_std_highSNR[1][i])*pbs for i in range(npch)]
linmod = LinearModel()
DM_linpars = linmod.guess(delmuarray, x=delnuarray)
#	DM_linout  = linmod.fit(delmuarray, DM_linpars, x=delnuarray, weights=1/(np.power(delmu_stdarray,2)))
DM_linout  = linmod.fit(delmuarray, DM_linpars, x=delnuarray)

DM_CCval = DM_linout.best_values['slope']
DM_CCvalstd = DM_linout.params['slope'].stderr

DMmodelfit = DM_linout.best_fit ##model gives deltime in seconds (used to shift data)

DMconstant = 4148.808
#uncertainty in the constant is 0.003 - only affects the Delta DM value in the 9th decimal
DMval = (DM_CCval/DMconstant)
DMvalstd = (DM_CCvalstd/DMconstant)
DMcheck = psr.DM_checker(freqmsMHz,bestpT_highSNR[1]*pbs)

    
"""Plotting starts"""

plt.close('all')

#plot onedim in blue dashed
#else plot in red
if meth == 'onedim':
    prof = 'b--'
    lcol='b'
else:
    prof = 'r-'
    lcol ='r'

##PLOT PROFILES##  
    
numplots = int(np.ceil(npch/8.))

"""Compute residuals"""

resdata = data_highsnr - model_highsnr
resnormed = (resdata-resdata.mean())/resdata.std()

"""Plot 1: Pulse profiles and fits"""
if taussec_highsnr[0] > 1:
    taulabel =  taussec_highsnr
    taulabelerr = lmfitstdssec_highsnr
    taustring = 'sec'
else:
    taulabel = taussec_highsnr*1000
    taulabelerr = lmfitstdssec_highsnr*1000
    taustring = 'ms'

for k in range(numplots):
    j = 8*k
    figg = plt.figure(k+1,figsize=(14,8))
    for i in range(8):
    		figg.subplots_adjust(left = 0.08, right = 0.98, wspace=0.35,hspace=0.35,bottom=0.15)  
    		#plt.rc('text', usetex=True)
    		#plt.rc('font', family='serif')              
    		plt.subplot(2,4,i+1)
    		plt.plot(profilexaxis,data_highsnr[j+i],'k',alpha = 0.30)
    		plt.plot(profilexaxis,model_highsnr[j+i],prof,lw = 2.0, alpha
                        = 0.7,label=r'$\tau: %.2f \pm %.2f$ %s'
                        %(taulabel[j+i], taulabelerr[j+i], taustring))
    		plt.title('%s at %.1f MHz' %(pulsar, freqMHz_highsnr[j+i]))
    		plt.ylim(ymax=1.3*np.max(data_highsnr[j+i]))
    		plt.xlim(xmax=pulseperiod)
    		plt.xticks(fontsize=12)
    		plt.yticks(fontsize=12)
    		plt.xlabel('time (s)',fontsize=14)
    		plt.legend(fontsize=12,numpoints=1)
    		plt.ylabel('normalized intensity',fontsize=14)
plt.show()


for i in range(nch):
    print'Tau (ms): %.2f' %(1000*taussec_highsnr[i])

#print9 = 'Tsub = %.5f' %tsub
#print eval('print{0}'.format(9))

#"""Save txt files"""

txtpath = r'./FitTxtfiles'
if not os.path.exists(txtpath):
    os.makedirs(txtpath)
txtfile = '%s_%s_%s_FreqTau.txt' %(pulsar,datac,meth)
FPtxt = os.path.join(txtpath,txtfile) 


#1. Freq, Tau, Tauerr
headr = 'Pulsar: %s Nch: %d Nbins: %d \n Freq (MHz), Tau (sec), Tau Err (sec)' %(pulsar,nch,nbins)
np.savetxt(FPtxt,np.column_stack((freqMHz_highsnr,taussec_highsnr,lmfitstdssec_highsnr)),header=headr)

print "%s saved in %s" %(txtfile,txtpath)


#"""Create folder to save plots to"""
#picpath = r'./Output_Plots'
#if not os.path.exists(picpath):
#    os.makedirs(picpath)
# 
#for i in range(numplots):
#    Summaryplot = '%s_%s_%s_%.2d.png'  % (pulsar,datac,meth, nch-i)
#    fileoutputtau = os.path.join(picpath,Summaryplot)
#    plt.savefig(fileoutputtau, dpi=150)
#    print 'Saved %s in %s' %(Summaryplot,picpath)
#    plt.close()
