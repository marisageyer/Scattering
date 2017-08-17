# -*- coding: utf-8 -*-
"""
Created on Tue Oct 28 09:49:12 2014

@author: Marisa
"""

import numpy as np
import matplotlib.pyplot as plt
from math import exp,sqrt,pi,erf,log
from scipy.special import iv, ive
from lmfit import Model, conf_interval, printfuncs
from lmfit.models import GaussianModel
#from sympy.functions.special.delta_functions import Heaviside
#import GPy
#import os

#plt.close()

def makeprofile(nbins = 2**9, ncomps = 1, amps = 1, means = 100, sigmas = 10):
    if ncomps == 1:
        npamps = np.array([amps])
        npmeans = np.array([means])
        npsigmas = np.array([sigmas])
    else:    
        npamps = np.array(amps)
        npmeans = np.array(means)
        npsigmas = np.array(sigmas)
   
    profile = np.zeros(nbins)
    x = np.linspace(1,nbins,nbins)
    
    for i in range(ncomps):
#        print npmeans[i]
        profile = profile + \
        npamps[i]*np.exp(-pow(x-npmeans[i],2)/(2*pow(npsigmas[i],2)))
    return x, profile            

def pulsetrain(npulses = 10, bins = np.linspace(1,512,512), profile = np.zeros(512)):
    nbins = np.max(bins)
    train = np.zeros(npulses*int(nbins))
    nbins = int(nbins)
    for i in range(npulses):
        startbin = int(i*nbins)
        endbin = int(startbin + nbins)
        train[startbin:endbin] = profile
    return train

def pulsetrain_bins(npulses, numberofbins, profile):
    binsrange = np.linspace(1,numberofbins,numberofbins)    
    nbins = np.max(binsrange)
    nbins = int(nbins)
#    print nbins
    train = np.zeros(npulses*int(nbins))

    for i in range(npulses):
        startbin = int(i*nbins)
        endbin = int(startbin + nbins)
        train[startbin:endbin] = profile
    return train

    
#def psrscatter(brfunc, profile):    
#    scattered = np.convolve(profile,brfunc)
#    profint = np.sum(profile)    
#    scint = np.sum(scattered)
#    scattered = scattered / scint * profint
#    bins = profile.shape[0]
#    out = scattered[0:bins]    
#    return out
    
    
def psrscatter(brfunc, profile):    
    scattered = np.convolve(profile,brfunc)
    profint = np.sum(profile)    
    scint = np.sum(scattered)
    scatterednorm = scattered / scint * profint
    bins = profile.shape[0]
    out = scatterednorm[0:bins]    
    return out

def psrscatter_noconserve(brfunc, profile):    
    scattered = np.convolve(profile,brfunc)
#    profint = np.sum(profile)    
#    scint = np.sum(scattered)
#    scattered = scattered / scint * profint
    bins = profile.shape[0]
    out = scattered[0:bins]    
    return out

def psrscatter_noconservefold(brfunc, profile,trainlength):    
    scattered = np.convolve(profile,brfunc)
#    profint = np.sum(profile)    
#    scint = np.sum(scattered)
#    scattered = scattered / scint * profint
    bins = profile.shape[0]
    out = scattered[0:bins*trainlength]    
    return out


def psrscatterpostfold(brfunc, profile):    
    scattered = np.convolve(profile,brfunc)
    profint = np.sum(profile)    
    scint = np.sum(scattered)
    scattered = scattered / scint * profint
    length = scattered.shape[0]
    out = scattered[0:length]    
    return out



#def extractpulse(train, pulsesfromend, binsperpulse):
#    start = -pulsesfromend*binsperpulse - 1
#    end = start + binsperpulse    
#    minbin = np.argmin(train[start:end])
#    start2 = start + minbin
#    end2 = start2 + binsperpulse
#    zerobpulse = train[start2:end2]-np.min(train[start2:end2])
#    rectangle = np.min(train[start2:end2])*binsperpulse
#    flux = np.sum(train[start2:end2]) - rectangle
#    return train[start2:end2], zerobpulse, flux

def profilespec(nu,specindex,profile):
    profilefreqspec = []
    for i in range(len(nu)):
        profilefreq = pow(nu[0],specindex)*pow(nu[i],-specindex)*profile
        profilefreqspec.append(profilefreq)
    return profilefreqspec   
    
def sigmaa(kappa,nu):                                       #nu is observing frequency in GHz
        scatteringstrength = kappa*nu**-2                   #kappa is the proportionality constant (radians) and the scattering strength at 1 GHz
        return scatteringstrength    

def tau(D,Ds,kappa,nu,light):                                     #D distance from observer to source in kpc, Ds distance from source to screen in kpc               
    scatteringtime = (D-Ds)*Ds*(sigmaa(kappa,nu))**2/(light*D)    #Tau in seconds
    return scatteringtime

def step(x):
    return 1 * (x >= 0)


#Create scattered profiles according to different broadening functions respectively and extract a single pulse

### Broadening functions
# 1. Isotropic scattering, isotropic screen

def broadfunc(x,tau):
    broadfunc = (1/tau)*np.exp(-x/tau)
    return broadfunc

def broadfuncstep(x,tau):
    tau = float(tau)
    broadfunc = (1/tau)*np.exp(-x/tau)*step(x)
    return broadfunc   

def broadfuncExt(x,tauE):
    p1 = np.sqrt(np.pi*tauE/(2*x**3))
    p22 = []
    for n in range(1,15,2):    
        p2a = (n**2*np.pi**2*tauE/(2*x))-1
        p2b = np.exp(-n**2*np.pi**2*tauE/(4*x))
        p2 = p2a*p2b
        p22.append(p2)
    p22sum = np.sum(p22,axis=0)
    return p1*p22sum
    
def tauE(D,kappa,nu,light):  ##tau for extended medium
    timeconst = 3*D*(sigmaa(kappa,nu))**2/(np.pi**2*light)
    return timeconst

def foldedexpclimb(nexp,period,tau):  #enter in bins (not seconds)
    expfunc2 = []
    tau = float(tau)
    period = float(period)
    periodspace = np.linspace((nexp-1)*period,nexp*period,period)   
    for i in range(nexp):
       expfunc = (1/tau)*np.exp(-(periodspace-i*period)/tau)
       expfunc2.append(expfunc)
       expfunc3 = np.array(expfunc2) 
       expadd = expfunc3.sum(axis=0)
    return expadd 

def foldedexp(nexp,period,tau):  #enter in bins (not seconds)
    expfunc2 = []
    tau = float(tau)
    period = float(period)
    periodspace = np.linspace((nexp-1)*period,nexp*period,period)   
    for i in range(nexp):
       expfunc = (1/tau)*np.exp(-(periodspace-i*period)/tau)
       expfunc2.append(expfunc)
       expfunc3 = np.array(expfunc2) 
       expadd = expfunc3.sum(axis=0)
       expadd = expadd - np.min(expadd)
    return expadd 

# 2. Anisotropic scattering, isotropic screen - Bessell function result
def broadfunc2(x,tau1,tau2):
    broadfunc2 = 1/(np.sqrt(tau1*tau2))*iv(0,x*(1/(2*tau2)-1/(2*tau1)))*np.exp(-x*(1/(2*tau1) + 1/(2*tau2)))*step(x) ##step is new 18 May 2017
    return broadfunc2  

def broadfunc2part(x,tau1,tau2):
    if tau1 != tau2:
        broadfunc2 = 1/(np.sqrt(tau1*tau2))*ive(0,x*(1/(2*tau2)-1/(2*tau1)))
    else:
        broadfunc2 = 1/(np.sqrt(tau1*tau2))*iv(0,x*(1/(2*tau2)-1/(2*tau1)))*np.exp(-x*(1/(2*tau1) + 1/(2*tau2)))
    return broadfunc2
    
def broadfunc2ive(x,tau1,tau2):
    broadfunc2ive = 1/(np.sqrt(tau1*tau2))*ive(0,x*(1/(2*tau2)-1/(2*tau1)))
    return broadfunc2ive

def broadfunc1D(x,tau):
    broadfunc1 = (1.0/np.sqrt(x*tau*pi))*np.exp(-x/tau)
    return broadfunc1 

def broadfunctrunc(x,tau,xmax):
    ##x variable is the time variable
    tau = float(tau)
    broadfunc = (1/tau)*np.exp(-x/tau)*step(xmax-x)
    return broadfunc

 
def foldedexp2(nexp,period,tau1,tau2):  #enter in bins (not seconds)
    expfunc2 = []
    tau1 = float(tau1)
    tau2 = float(tau2)
    period = float(period)
    periodspace = np.linspace((nexp-1)*period,nexp*period,period)   
    for i in range(nexp):
       expfunc =  broadfunc2(periodspace-i*period,tau1,tau2)
       expfunc2.append(expfunc)
       expfunc3 = np.array(expfunc2) 
       expadd = expfunc3.sum(axis=0)
       expadd = expadd - np.min(expadd)
    return expadd

def extractpulse(train, pulsesfromend, binsperpulse):
    if pulsesfromend == 0:
        start = 0
        end = binsperpulse
        zerobpulse = train[start:end]-np.min(train[start:end])
        rectangle = np.min(train[start:end])*binsperpulse
        flux = np.sum(train[start:end]) - rectangle
        return train[start:end], zerobpulse, rectangle, flux
    
    else:     
        start = -pulsesfromend*binsperpulse
        end = start + binsperpulse 
        zerobpulse = train[start:end]-np.min(train[start:end])
        rectangle = np.min(train[start:end])*binsperpulse
        flux = np.sum(train[start:end]) - rectangle
        return train[start:end], zerobpulse, rectangle, flux
      
def peaksnr(x, profile, snr):
    bins = profile.shape[0]
    noise = np.random.random(bins)
    peak = np.max(profile)
    out = profile/peak * snr + noise
    return out
    

def GxETrain(x,mu,sigma, A, tau, dc, nbins):
#This model convolves a pulsetrain with a broadening function
#It extracts one of the last convolved profiles, subtracts the climbed baseline and then adds noise to it
    mu, sigma, A, tau = float(mu),float(sigma), float(A), float(tau)    
    bins, profile = makeprofile(nbins = nbins, ncomps = 1, amps = A, means = mu, sigmas = sigma)
    binstau = np.linspace(1,nbins,nbins)
#    binstau = x
    scat = psrscatter_noconserve(broadfuncstep(binstau,tau),pulsetrain_bins(3, nbins, profile))  
#    scat = psrscatter(broadfuncstep(binstau,tau),pulsetrain_bins(3, nbins, profile))  
    climb, observed_nonoise, rec, flux = extractpulse(scat, 2, nbins)
#    print np.min(climb)
    return observed_nonoise + dc
    
def returnclimb(x,mu,sigma,A,tau,dc,nbins):
    mu, sigma, A, tau = float(mu),float(sigma), float(A), float(tau)    
    bins, profile = makeprofile(nbins = nbins, ncomps = 1, amps = A, means = mu, sigmas = sigma)
    binstau = np.linspace(1,nbins,nbins)
    scat = psrscatter_noconserve(broadfuncstep(binstau,tau),pulsetrain_bins(3, nbins, profile))  
    climb, observed_nonoise, rec, flux = extractpulse(scat, 2, nbins)
    return np.min(climb)
    


def GxETrainClimb(x,mu,sigma, A, tau, dc, nbins):
#This model convolves a pulsetrain with a broadening function
#It extracts one of the last convolved profiles, subtracts the climbed baseline and then adds noise to it
    mu, sigma, A, tau = float(mu),float(sigma), float(A), float(tau)    
    bins, profile = makeprofile(nbins = nbins, ncomps = 1, amps = A, means = mu, sigmas = sigma)
    binstau = np.linspace(1,nbins,nbins)
#    binstau = x
    scat = psrscatter_noconserve(broadfuncstep(binstau,tau),pulsetrain_bins(3, nbins, profile))  
#    scat = psrscatter(broadfuncstep(binstau,tau),pulsetrain_bins(3, nbins, profile))  
    climb, observed_nonoise, rec, flux = extractpulse(scat, 2, nbins)
    return climb + dc


def GxETrainTemp(x, tau, dc, nbins, A, templatedata):
#This model fits a Guassian function to a high freq template and then
#convolves a pulsetrain made up of these templates with a broadening function
#It extracts one of the last convolved profiles, subtracts the climbed baseline and then adds noise to it     
    binstau = np.linspace(1,nbins,nbins)    
    Gmod = GaussianModel()
    gpars = Gmod.guess(templatedata, x=binstau)
    gout = Gmod.fit(templatedata,gpars, x=binstau)
#    fitA = gout.best_values['amplitude']
#    fitS = gout.best_values['sigma']
#    fitC = gout.best_values['center']
    profile = A*gout.best_fit    
#    bins, profile = A*makeprofile(nbins = nbins, ncomps = 1, amps = fitA, means = fitC, sigmas = fitS)
    pulsetrain = pulsetrain_bins(3, nbins, profile)
    scat = psrscatter(broadfunc(binstau,tau),pulsetrain)   
    climb, observed_nonoise, rec, flux = extractpulse(scat, 2, nbins)
    return observed_nonoise + dc

    

def GxETrainMulti(x, ncomps, mus1, mus2, sigmas1, sigmas2, As1, As2, tau, dc, nbins):
#    As = np.array(As)
#    mus = np.array(mus)
#    sigmas = np.array(sigmas)
#    As1,As2,mus1,mus2,sigmas1,sigmas2 = As[0],As[1],mus[0],mus[1],sigmas[0],sigmas[1]    
    bins, profile = makeprofile(nbins = nbins, ncomps = 2, amps = [As1,As2], means = [mus1,mus2], sigmas = [sigmas1,sigmas2])
    binstau = np.linspace(1,nbins,nbins)
#    binstau = x
    scat = psrscatter(broadfunc(binstau,tau),pulsetrain_bins(3, nbins, profile))   
    climb, observed_nonoise, rec, flux = extractpulse(scat, 2, nbins)
    return observed_nonoise + dc

def GxETrainMulti1D(x, ncomps, mus1, mus2, sigmas1, sigmas2, As1, As2, tau, dc, nbins):
#    As = np.array(As)
#    mus = np.array(mus)
#    sigmas = np.array(sigmas)
#    As1,As2,mus1,mus2,sigmas1,sigmas2 = As[0],As[1],mus[0],mus[1],sigmas[0],sigmas[1]    
    bins, profile = makeprofile(nbins = nbins, ncomps = 2, amps = [As1,As2], means = [mus1,mus2], sigmas = [sigmas1,sigmas2])
    binstau = np.linspace(1,nbins,nbins)
#    binstau = x
    scat = psrscatter(broadfunc1D(binstau,tau),pulsetrain_bins(3, nbins, profile))   
    climb, observed_nonoise, rec, flux = extractpulse(scat, 2, nbins)
    return observed_nonoise + dc



def GxETrainLong(x,mu,sigma, A, tau, trainlength, dc, nbins):
#This model convolves a pulsetrain with a broadening function
#It extracts one of the last convolved profiles, subtracts the climbed baseline and then adds noise to it
    bins, profile = makeprofile(nbins = nbins, ncomps = 1, amps = A, means = mu, sigmas = sigma)
    binstau = np.linspace(1,trainlength*nbins,trainlength*nbins)
    scat = psrscatter(broadfunc(binstau,tau),pulsetrain_bins(20, nbins, profile))   
    climb, observed_nonoise, rec, flux = extractpulse(scat, 2, nbins)
    return observed_nonoise + dc


def GxETrainAniDC(x,mu, sigma, A, tau1, tau2, dc, nbins):
    bins, profile = makeprofile(nbins = nbins, ncomps = 1, amps = A, means = mu, sigmas = sigma)  
    binstau = np.linspace(1,nbins,nbins)
    scat = psrscatter(broadfunc2(binstau,tau1,tau2),pulsetrain(3, bins, profile))   
    climb, observed_nonoise, rec,flux = extractpulse(scat, 2, nbins)
    return observed_nonoise + dc

def GxETrain1D(x,mu, sigma, A, tau1, dc, nbins):
    mu, sigma, A, tau1 = float(mu),float(sigma), float(A), float(tau1)
    bins, profile = makeprofile(nbins = nbins, ncomps = 1, amps = A, means = mu, sigmas = sigma)  
    binstau = np.linspace(1,nbins,nbins)
    scat = psrscatter(broadfunc1D(binstau,tau1),pulsetrain(3, bins, profile))   
    climb, observed_nonoise, rec,flux = extractpulse(scat, 2, nbins)
    return observed_nonoise + dc

def returnclimb1D(x,mu,sigma,A,tau1,dc,nbins):
    mu, sigma, A, tau1 = float(mu),float(sigma), float(A), float(tau1)
    bins, profile = makeprofile(nbins = nbins, ncomps = 1, amps = A, means = mu, sigmas = sigma)  
    binstau = np.linspace(1,nbins,nbins)
    scat = psrscatter(broadfunc1D(binstau,tau1),pulsetrain(3, bins, profile))   
    climb, observed_nonoise, rec,flux = extractpulse(scat, 2, nbins)
    return np.min(climb)
    

def GxESingleFoldDC(x,mu,sigma,A,tau,trainlength,dc,nbins):
    #This model takes a single Guassian pulse with mean mu and sigma
    #Convolves it with a broadening function
    #It extracts one of the last convolved profiles subtracts the climbed baseline and then adds noise to it       
    observed_postfold = np.zeros(nbins)      
    bins, profile = makeprofile(nbins = nbins, ncomps = 1, amps = A, means = mu, sigmas = sigma)
    binstau = np.linspace(1,trainlength*nbins,trainlength*nbins)
    scat = psrscatterpostfold(broadfunc(binstau,tau),pulsetrain(1, bins, profile))
    climb, observed_nonoise, rec, flux = extractpulse(scat, 0, trainlength*nbins)
    for i in range(trainlength*nbins):    
        observed_postfold[np.mod(i,nbins)] += observed_nonoise[i]         
        GxESingleFold = observed_postfold[x]-np.min(observed_postfold[0:nbins])
    return GxESingleFold +dc

def GxESingleFoldnoDC(x,mu,sigma,A,tau,trainlength,nbins):
    #This model takes a single Guassian pulse with mean mu and sigma
    #Convolves it with a broadening function
    #It extracts one of the last convolved profiles subtracts the climbed baseline and then adds noise to it       
    observed_postfold = np.zeros(nbins)      
    bins, profile = makeprofile(nbins = nbins, ncomps = 1, amps = A, means = mu, sigmas = sigma)
    binstau = np.linspace(1,trainlength*nbins,trainlength*nbins)
    scat = psrscatterpostfold(broadfunc(binstau,tau),pulsetrain(1, bins, profile))
    climb, observed_nonoise, rec, flux = extractpulse(scat, 0, trainlength*nbins)
    for i in range(trainlength*nbins):    
        observed_postfold[np.mod(i,nbins)] += observed_nonoise[i]         
        GxESingleFold = observed_postfold[x]-np.min(observed_postfold[0:nbins])
    return GxESingleFold


def GxETrainAniDC2part(x,mu, sigma, A, tau1, tau2, dc, nbins):
    bins, profile = makeprofile(nbins = nbins, ncomps = 1, amps = A, means = mu, sigmas = sigma)  
    binstau = np.linspace(1,nbins,nbins)
    scat = psrscatter(broadfunc2part(binstau,tau1,tau2),pulsetrain(3, bins, profile))   
    climb, observed_nonoise, rec,flux = extractpulse(scat, 2, nbins)
    return observed_nonoise + dc

def GxETrainAniDCive(x,mu, sigma, A, tau1, tau2, dc, nbins):
    bins, profile = makeprofile(nbins = nbins, ncomps = 1, amps = A, means = mu, sigmas = sigma)  
    binstau = np.linspace(1,nbins,nbins)
    scat = psrscatter(broadfunc2ive(binstau,tau1,tau2),pulsetrain(3, bins, profile))   
    climb, observed_nonoise, rec,flux = extractpulse(scat, 2, nbins)
    return observed_nonoise + dc



def smooth12(y, box_pts):
    gauss = np.ones(box_pts)
    #    box = np.ones(box_pts)/box_pts
    sigma = (1./12.)*box_pts
    mean = box_pts/2.
    for i in range(box_pts):
        gauss[i] = (1/(sigma*np.sqrt(2*np.pi)))*np.exp(-(i-mean)**2/(2*sigma**2))
    y_smooth = np.convolve(y, gauss, mode='same')
    return y_smooth


def smooth(y, box_pts):
    gauss = np.ones(box_pts)    
#    box = np.ones(box_pts)/box_pts
    sigma = (1./6.)*box_pts
    mean = box_pts/2.
    for i in range(box_pts):    
        gauss[i] = (1/(sigma*np.sqrt(2*np.pi)))*np.exp(-(i-mean)**2/(2*sigma**2))
    y_smooth = np.convolve(y, gauss, mode='same')
    return y_smooth


    
def PLM(x,K,k):
    return K*pow(x,-k)

def TemplateBasedModel(x,nbins,tau,A1,A2,mu1,mu2,sig1,sig2,bdc,AA,DC):
    ## These: bP1, bP1_std, bP2, bP2_std,bdc are results from FitTemplate_TwoPeaks
    ## written out here as A1,A2,mu1,mu2,sig1,sig2,bdc
    bP1 = np.array([sig1,mu1,A1])
    bP2 = np.array([sig2,mu2,A2])
    xbr = np.linspace(1,3*nbins,3*nbins)
    profile = TwoPeaksModel(1,nbins,bP1[2],bP2[2],bP1[1],bP2[1],bP1[0],bP2[0],bdc)
    scattemp = AA*psrscatter(broadfunc(xbr,tau),profile) + DC
    scatbins =  scattemp[0:nbins]
    return scatbins
   
   
def tau_fitdatatotemp(data,template,bP1,bP2,bP1_std,bP2_std,bdc,bdc_std):
    normdata = data/np.max(data)
    max_data, max_temp = np.argmax(normdata), np.argmax(template)
    shift = max_temp-max_data
    rollnorm = np.roll(normdata, shift)
    A1, A2, A1_std, A2_std = bP1[2], bP2[2], bP1_std[2], bP2_std[2]
    mu1, mu2, mu1_std, mu2_std = bP1[1], bP2[1],bP1_std[1], bP2_std[1]
    sig1, sig2,sig1_std, sig2_std = bP1[0], bP2[0],bP1_std[0], bP2_std[0]
    ## These: bP1, bP1_std, bP2, bP2_std,bdc are results from FitTemplate_TwoPeaks
    ## written out here as A1,A2,mu1,mu2,sig1,sig2,bdc
    nbins = len(data)
    model = Model(TemplateBasedModel)
    k = 5.0
    model.set_param_hint('nbins', value=nbins, vary=False)
    model.set_param_hint('A1', value=A1, vary=True, min = A1-k*A1_std, max=A1+k*A1_std)
    model.set_param_hint('A2', value=A2, vary=True,min = A2-k*A2_std, max=A2+k*A2_std)
    model.set_param_hint('sig1', value=sig1, vary=True, min = sig1-k*sig1_std, max=sig1+k*sig1_std)
    model.set_param_hint('sig2', value=sig2, vary=True, min = sig2-k*sig2_std, max=sig2+k*sig2_std)
    model.set_param_hint('mu1', value=mu1, vary=True, min = mu1-k*mu1_std, max=mu1+k*mu1_std)
    model.set_param_hint('mu2', value=mu2, vary=True, min = mu2-k*mu2_std, max=mu2+k*mu2_std)
    model.set_param_hint('bdc', value=bdc, vary=True)
    model.set_param_hint('DC', value=0.1, vary=True, min =0, max=np.max(data)/10.)
    model.set_param_hint('AA', value=1.2, vary=True)  
    model.set_param_hint('tau', value=3.0, vary=True,min=0,max=nbins)
    
    pars = model.make_params()
    xax=np.linspace(1,nbins,nbins)
    #"""Fit data"""
    
    result = model.fit(rollnorm,pars,x=xax)
    print(result.fit_report(show_correl = False))
    modfit = result.best_fit
    return rollnorm, modfit
     

   
#    for i in range(trainlength*nbins):    
#        observed_postfold[np.mod(i,nbins)] += observed_nonoise[i]         
#        GxESingleFold = observed_postfold[x]-np.min(observed_postfold[0:nbins])
#    return GxESingleFold +dc

def TwoPeaksModel(x,nbins,A1,A2,mu1,mu2,sig1,sig2,dc):
    bins, twoprofile = makeprofile(nbins = nbins, ncomps = 2, amps = [A1,A2], means = [mu1,mu2], sigmas = [sig1,sig2])
    return twoprofile + dc
    
    
    
def FitTemplate_TwoPeaks(templatedata):
    nbins = 1024
    model = Model(TwoPeaksModel)
    model.set_param_hint('nbins', value=nbins, vary=False)       
    model.set_param_hint('sig1', value=20, vary=True, min=0)
    model.set_param_hint('sig2', value=8, vary=True,min=0)
    model.set_param_hint('mu1', value=700, vary=True,min =500, max = nbins)
    model.set_param_hint('mu2', value=730, vary=True,min =500, max = nbins)
    model.set_param_hint('A1',value=np.max(templatedata), vary=True)
    model.set_param_hint('A2',value=np.max(templatedata)/2., vary=True)
    model.set_param_hint('dc',value=0.01,vary=True)
    
    pars = model.make_params()
    xax=np.linspace(1,nbins,nbins)
    #"""Fit data"""
    result = model.fit(templatedata,pars,x=xax)
    print(result.fit_report(show_correl = False))
    modfit = result.best_fit
    bsig1 = result.best_values['sig1']
    bsig2 = result.best_values['sig2']
    bmu1 = result.best_values['mu1']
    bmu2 = result.best_values['mu2']
    bA1= result.best_values['A1']
    bA2 = result.best_values['A2']
    bdc = result.best_values['dc']
    
    bsig1_std = result.params['sig1'].stderr
    bsig2_std = result.params['sig2'].stderr
    bmu1_std = result.params['mu1'].stderr
    bmu2_std = result.params['mu2'].stderr
    bA1_std = result.params['A1'].stderr
    bA2_std = result.params['A2'].stderr
    bdc_std = result.params['dc'].stderr

    
    bP1 = np.array([bsig1,bmu1,bA1])
    bP1_std = np.array([bsig1_std,bmu1_std,bA1_std])
    bP2 = np.array([bsig2,bmu2,bA2])
    bP2_std = np.array([bsig2_std,bmu2_std,bA2_std])
  

    return modfit, bP1, bP1_std, bP2, bP2_std,bdc,bdc_std
     

def FitTemplate(templatedata,nbins):
    binstau = np.linspace(1,len(templatedata),len(templatedata))   
    Gmod = GaussianModel()
    gpars = Gmod.guess(templatedata, x=binstau)
    gout = Gmod.fit(templatedata,gpars, x=binstau)
    fitA = gout.best_values['amplitude']
    fitS = gout.best_values['sigma']
#    fitC = gout.best_values['center']
    bf = gout.best_fit  
    scaleS = fitS*nbins/len(templatedata)
    return bf, fitA, scaleS



def tau_tempfitter(data,nbins,templatedata):
    ##This template fitter uses the Gaussian fit from FitTemplate
    ## From FitTemplate it extracts a FIXED WIDTH
    ## I.e this fitter is used when I want to compare my fits to other values in the literature
    ## where I know authors have used a high frequency EPN profile Gaussian fit with a fixed width
    bf, A, S = FitTemplate(templatedata,nbins)
    binpeak = np.argmax(data)    
    profile_peak = np.max(data)
#    A_init = profile_peak*S*np.sqrt(2*np.pi)
    modelname = GxETrain
    model = Model(modelname)

    model.set_param_hint('nbins', value=nbins, vary=False)       
#    model.set_param_hint('sigma', value=S, vary=True,min =0, max = nbins)
    model.set_param_hint('sigma', value=S, vary=False)
    model.set_param_hint('mu', value=binpeak, vary=True,min =0, max = nbins)
    model.set_param_hint('A',value=profile_peak, vary=True)
    model.set_param_hint('tau',value=200, vary=True, min=0)
    model.set_param_hint('dc',value = 0, vary = True)
    pars = model.make_params()
    xax=np.linspace(1,nbins,nbins)
    #"""Fit data"""
    result = model.fit(data,pars,x=xax)
    print(result.fit_report(show_correl = False))
    
    noiselessmodel = result.best_fit
    besttau = result.best_values['tau']
    taustd = result.params['tau'].stderr  ##estimated 1 sigma error

    bestsig = result.best_values['sigma']
    bestmu = result.best_values['mu']
    bestA = result.best_values['A']
    bestdc = result.best_values['dc']
    
    bestsig_std = result.params['sigma'].stderr
    bestmu_std = result.params['mu'].stderr
    bestA_std = result.params['A'].stderr
    bestdc_std = result.params['dc'].stderr    
    
    bestparams = np.array([bestsig,bestmu,bestA,bestdc])
    bestparams_std = np.array([bestsig_std,bestmu_std,bestA_std,bestdc_std])
    
    rchi = result.redchi    
    #return best values and std errors on the other parameters as well    
    
    return noiselessmodel, besttau, taustd, bestparams, bestparams_std, rchi


def tau_fitterClimb(data,nbins):
    ##for PSR 1848 and PSR 1914-0440 use 0.805*binpeak
    profile_peak = np.max(data)
    binpeak = np.argmax(data)  
    modelname = GxETrainClimb
    model = Model(modelname)
            
    model.set_param_hint('nbins', value=nbins, vary=False)       
    model.set_param_hint('sigma', value=15, vary=True, min =0, max = nbins)
    model.set_param_hint('mu', value=0.8*binpeak, vary=True, min=0, max = nbins)
    model.set_param_hint('A',value=profile_peak, vary=True, min=0)
    model.set_param_hint('tau',value=200, vary=True, min=0)
    model.set_param_hint('dc',value = 0, vary = True)
    pars = model.make_params()
    xax=np.linspace(1,nbins,nbins)
    #"""Fit data"""
    result = model.fit(data,pars,x=xax)
    print(result.fit_report(show_correl = False))
    
    noiselessmodel = result.best_fit
    besttau = result.best_values['tau']
    taustd = result.params['tau'].stderr  ##estimated 1 sigma error

    bestsig = result.best_values['sigma']
    bestmu = result.best_values['mu']
    bestA = result.best_values['A']
    bestdc = result.best_values['dc']
    
    bestsig_std = result.params['sigma'].stderr
    bestmu_std = result.params['mu'].stderr
    bestA_std = result.params['A'].stderr
    bestdc_std = result.params['dc'].stderr    
    
    bestparams = np.array([bestsig,bestmu,bestA,bestdc])
    bestparams_std = np.array([bestsig_std,bestmu_std,bestA_std,bestdc_std])
    
    rchi = result.redchi    
    
    #return best values and std errors on the other parameters as well    
    
    return result, noiselessmodel, besttau, taustd, bestparams, bestparams_std, rchi




def tau_fitter(data,nbins):
    ##for PSR 1848 and PSR 1914-0440 use 0.805*binpeak
    profile_peak = np.max(data)
    binpeak = np.argmax(data)  
    modelname = GxETrain
    model = Model(modelname)
            
    model.set_param_hint('nbins', value=nbins, vary=False)       
    model.set_param_hint('sigma', value=15, vary=True, min =0, max = nbins)
    model.set_param_hint('mu', value=0.805*binpeak, vary=True, min=0, max = nbins)
    model.set_param_hint('A',value=profile_peak, vary=True, min=0)
    model.set_param_hint('tau',value=200, vary=True, min=0)
    model.set_param_hint('dc',value = 0, vary = True)
    pars = model.make_params()
    xax=np.linspace(1,nbins,nbins)
    #"""Fit data"""
    result = model.fit(data,pars,x=xax)
#    print(result.fit_report(show_correl = True))
    
    noiselessmodel = result.best_fit
    besttau = result.best_values['tau']
    taustd = result.params['tau'].stderr  ##estimated 1 sigma error

    bestsig = result.best_values['sigma']
    bestmu = result.best_values['mu']
    bestA = result.best_values['A']
    bestdc = result.best_values['dc']
    
    bestsig_std = result.params['sigma'].stderr
    bestmu_std = result.params['mu'].stderr
    bestA_std = result.params['A'].stderr
    bestdc_std = result.params['dc'].stderr    
    
    bestparams = np.array([bestsig,bestmu,bestA,bestdc])
    bestparams_std = np.array([bestsig_std,bestmu_std,bestA_std,bestdc_std])
    
    """correlations with sigma"""    
    corsig = result.params['sigma'].correl
    #corA = result.params['A'].correl
    #corlist = [corsig,corA]
    
    
    rchi = result.redchi
    #return best values and std errors on the other parameters as well    
    
    return result, noiselessmodel, besttau, taustd, bestparams, bestparams_std, rchi, corsig


def tau_fitter_fix(data,nbins,fixval):
    #Fixed width fitter, need to pass it the fixed values
    #Also used to fix tau
    
    profile_peak = np.max(data)
    binpeak = np.argmax(data)  
    modelname = GxETrain
    model = Model(modelname)
            
    model.set_param_hint('nbins', value=nbins, vary=False)       
    model.set_param_hint('sigma', value=fixval, vary=False)
#    model.set_param_hint('sigma', value=15, vary=True, min =0, max = nbins)
    model.set_param_hint('mu', value=binpeak, vary=True, min=0, max = nbins)
    model.set_param_hint('A',value=profile_peak, vary=True)
    model.set_param_hint('tau',value=200, vary=True, min=0)
#    model.set_param_hint('tau',value=fixval, vary=False)
#    model.set_param_hint('dc',value = 0, vary = False)
    model.set_param_hint('dc',value = 0, vary =True)
    pars = model.make_params()
    xax=np.linspace(1,nbins,nbins)
    #"""Fit data"""
    result = model.fit(data,pars,x=xax)
    print(result.fit_report(show_correl = False))
    
    noiselessmodel = result.best_fit
    besttau = result.best_values['tau']
    taustd = result.params['tau'].stderr  ##estimated 1 sigma error

    bestsig = result.best_values['sigma']
    bestmu = result.best_values['mu']
    bestA = result.best_values['A']
    bestdc = result.best_values['dc']
    
    bestsig_std = result.params['sigma'].stderr
    bestmu_std = result.params['mu'].stderr
    bestA_std = result.params['A'].stderr
    bestdc_std = result.params['dc'].stderr    
    
    bestparams = np.array([bestsig,bestmu,bestA,bestdc])
    bestparams_std = np.array([bestsig_std,bestmu_std,bestA_std,bestdc_std])
    
    rchi = result.redchi    
    
    #return best values and std errors on the other parameters as well    
    
    return noiselessmodel, besttau, taustd, bestparams, bestparams_std, rchi



def tau_fitter_fix1D(data,nbins,fixval):
    #Fixed width fitter, need to pass it the fixed values
    #Also used to fix tau
    
    profile_peak = np.max(data)
    binpeak = np.argmax(data)  
    modelname = GxETrain1D
    model = Model(modelname)
            
    model.set_param_hint('nbins', value=nbins, vary=False)       
#    model.set_param_hint('sigma', value=fixval, vary=False)
    model.set_param_hint('sigma', value=15, vary=True, min =0, max = nbins)
    model.set_param_hint('mu', value=binpeak, vary=True, min=0, max = nbins)
    model.set_param_hint('A',value=profile_peak, vary=True)
#    model.set_param_hint('tau',value=200, vary=True, min=0)
    model.set_param_hint('tau1',value=fixval, vary=False)
#    model.set_param_hint('dc',value = 0, vary = False)
    model.set_param_hint('dc',value = 0, vary =True)
    pars = model.make_params()
    xax=np.linspace(1,nbins,nbins)
    #"""Fit data"""
    result = model.fit(data,pars,x=xax)
    print(result.fit_report(show_correl = False))
    
    noiselessmodel = result.best_fit
    besttau = result.best_values['tau1']
    taustd = result.params['tau1'].stderr  ##estimated 1 sigma error

    bestsig = result.best_values['sigma']
    bestmu = result.best_values['mu']
    bestA = result.best_values['A']
    bestdc = result.best_values['dc']
    
    bestsig_std = result.params['sigma'].stderr
    bestmu_std = result.params['mu'].stderr
    bestA_std = result.params['A'].stderr
    bestdc_std = result.params['dc'].stderr    
    
    bestparams = np.array([bestsig,bestmu,bestA,bestdc])
    bestparams_std = np.array([bestsig_std,bestmu_std,bestA_std,bestdc_std])
    
    rchi = result.redchi    
    
    #return best values and std errors on the other parameters as well    
    
    return noiselessmodel, besttau, taustd, bestparams, bestparams_std, rchi



    
def tau_fittermulti(data,nbins,guess_mean1, guess_mean2, guess_peak1, guess_peak2, maxmean):
#    binlen = len(data)
#    profile_peak1 = np.max(data[0:binlen/2])
#    profile_peak2 = np.max(data[binlen/2:binlen])
#    profile_peak1 = np.max(data[40:55])
#    profile_peak2 = np.max(data[55:60])
#    print profile_peak1
#    print profile_peak2
#    binpeak = np.argmax(data)
       
    modelname = GxETrainMulti
    model = Model(modelname)
            
    model.set_param_hint('nbins', value=nbins, vary=False)
    model.set_param_hint('ncomps', value=2, vary=False)       
    model.set_param_hint('sigmas1', value=15, vary=True, min =0, max = nbins/4)
    model.set_param_hint('sigmas2', value=15, vary=True, min =0, max = nbins/4)
    model.set_param_hint('mus1', value=guess_mean1, vary=True, min=660, max =720)
    model.set_param_hint('mus2', value=guess_mean2, vary=True, min=725, max = 750)
    model.set_param_hint('As1',value=guess_peak1, vary=True,min=0)
    model.set_param_hint('As2',value=guess_peak2, vary=True,min=0)
    model.set_param_hint('tau',value=10, vary=True, min=0)
    model.set_param_hint('dc',value = 0.1, vary = True)
#    val1, val2, val3, val4, val5, val6 = 20.25, 7.58, 689.5, 739.8, 1.05, 0.412  ##These are the mean values from the free fit
#    kt = 0.05
#    min1, min2, min3, min4, min5, min6 = (1-kt)*val1, (1-kt)*val2, (1-kt)*val3, (1-kt)*val4, (1-kt)*val5, (1-kt)*val6 
#    max1, max2, max3, max4, max5, max6 = (1+kt)*val1, (1+kt)*val2, (1+kt)*val3, (1+kt)*val4, (1+kt)*val5, (1+kt)*val6 
#    model.set_param_hint('nbins', value=nbins, vary=False)
#    model.set_param_hint('ncomps', value=2, vary=False)       
#    model.set_param_hint('sigmas1', value=val1*1.2, vary=True, min =min1, max = max1)
#    model.set_param_hint('sigmas2', value=val2*0.8, vary=True, min =min2, max = max2)
#    model.set_param_hint('mus1', value=val3*1.01, vary=True, min=min3, max =max3)
#    model.set_param_hint('mus2', value=val4*0.9, vary=True, min=min4, max = max4)
#    model.set_param_hint('As1',value=val5*0.6, vary=True,min=min5, max=max5)
#    model.set_param_hint('As2',value=val6*0.3, vary=True,min=min6, max=max6)
#    model.set_param_hint('tau',value=10, vary=True, min=0)
#    model.set_param_hint('dc',value = 0.1, vary = True)
    pars = model.make_params()
    xax=np.linspace(1,nbins,nbins)
    #"""Fit data"""
    result = model.fit(data,pars,x=xax)
    print(result.fit_report(show_correl = False))
    
    noiselessmodel = result.best_fit
    besttau = result.best_values['tau']
    taustd = result.params['tau'].stderr  ##estimated 1 sigma error

    bestsig1 = result.best_values['sigmas1']
    bestsig2 = result.best_values['sigmas2']
    bestmu1 = result.best_values['mus1']
    bestmu2 = result.best_values['mus2']
    bestA1 = result.best_values['As1']
    bestA2 = result.best_values['As2']
    
    bestdc = result.best_values['dc']    
    
    bestsig_std1 = result.params['sigmas1'].stderr
    bestsig_std2 = result.params['sigmas2'].stderr
    bestmu_std1 = result.params['mus1'].stderr
    bestmu_std2 = result.params['mus2'].stderr
    bestA_std1 = result.params['As1'].stderr
    bestA_std2 = result.params['As2'].stderr
    
    bestdc_std = result.params['dc'].stderr    
    
    bestparams = np.array([bestsig1,bestsig2,bestmu1,bestmu2, bestA1, bestA2,bestdc])
    bestparams_std = np.array([bestsig_std1,bestsig_std2, bestmu_std1,bestmu_std2, bestA_std1,bestA_std2,bestdc_std])
    
    rchi = result.redchi    
    
    return noiselessmodel, besttau, taustd, bestparams, bestparams_std, rchi


def tau_fittermultiloop(data,nbins,guess_mean1, guess_mean2):
    
    profile_peak1 = np.max(data[0:32])
    profile_peak2 = np.max(data[32:64])
    print profile_peak1
    print profile_peak2
#    binpeak = np.argmax(data)
       
    modelname = GxETrainMulti
    model = Model(modelname)
            
    model.set_param_hint('nbins', value=nbins, vary=False)
    model.set_param_hint('ncomps', value=2, vary=False)       
    model.set_param_hint('sigmas1', value=15, vary=True, min =0, max = nbins)
    model.set_param_hint('sigmas2', value=15, vary=True, min =0, max = nbins)
    model.set_param_hint('mus1', value=guess_mean1, vary=True, min=0, max = nbins/2)
    model.set_param_hint('mus2', value=guess_mean2, vary=True, min=0, max = nbins/2)
    model.set_param_hint('As1',value=profile_peak1, vary=True,min=0)
    model.set_param_hint('As2',value=profile_peak2, vary=True,min=0)
    model.set_param_hint('tau',value=200, vary=True, min=0)
    model.set_param_hint('dc',value = 0.1, vary = True)
    pars = model.make_params()
    xax=np.linspace(1,nbins,nbins)
    #"""Fit data"""
    result = model.fit(data,pars,x=xax)
    print(result.fit_report(show_correl = False))
    
    noiselessmodel = result.best_fit
    besttau = result.best_values['tau']
    taustd = result.params['tau'].stderr  ##estimated 1 sigma error

    return noiselessmodel, besttau, taustd

    
def tau_fittermulti1D(data,nbins,guess_mean1, guess_mean2):
    binlen = len(data)
    profile_peak1 = np.max(data[0:binlen/2])
    profile_peak2 = np.max(data[binlen/2:binlen])
    print profile_peak1
    print profile_peak2
#    binpeak = np.argmax(data)
       
    modelname = GxETrainMulti1D
    model = Model(modelname)
            
    model.set_param_hint('nbins', value=nbins, vary=False)
    model.set_param_hint('ncomps', value=2, vary=False)       
    model.set_param_hint('sigmas1', value=15, vary=True, min =0, max = nbins)
    model.set_param_hint('sigmas2', value=15, vary=True, min =0, max = nbins)
    model.set_param_hint('mus1', value=guess_mean1, vary=True, min=0, max = nbins)
    model.set_param_hint('mus2', value=guess_mean2, vary=True, min=0, max = nbins)
    model.set_param_hint('As1',value=profile_peak1, vary=True,min=0)
    model.set_param_hint('As2',value=profile_peak2, vary=True,min=0)
    model.set_param_hint('tau',value=200, vary=True, min=0)
    model.set_param_hint('dc',value = 0.1, vary = True)
    pars = model.make_params()
    xax=np.linspace(1,nbins,nbins)
    #"""Fit data"""
    result = model.fit(data,pars,x=xax)
    print(result.fit_report(show_correl = False))
    
    noiselessmodel = result.best_fit
    besttau = result.best_values['tau']
    taustd = result.params['tau'].stderr  ##estimated 1 sigma error

    return noiselessmodel, besttau, taustd



def tau_fitter_postfold(data,nbins,trainlength):
    
    profile_peak = np.max(data)
    binpeak = np.argmax(data)  
    modelname = GxESingleFoldDC
    model2 = Model(modelname)
            
    model2.set_param_hint('nbins', value=nbins, vary=False) 
    model2.set_param_hint('trainlength', value=trainlength, vary=False)       
    model2.set_param_hint('sigma', value=15, vary=True, min =0, max = nbins)
    model2.set_param_hint('mu', value=binpeak, vary=True, min=0, max = nbins)
    model2.set_param_hint('A',value=profile_peak, vary=True)
    model2.set_param_hint('tau',value=200, vary=True, min=0)
    model2.set_param_hint('dc',value = 0, vary = True)
    pars = model2.make_params()
    xax = np.arange(0,nbins,1)
    #"""Fit data"""
    result = model2.fit(data,pars,x=xax)
#    print(result.fit_report(show_correl = False))
    
    noiselessmodel = result.best_fit
    besttau = result.best_values['tau']
    taustd = result.params['tau'].stderr  ##estimated 1 sigma error

    return noiselessmodel, besttau, taustd


def tau_fitter_postfold_noDC(data,nbins,trainlength):
    
    profile_peak = np.max(data)
    binpeak = np.argmax(data)  
    modelname = GxESingleFoldnoDC
    model2 = Model(modelname)
            
    model2.set_param_hint('nbins', value=nbins, vary=False) 
    model2.set_param_hint('trainlength', value=trainlength, vary=False)       
    model2.set_param_hint('sigma', value=15, vary=True, min =0, max = nbins)
    model2.set_param_hint('mu', value=binpeak, vary=False, min=0, max = nbins)
    model2.set_param_hint('A',value=profile_peak, vary=True)
    model2.set_param_hint('tau',value=200, vary=True, min=0)
    pars = model2.make_params()
    xax = np.arange(0,nbins,1)
    #"""Fit data"""
    result = model2.fit(data,pars,x=xax)
#    print(result.fit_report(show_correl = False))
    
    noiselessmodel = result.best_fit
    besttau = result.best_values['tau']
    taustd = result.params['tau'].stderr  ##estimated 1 sigma error

    return noiselessmodel, besttau, taustd




    

def tau_1D_fitter(data,nbins):
    
    profile_peak = np.max(data)
    binpeak = np.argmax(data)  
    modelname = GxETrain1D
    model = Model(modelname)
            
    model.set_param_hint('nbins', value=nbins, vary=False)       
    model.set_param_hint('sigma', value=15, vary=True, min =0, max = nbins)
    model.set_param_hint('mu', value=binpeak, vary=True, min=0, max = nbins)
    model.set_param_hint('A',value=profile_peak, vary=True,min=0)
    model.set_param_hint('tau1',value=200, vary=True, min=0)
#    model.set_param_hint('tau1',value=166.792877, vary=False)
    model.set_param_hint('dc',value = 0, vary = True)
    pars = model.make_params()
    
    result = model.fit(data,pars,x=np.linspace(1,nbins,nbins))
#    print(result.fit_report(show_correl = False))
   
    noiselessmodel = result.best_fit
    besttau = result.best_values['tau1']
    taustd = result.params['tau1'].stderr  ##estimated 1 sigma error
    
    bestsig = result.best_values['sigma']
    bestmu = result.best_values['mu']
    bestA = result.best_values['A']
    bestdc = result.best_values['dc']
    
    bestsig_std = result.params['sigma'].stderr
    bestmu_std = result.params['mu'].stderr
    bestA_std = result.params['A'].stderr
    bestdc_std = result.params['dc'].stderr    
    
    bestparams = np.array([bestsig,bestmu,bestA,bestdc])
    bestparams_std = np.array([bestsig_std,bestmu_std,bestA_std,bestdc_std])
    
    """correlations with sigma"""    
    corsig = result.params['sigma'].correl
    #corA = result.params['A'].correl
    #corlist = [corsig,corA]
   
    
    rchi = result.redchi    
    
    #return best values and std errors on the other parameters as well    
    
    return result, noiselessmodel, besttau, taustd, bestparams, bestparams_std, rchi, corsig
#    return noiselessmodel, besttau, taustd



def tau_ani_fitter(data,nbins):
    
    profile_peak = np.max(data)
    print 'PP: %.3f' %profile_peak
    binpeak = np.argmax(data)
    print 'Bin P: %.3f' %binpeak    
    modelname = GxETrainAniDC
    modelAni = Model(modelname)    

    nr_cond = 6
    start1 = nbins/2.    
    incr = int(nbins/5.)
#    incr = 80.
   
    besttau1, tau1std, besttau2, tau2std = np.zeros(nr_cond),np.zeros(nr_cond),np.zeros(nr_cond),np.zeros(nr_cond)
    corsigs = []
    tau_percerrors, tau2_percerrors = np.zeros(nr_cond), np.zeros(nr_cond)
    bestsig, bestmu, bestA, bestdc, rchi = np.zeros(nr_cond),np.zeros(nr_cond),np.zeros(nr_cond),np.zeros(nr_cond),np.zeros(nr_cond)
    bestsig_std, bestmu_std, bestA_std, bestdc_std = np.zeros(nr_cond),np.zeros(nr_cond),np.zeros(nr_cond),np.zeros(nr_cond)

    resultAnis= []
    ##Run the fit over several initial conditions of tau values, and choose the one with the min. percentage error.
    ##This is done to avoid the instances in which some of the errors in the tau values blow up, based on the choice of initial condition
#    count = 0    
    for k in range(nr_cond):
            print "k: %d/%d" %(k,nr_cond)
#        for j in range(nr_cond):
#            starts1 = [nbins/8.,nbins/4.,nbins/2.,nbins/1.5] 
            modelAni.set_param_hint('nbins', value=nbins, vary=False)       
            modelAni.set_param_hint('sigma', value=15., vary=True, min =0, max = nbins)
            modelAni.set_param_hint('mu', value=binpeak, vary=True, min=0, max = nbins)
            modelAni.set_param_hint('A',value=profile_peak, vary=True)
            modelAni.set_param_hint('tau1',value=start1, vary=True, min=0, max = 1*nbins)
            modelAni.set_param_hint('tau2',value=(k+1)*incr, vary=True, min=0, max = 1*nbins)
            modelAni.set_param_hint('dc',value = 0.01, vary = True)
            parsAni = modelAni.make_params()
    #        print parsAni
    
            resultAni = modelAni.fit(data,parsAni,x=np.linspace(1,nbins,nbins))
            resultAnis.append(resultAni)
            print(resultAni.fit_report(show_correl = False))    
            noiselessmodel = resultAni.best_fit                
#            besttau1[count] = resultAni.best_values['tau1']
#            tau1std[count] = resultAni.params['tau1'].stderr  ##estimated 1 sigma error
#            besttau2[count] = resultAni.best_values['tau2']
#            tau2std[count] = resultAni.params['tau2'].stderr 
#            count += 1
            besttau1[k] = resultAni.best_values['tau1']
            tau1std[k] = resultAni.params['tau1'].stderr  ##estimated 1 sigma error
            besttau2[k] = resultAni.best_values['tau2']
            tau2std[k] = resultAni.params['tau2'].stderr 
 
            bestsig[k] = resultAni.best_values['sigma']
            bestmu[k] = resultAni.best_values['mu']
            bestA[k] = resultAni.best_values['A']
            bestdc[k] = resultAni.best_values['dc']
    
            bestsig_std[k] = resultAni.params['sigma'].stderr
            bestmu_std[k] = resultAni.params['mu'].stderr
            bestA_std[k] = resultAni.params['A'].stderr
            bestdc_std[k] = resultAni.params['dc'].stderr  
          
            rchi[k] = resultAni.redchi
            corsig = resultAni.params['sigma'].correl
            corsigs.append(corsig)

    tau_percerrors = tau1std/besttau1*100
    tau2_percerrors = tau2std/besttau2*100
    ave_err = (tau_percerrors+tau2_percerrors)/2.
    ave_err[ave_err == 0] = np.max(ave_err)
    ind_minerr = np.nanargmin(ave_err)
    print ave_err
    print ind_minerr    
    
    finaltau1,finaltau2 = besttau1[ind_minerr],besttau2[ind_minerr]
    finaltau1std, finaltau2std = tau1std[ind_minerr],tau2std[ind_minerr]
    noiselessmodel = resultAnis[ind_minerr].best_fit
    
    finalsig, finalmu, finalA, finaldc = bestsig[ind_minerr],bestmu[ind_minerr],bestA[ind_minerr],bestdc[ind_minerr]
    finalsig_std, finalmu_std, finalA_std, finaldc_std = bestsig_std[ind_minerr],bestmu_std[ind_minerr],bestA_std[ind_minerr],bestdc_std[ind_minerr]
    
    finalrchi = rchi[ind_minerr]    
    
    finalparams = np.array([finalsig,finalmu,finalA,finaldc])
    finalparams_std = np.array([finalsig_std,finalmu_std,finalA_std,finaldc_std])
           
    finalcorsig = corsigs[ind_minerr]
    
    return resultAnis[ind_minerr], noiselessmodel, finaltau1, finaltau1std, finaltau2, finaltau2std, finalparams, finalparams_std, finalrchi, finalcorsig
    

def make_mono(observingfreq,bandwidth):
    mono_freq = 0.5*np.sqrt(bandwidth**2 + 4*observingfreq**2)
    return mono_freq



def simulate(pulseperiod, tausec,dutycycle,specindex,freq,freqlow,nbins,snr,simu):
        # tausec will either be a single value (iso, 1D) or a 2 column array (aniso)  
        # assume tausec in seconds
        # make sure freqs are entered in GHz
        m = float(nbins/3)                      # Let the mean of the pulse be at a 1/4 overall bins
        w50 = float((dutycycle/100)*nbins)      # FWHM
        s = w50/(2*np.sqrt(2*np.log(2)))    # sigma calculated through the dutycycle
        a = 1                               # amplitude. at some point this will be replaced by the intrinsic spectrum. 
        bs = pulseperiod/nbins
        sb = 1./bs
    
        bins, profile = makeprofile(nbins = nbins, ncomps = 1, amps = a, means = m, sigmas = s)
        xax = np.linspace(1,nbins,nbins)
        profile_intr = pow(freqlow,specindex)*pow(freq,-specindex)*profile
        profile_intr_norm = profile_intr/np.sum(profile)     
        if simu in  ('iso','ISO', 'Iso'):
            scat = psrscatter(broadfunc(xax,sb*tausec),pulsetrain(3, bins, profile_intr_norm))       
            climb, observed_nonoise, rec, flux = extractpulse(scat, 2, nbins)
            peak = np.max(observed_nonoise)
            noise = np.random.normal(0,peak/snr,nbins)
            observedadd = climb + noise
        elif simu in ('onedim','1D','Onedim'):
            scat = psrscatter(broadfunc1D(xax,sb*tausec),pulsetrain(3, bins, profile_intr_norm))       
            climb, observed_nonoise, rec, flux = extractpulse(scat, 2, nbins)
            peak = np.max(observed_nonoise)
            noise = np.random.normal(0,peak/snr,nbins)
            observedadd = climb + noise
        elif simu in ('aniso','ANISO','Aniso'):
            scat = psrscatter(broadfunc2part(xax,sb*tausec[0],sb*tausec[1]),pulsetrain(3, bins, profile_intr_norm))       
            climb, observed_nonoise, rec, flux = extractpulse(scat, 2, nbins)
            peak = np.max(observed_nonoise)
            noise = np.random.normal(0,peak/snr,nbins)
            observedadd = climb + noise
        else:
            print "Invalid simu string. Choose from: iso, onedim, aniso"
        return profile_intr_norm, observedadd
           


def simulatetrain(pulseperiod, tausec,dutycycle,nbins,snr,simu,tl):
        # tausec will either be a single value (iso, 1D) or a 2 column array (aniso)  
        # assume tausec in seconds
        # make sure freqs are entered in GHz
        m = float(nbins/3)                      # Let the mean of the pulse be at a 1/4 overall bins
        w50 = float((dutycycle/100)*nbins)      # FWHM
        s = w50/(2*np.sqrt(2*np.log(2)))    # sigma calculated through the dutycycle
        a = 1                               # amplitude. at some point this will be replaced by the intrinsic spectrum. 
        bs = pulseperiod/nbins
        sb = 1./bs    
        bins, profile = makeprofile(nbins = nbins, ncomps = 1, amps = a, means = m, sigmas = s)
#        peak = np.max(profile)
#        noise = np.random.normal(0,peak/snr,nbins)        
#        profile_intr = pow(freqlow,specindex)*pow(freq,-specindex)*profile
#        profile_intr_norm = profile_intr/np.sum(profile)
#        profile_intr_normN = profile_intr_norm + noise
#        profileN = profile + noise
        xax = np.linspace(1,nbins,nbins)
        
        if tausec == 0:
            scat = pulsetrain(tl, bins, profile)  
#            scatN = pulsetrain(tl, bins, profileN)
            climb, observed_nonoise, rec, flux = extractpulse(scat, 2, nbins)
            peak = np.max(observed_nonoise)
            noise = np.random.normal(0,peak/snr,nbins)
            noisetr = np.random.normal(0,peak/snr,tl*nbins)
            climbN = climb + noise
            scatN = scat + noisetr
#            climbN, observed_nonoise, rec, flux = extractpulse(scatN, 2, nbins)
            #peak = np.max(observed_nonoise)
            #noise = np.random.normal(0,peak/snr,nbins)
            #noisetr = np.random.normal(0,peak/snr,tl*nbins)
        else:   
            if simu in  ('iso','ISO', 'Iso'):
                scat = psrscatter(broadfunc(xax,sb*tausec),pulsetrain(tl, bins, profile))
                climb, observed_nonoise, rec, flux = extractpulse(scat, 2, nbins)           
                peak = np.max(observed_nonoise)
                noise = np.random.normal(0,peak/snr,nbins)
                noisetr = np.random.normal(0,peak/snr,tl*nbins)
                climbN = climb + noise
                scatN = scat + noisetr
#                #observedadd = climb + noise
#            elif simu in ('onedim','1D','Onedim'):
#                scat = psrscatter(broadfunc1D(xax,sb*tausec),pulsetrain(tl, bins, profile_intr_norm))       
#                climb, observed_nonoise, rec, flux = extractpulse(scat, 2, nbins)
#                peak = np.max(observed_nonoise)
#                noise = np.random.normal(0,peak/snr,nbins)
##                noisetr = np.random.normal(0,peak/snr,tl*nbins)
#                observedadd = climb + noise
#            elif simu in ('aniso','ANISO','Aniso'):
#                scat = psrscatter(broadfunc2part(xax,sb*tausec[0],sb*tausec[1]),pulsetrain(tl, bins, profile_intr_norm))       
#                climb, observed_nonoise, rec, flux = extractpulse(scat, 2, nbins)
#                peak = np.max(observed_nonoise)
#                noise = np.random.normal(0,peak/snr,nbins)
##                noisetr = np.random.normal(0,peak/snr,tl*nbins)
#                observedadd = climb + noise
            else:
                print "Invalid simu string. Choose from: iso, onedim, aniso"
        return profile, scat, scatN, climb, climbN



#def tau_ani_fitter_2part(data,nbins):
#    
#    profile_peak = np.max(data)
#    print 'PP: %.3f' %profile_peak
#    binpeak = np.argmax(data)
#    print 'Bin P: %.3f' %binpeak    
#    modelname = GxETrainAniDC2part
#    modelAni = Model(modelname)
#            
#    modelAni.set_param_hint('nbins', value=nbins, vary=False)       
#    modelAni.set_param_hint('sigma', value=15., vary=True, min =0, max = nbins)
#    modelAni.set_param_hint('mu', value=binpeak, vary=True, min=0, max = nbins)
#    modelAni.set_param_hint('A',value=profile_peak, vary=True)
#    modelAni.set_param_hint('tau1',value=5.0, vary=True, min=0, max = 1*nbins)
#    modelAni.set_param_hint('tau2',value=100.0, vary=True, min=0, max = 1*nbins)
#    modelAni.set_param_hint('dc',value = 0, vary = True)
#    parsAni = modelAni.make_params()
#    
#    #"""Fit data"""
#    resultAni = modelAni.fit(data,parsAni,x=np.linspace(1,nbins,nbins))
#    print(resultAni.fit_report(show_correl = False))
#    
#    noiselessmodel = resultAni.best_fit
#    besttau1 = resultAni.best_values['tau1']
#    tau1std = resultAni.params['tau1'].stderr  ##estimated 1 sigma error
#    besttau2 = resultAni.best_values['tau2']
#    tau2std = resultAni.params['tau2'].stderr  
#    
#    return noiselessmodel, besttau1, tau1std, besttau2, tau2std
#
#
#def tau_ani_fitter_ive(data,nbins):
#    
#    profile_peak = np.max(data)
#    print 'PP: %.3f' %profile_peak
#    binpeak = np.argmax(data)
#    print 'Bin P: %.3f' %binpeak    
#    modelname = GxETrainAniDCive
#    modelAni = Model(modelname)
#            
#    modelAni.set_param_hint('nbins', value=nbins, vary=False)       
#    modelAni.set_param_hint('sigma', value=15., vary=True, min =0, max = nbins)
#    modelAni.set_param_hint('mu', value=binpeak, vary=True, min=0, max = nbins)
#    modelAni.set_param_hint('A',value=profile_peak, vary=True)
#    modelAni.set_param_hint('tau1',value=5.0, vary=True, min=0, max = 1*nbins)
#    modelAni.set_param_hint('tau2',value=100.0, vary=True, min=0, max = 1*nbins)
#    modelAni.set_param_hint('dc',value = 0, vary = True)
#    parsAni = modelAni.make_params()
#    
#    #"""Fit data"""
#    resultAni = modelAni.fit(data,parsAni,x=np.linspace(1,nbins,nbins))
#    print(resultAni.fit_report(show_correl = False))
#    
#    noiselessmodel = resultAni.best_fit
#    besttau1 = resultAni.best_values['tau1']
#    tau1std = resultAni.params['tau1'].stderr  ##estimated 1 sigma error
#    besttau2 = resultAni.best_values['tau2']
#    tau2std = resultAni.params['tau2'].stderr  
#    
#    return noiselessmodel, besttau1, tau1std, besttau2, tau2std
#


    
def find_rms(data,nbins):
    windowsize = 32
    windows = int(nbins/windowsize)
    rms_loc = np.zeros(windows)
    for i in range(windows):
        start = i*windowsize
        end = start + windowsize
        rms_loc[i] = np.std(data[start:end])
    return np.min(rms_loc)
    
    
def find_modelflux(model,nbins):
#note: this function is set up to find flux from model and not from noisy data.
        start = 0
        end = nbins
        zerobpulse = model[start:end]-np.min(model[start:end])
#        rectangle = np.min(model[start:end])*nbins
        flux = np.sum(zerobpulse[start:end])
        return flux/nbins

def find_peaksnr_smooth(data,rms):
    boxsize = int(0.05*len(data))
    smootheddata = smooth(data,boxsize)
    peak = np.max(smootheddata)
    snr = peak/rms   
    return snr

def find_peaksnr(data,rms):
    peak = np.max(data)
    snr = peak/rms   
    return snr
        

def generate_scattered(nbins,period,dutycycle,tau,freq,snr):   
#   period is pulseperiod in seconds
#   dutycycle as a % of the overall pulseperiod
#   tau is characteristic scattering time in seconds
#   freq in GHz
    m = float(nbins/4)               # Let the mean of the pulse be at a 1/4 overall bins
    w50 = float((dutycycle/100)*nbins)   # FWHM
    s = w50/(2*np.sqrt(2*np.log(2))) # sigma calculated through the dutycycle
    a = 1                            # amplitude. at some point this will be replaced by the intrinsic spectrum. 
    taubins = tau*nbins/period
    xax = np.arange(0,nbins,1)
    
    bins, profile = makeprofile(nbins = nbins, ncomps = 1, amps = a, means = m, sigmas = s)
    scat = psrscatter(broadfunc(xax,taubins),pulsetrain(3, bins, profile))      
    climb, observed_nonoise, rec, flux = extractpulse(scat, 2, nbins)
    peak = np.max(observed_nonoise)
    noise = np.random.normal(0,peak/snr,nbins)
    noisyscatteredprofile = climb + noise
    return noisyscatteredprofile
    
    
def autocorr(timeseries):
    timeseries -= np.mean(timeseries)
    autocorr_f = np.correlate(timeseries, timeseries, mode='full')
    temp = autocorr_f[autocorr_f.size/2:]/autocorr_f[autocorr_f.size/2]
    return temp


def second_smallest(numbers):
    m1, m2 = float('inf'), float('inf')
    for x in numbers:
        if x <= m1:
            m1, m2 = x, m1
        elif x < m2:
            m2 = x
    return m2

def DM_Model(x,C,k):
    highest_fr = x[-1]
    deltime = -C*((1/highest_fr**2) + (1/np.power(x,2))) + k
    return deltime
  

def DM_checker(freqMHz, time): ##(freq in MHz, time in sec)
    delta_t = (time[-1] - time[0])
    delta_nu = -(1/(freqMHz[-1]**2) - 1/(freqMHz[0]**2))
    DMconst = 4148.808
    DM = float(delta_t)/float(delta_nu)/DMconst
    return DM
#    
    
def tauatfreq(oldfreq,oldtau,newfreq,specindex):
    newtau = oldtau*(newfreq)**(-specindex)/(oldfreq**(-specindex))
    return newtau


def fluxatfreq(oldfreq,oldflux,newfreq,specindex):
    newflux = oldflux*(newfreq)**(-specindex)/(oldfreq**(-specindex))
    return newflux    


def linking_specind(tau1, freq1, tau2, freq2):
    specind = np.log(tau1/tau2)/np.log(freq2/freq1)
    return specind
#def timeshift_fromDM(DM,freqshifts):
#    delta_nu = -(1/(freqMHz[-1]**2) - 1/(freqMHz[0]**2))
#    delta_t = (time[-1] - time[0])
#    DMconst = 4148.808
#    DM = float(delta_t)/float(delta_nu)/DMconst
#    return DM



def Ramanch_model(x,c1,c2,gamma):
    tauarray_ms = c1*np.power(x,2.2) + c2*np.power(x,gamma)
    return tauarray_ms


def Ramanch_fit(x,data):

    model = Model(Ramanch_model)
            
#    model.set_param_hint('DMarray', value=DMarray, vary=False)       
    model.set_param_hint('c1', value=2.26e-7, vary=True,min=1e-8,max=1e-6)
    model.set_param_hint('c2',value=2.26e-7*0.00205, vary=True)
    model.set_param_hint('gamma',value = 2.2+1.74, vary = True,min=3.0, max=5.5)
#    model.set_param_hint('k',value = 0.001, vary = True)
#    model.set_param_hint('gamma2',value = 2.2, vary = True)
    pars = model.make_params()
    

    result = model.fit(data,pars,x=x,weights=1.0)
    print(result.fit_report(show_correl = False))
    
    Rmodel = result.best_fit
    bestc1 = result.best_values['c1']
    bestc2 = result.best_values['c2']
    bestgam = result.best_values['gamma']
#    bestgam2 = result.best_values['k']

    bestparams = np.array([bestgam,bestc1,bestc2])
    #bestparams_std = np.array([bestsig_std,bestmu_std,bestA_std,bestdc_std])
    
    rchi = result.redchi    
    
    #return best values and std errors on the other parameters as well    
    
    return Rmodel, bestparams, rchi


def Bhat_model(x,k1,k2,k3):
    #x here is the log(DM)
    logtau_ms = k1 + k2*x + k3*np.power(x,2)
    return logtau_ms


def Bhat_fit(xdat,data):

    model = Model(Bhat_model)
    model.set_param_hint('k1', value=-7.0, vary=True,min=-10.0, max=-5.0)
    model.set_param_hint('k2',value=2.0, vary=True,min=0.0, max=4.0)
    model.set_param_hint('k3',value=0.8, vary=True,min=0.0,max=10.)

    pars = model.make_params()  

    result = model.fit(data,pars,x=xdat)
    print(result.fit_report(show_correl = False))
    
    Bmodel = result.best_fit
    bestk1 = result.best_values['k1']
    bestk2 = result.best_values['k2']
    bestk3 = result.best_values['k3']

    bestparams = np.array([bestk1,bestk2,bestk3])
    #bestparams_std = np.array([bestsig_std,bestmu_std,bestA_std,bestdc_std])
    
    rchi = result.redchi    
    
    #return best values and std errors on the other parameters as well    
    
    return Bmodel, bestparams, rchi




#def autocorr(x):
#    result = np.correlate(x, x, mode = 'full')
#    maxcorr = np.argmax(result)
#    result = result / result[maxcorr]     # <=== normalization
#    return result[result.size/2:]    
