#!/Users/mgeyer/anaconda2/bin/python
"""
Utility functions for scatter fitting routines
"""
import os
import argparse
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from lmfit import Model, conf_interval, printfuncs
from lmfit import minimize, Parameter, Parameters, fit_report
from lmfit.models import LinearModel, PowerLawModel, ExponentialModel, QuadraticModel, GaussianModel 


def writePDVtimeSeries(dspec, freqs, tspan, nch=1, src='FRB', ofn='timeseries.ascii'): 
    """A function to mimic a simple time series output of running pdV -kAT,
    assumes nsub=1, npol=1, nbin based on input dynamic spectrum
    This is used to input into Marisa's code
    dspec: 2-d numpy array, dynamic spectrum (time, freq)
    freqs: 1-d numpy array, frequencies in MHz
    tspan: float, obs time in seconds
    nch: int, number of channels to output
    src: string, source name
    ofn: string, output file name
    """
    fDecFactor = int(freqs.shape[0] / nch) # frequency decimation factor
    decdspec = np.zeros((dspec.shape[0], nch))
    if dspec.shape[1] % fDecFactor > 0:
        print 'WARNING: the number of frequency channels %i is not an integer multiple of the number of subband channels %i, %i frequency channels will be dropped'%(dspec.shape[1], nch, dspec.shape[1] % fDecFactor)
        dspecClip = dspec[:,:dspec.shape[1] - (dspec.shape[1] % fDecFactor)]
        freqsClip = freqs[:dspec.shape[1] - (dspec.shape[1] % fDecFactor)]
    else:
        dspecClip = dspec.copy()
        freqsClip = freqs.copy()
    dspecClip -= dspecClip.mean()

    # average down in frequency
    decdspec = np.mean(dspecClip.reshape(dspecClip.shape[0], nch, fDecFactor), axis=2)
    decfreqs = np.mean(freqsClip.reshape(nch, fDecFactor), axis=1)
    bw = freqs[fDecFactor] - freqs[0]

    # overall header
    # # File: B0611+22_L116889 Src: J0614+2229 Nsub: 1 Nch: 8 Npol: 1 Nbin: 1024 RMS: 93.8508
    hdrStr = '# File: %s Src: %s Nsub: 1 Nch: %i Npol: 1 Nbin: %i RMS: %f'%(src, src, nch, decdspec.shape[0], np.std(decdspec))

    for chIdx in np.arange(nch):
        # write channel header
        # # MJD(mid): 56384.680540767529 Tsub: 599.953 Freq: 114.752 BW: 9.76562
        # NOTE: MJD is fixed, and not the real MJD
        hdrStr += '\n# MJD(mid): 56384.680540767529 Tsub: %f Freq: %f BW: %f'%(tspan, decfreqs[chIdx], bw)
        
        #write ascii profile
        # nsub nch nbin stokesI
        # 0 1 436 62.855
        for tIdx in np.arange(decdspec.shape[0]):
            hdrStr += '\n0 %i %i %f'%(chIdx, tIdx, decdspec[tIdx, chIdx])

    ofh = open(ofn, 'w')
    ofh.write(hdrStr)
    ofh.close()
    

    
### Read ascii files

def read_headerfull(filepath):
    """Read asciis as produced by pdv -KAt"""
    f = open(filepath)
    lines = f.readlines()
    header0 = lines[0]
    header1 = lines[1]
    h0_lines = header0.split()
    if h0_lines[0] == '#':
        h0_lines = h0_lines[1:len(h0_lines)]
    else:
        h0_lines = h0_lines    
    file_name = h0_lines[1]
    pulsar_name = h0_lines[3]
    nsub = int(h0_lines[5])
    nch = int(h0_lines[7])
    npol = int(h0_lines[9])
    nbins = int(h0_lines[11])
    rms = float(h0_lines[13])
    h1_lines = header1.split()
    tsub = float(h1_lines[4])  
#    return file_name, pulsar_name, nsub, nch, npol, nbins, rms
    return pulsar_name, nch, nbins, nsub, rms, tsub



def read_data(filepath, profilenumber, nbins):
    """Read data as produced by pdv -KAt"""
    d = open(filepath)
    lines = d.readlines()
    
    profile_start = 2+profilenumber*(nbins+1)
    profile_end = profile_start + nbins
    
    lines_block = lines[profile_start:profile_end]

    if lines[profile_start-1].split()[0] == '#':
        freqc = float(lines[profile_start-1].split()[6])
        bw = float(lines[profile_start-1].split()[8])
        freqm = 10**((np.log10(freqc+ bw/2.)+ np.log10(freqc - bw/2.))/2)
    else:
        freqc = float(lines[profile_start-1].split()[5])
        bw = float(lines[profile_start-1].split()[7])
        freqm = 10**((np.log10(freqc+ bw/2.)+ np.log10(freqc - bw/2.))/2)
    datalist = []
    for i in range(nbins):
        data= float(lines_block[i].split()[3])
        datalist.append(data)

    return np.array(datalist), freqc, freqm

    
    
    

def find_rms(data,nbins):
    windowsize = 32 
    windows = int(nbins/windowsize)
    rms_loc = np.zeros(windows)
    for i in range(windows):
        start = i*windowsize
        end = start + windowsize
        rms_loc[i] = np.std(data[start:end])
    return np.min(rms_loc)

def smooth(y, box_pts):
    gauss = np.ones(box_pts)
#    box = np.ones(box_pts)/box_pts
    sigma = (1./6.)*box_pts
    mean = box_pts/2.
    for i in range(box_pts):
        gauss[i] = (1/(sigma*np.sqrt(2*np.pi)))*np.exp(-(i-mean)**2/(2*sigma**2))
    y_smooth = np.convolve(y, gauss, mode='same')
    return y_smooth


def find_peaksnr_smooth(data,rms):
    if rms == 0:
        snr = float('nan')
    else:
        boxsize = int(0.05*len(data))
        smootheddata = smooth(data,boxsize)
        peak = np.max(smootheddata)
        snr = peak/rms
    return snr

def find_peaksnr(data,rms):
    if rms == 0:
        snr = float('nan')
    else:
        peak = np.max(data)
        snr = peak/rms
    return snr 

def tauatfreq(oldfreq,oldtau,newfreq,specindex):
    newtau = oldtau*(newfreq)**(-specindex)/(oldfreq**(-specindex))
    return newtau

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
    bins = profile.shape[0]
    out = scattered[0:bins]
    return out
    
def step(x):
    return 1 * (x >= 0)

### Broadening functions
# 1. Isotropic scattering

def broadfunc(x,tau):
    tau = float(tau)
    broadfunc = (1/tau)*np.exp(-x/tau)*step(x)
    return broadfunc   

# 2. Extremely anisotropic scattering

def broadfunc1D(x,tau):
    broadfunc1 = (1.0/np.sqrt(x*tau*np.pi))*np.exp(-x/tau)
    return broadfunc1


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


def returnclimb(x,mu,sigma,A,tau,dc,nbins):
    mu, sigma, A, tau = float(mu),float(sigma), float(A), float(tau)
    bins, profile = makeprofile(nbins = nbins, ncomps = 1, amps = A, means = mu, sigmas = sigma)
    binstau = np.linspace(1,nbins,nbins)
    scat = psrscatter_noconserve(broadfunc(binstau,tau),pulsetrain_bins(3, nbins, profile))
    climb, observed_nonoise, rec, flux = extractpulse(scat, 2, nbins)
    return np.min(climb)


def find_modelflux(model,nbins):
    #note: this function is set up to find flux from model and not from noisy data.
    start = 0
    end = nbins
    zerobpulse = model[start:end]-np.min(model[start:end])
    flux = np.sum(zerobpulse[start:end])
    return flux/nbins



def GxETrain(x,mu,sigma, A, tau, dc, nbins):
#This model convolves a pulsetrain with a broadening function
#It extracts one of the last convolved profiles, subtracts the climbed baseline and then adds noise to it
    mu, sigma, A, tau = float(mu),float(sigma), float(A), float(tau)
    bins, profile = makeprofile(nbins = nbins, ncomps = 1, amps = A, means = mu, sigmas = sigma)
    binstau = np.linspace(1,nbins,nbins)
    scat = psrscatter_noconserve(broadfunc(binstau,tau),pulsetrain_bins(3, nbins, profile))
#    scat = psrscatter(broadfunc(binstau,tau),pulsetrain_bins(3, nbins, profile))  
    climb, observed_nonoise, rec, flux = extractpulse(scat, 2, nbins)
    return observed_nonoise + dc
    
def GxETrain1D(x,mu, sigma, A, tau1, dc, nbins):
    mu, sigma, A, tau1 = float(mu),float(sigma), float(A), float(tau1)
    bins, profile = makeprofile(nbins = nbins, ncomps = 1, amps = A, means = mu, sigmas = sigma)
    binstau = np.linspace(1,nbins,nbins)
    scat = psrscatter(broadfunc1D(binstau,tau1),pulsetrain(3, bins, profile))
    climb, observed_nonoise, rec,flux = extractpulse(scat, 2, nbins)
    return observed_nonoise + dc
    
    
    
def tau_fitter(data,nbins, verbose=True):
    profile_peak = np.max(data)
    binpeak = np.argmax(data)  
    modelname = GxETrain
    model = Model(modelname)
                 
    model.set_param_hint('nbins', value=nbins, vary=False)            
    model.set_param_hint('sigma', value=15, vary=True, min =0, max = nbins)
    model.set_param_hint('mu', value=binpeak, vary=True, min=0, max = nbins)
    model.set_param_hint('A',value=profile_peak, vary=True, min=0)
    model.set_param_hint('tau',value=200, vary=True, min=0)
    model.set_param_hint('dc',value = 0, vary = True)
    pars = model.make_params()
    xax=np.linspace(1,nbins,nbins)

    #"""Fit data"""
    result = model.fit(data,pars,x=xax)
    if verbose == True:
        print(result.fit_report(show_correl = True))
    else:
        print "To see fit report, use verbose=True"
    
    noiselessmodel = result.best_fit
    besttau = result.best_values['tau']
    taustd = result.params['tau'].stderr  ##estimated 1 sigma error
    if taustd == None:
       taustd = 0

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
    if taustd == None:
       taustd = 0
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

    return result, noiselessmodel, besttau, taustd, bestparams, bestparams_std, rchi, corsig
       
    
    
    
def produce_taufits(filepath,meth='iso',pulseperiod=None,snr_cut=None,
        verbose=True, plotparams=False, plotflux=False, savefigure=False, noshow=False):
        pulsar, nch, nbins,nsub, lm_rms, tsub = read_headerfull(filepath)

        if verbose == True:
            verboseTag = True
        else:
            verboseTag = False
        
        print0 = "Pulsar name: %s" %pulsar
        print1 = "Number of freq. channels: %d \nFreq channels will be labeled 0 - %d" %(nch, nch-1)
        print2 = "Number of bins: %d" %nbins
        print3 = "RMS: %.2f" %lm_rms
        print4 = "Tsub: %.2f sec" %tsub 
        for k in range(5):
              print eval('print{0}'.format(k))
        print"--------------------------------------------------------"
        
        if pulseperiod==None:
            ## Define time axis, and time/bins conversions
            print ("Using Tsub in header to convert bins to time. Assumption is that tsub corresponds to 1.0 phase, corresponding to nbins.  This should be adapted for search data.")
            pulseperiod = tsub
        else:
            pulseperiod = pulseperiod #set to provided pulseperiod in seconds
        
        profilexaxis = np.linspace(0,pulseperiod,nbins)
        pbs = pulseperiod/nbins
        tbs = tsub/nbins
        """Initialise vector outputs"""
        obtainedtaus, lmfittausstds = [], []
        """freqmsMHz will correctly associate scattering time values 
        (tau) with frequency, taking into account frequency integration 
        across a sub-band. Whereas freqcsMHz is the centre freq. to the subband"""
        
        freqmsMHz, freqcsMHz = [], []    
        noiselessmodels =[]
        results, datas, comp_SNRs, comp_rmss = [], [], [], []
        redchis, paramset, paramset_std, correls = [], [], [], []
        
        halfway = nbins/2.

        for i in range(nch):
            print"--------------------------------------------------------"
            print "Channel %d" %i
            """Read in (pdv) data""" 
            data, freqc, freqm = read_data(filepath,i,nbins)
            freqmsMHz.append(freqm)
            freqcsMHz.append(freqc)
            # roll the data of lowest freq channel to middle of bins 
            if i ==0:
                peakbin = np.argmax(data)
                shift = int(halfway -int(peakbin))
                if verboseTag:
                    print 'peak bin at lowest freq channel:%d' %peakbin
            else:
                peakbin = peakbin
                shift = int(halfway - int(peakbin))
            data = np.roll(data,shift)
            if verboseTag:
                print "Rolling data by -%d bins" %shift
            comp_rms = find_rms(data,nbins)

            if meth is None:
                        print "No fitting method was chosen. Will default to an isotropic fitting model. \n Use option -m with 'onedim' to change."
                        result, noiselessmodel, besttau, taustd, bestparams, bestparams_std, redchi, corsig = tau_fitter(data,nbins,verbose=verboseTag)

            elif meth == 'iso':
                        result, noiselessmodel, besttau, taustd, bestparams, bestparams_std, redchi, corsig = tau_fitter(data,nbins,verbose=verboseTag)

            elif meth == 'onedim':
                        result, noiselessmodel, besttau, taustd, bestparams, bestparams_std, redchi, corsig = tau_1D_fitter(data,nbins)         

            comp_SNR_model = find_peaksnr(noiselessmodel,comp_rms)

            if verboseTag:
                print 'Estimated SNR (from model peak and data rms): %.2f' % comp_SNR_model
            comp_SNR =  find_peaksnr_smooth(data,comp_rms)
            print 'Estimated SNR (from data peak and rms): %.2f' % comp_SNR
            print 'Channel Tau (ms): %.2f \pm %.2f ms' %(besttau,taustd)
            
           
            obtainedtaus.append(besttau)
            lmfittausstds.append(taustd)
            noiselessmodels.append(noiselessmodel)
            results.append(result)
            datas.append(data)
            comp_SNRs.append(comp_SNR)
            #new:
            comp_rmss.append(comp_rms)
            redchis.append(redchi)
            paramset.append(bestparams)
            paramset_std.append(bestparams_std)
        #    if plotflux == True:
        #        correls.append(corsig)
        
        
        #if plotflux == True:
        #    cor_sigA = np.zeros(len(correls))
        #    for i in range(len(correls)):
        #        cor_sigA[i] = correls[i]['A']

   
        paramset = np.transpose(paramset)
        paramset_std = np.transpose(paramset_std)
         
        """Next check if any of the subbands contain only zeros. This happens with high RFI excision in LOFAR bands"""
        zero_ch = []
        for i in range(nch):
            all_zeros = not np.any(datas[i])
            if all_zeros:
                zero_ch.append(i)
        
        print"--------------------------------------------------------"

        if zero_ch:
            print "\n"
            print "%d channels have all zeroes (channels(s):" %len(zero_ch), zero_ch,  ") and will be removed."
            if verboseTag:
                print "All zero channels are assigned SNR of 0"

          
        if snr_cut: 
            print "Using SNR cutoff of %.2f" %snr_cut
            comp_SNRs = np.nan_to_num(comp_SNRs)
            (ind_lowSNR,) = np.where(np.array(comp_SNRs) < snr_cut)
            print "SNR cutoff will exclude %d channels (channel(s): %s)" %(len(ind_lowSNR), ind_lowSNR)

            
            data_highsnr = np.delete(np.array(datas),ind_lowSNR,0)
            model_highsnr = np.delete(np.array(noiselessmodels),ind_lowSNR,0)
            taus_highsnr = np.delete(np.array(obtainedtaus),ind_lowSNR)
            lmfitstds_highsnr = np.delete(np.array(lmfittausstds),ind_lowSNR)
            freqMHz_highsnr = np.delete(np.array(freqmsMHz),ind_lowSNR)
            #New:
            comp_rmss_highsnr = np.delete(np.array(comp_rmss),ind_lowSNR)
            redchis_highsnr = np.delete(np.array(redchis),ind_lowSNR)
            #corsigA_highsnr = np.delete(np.array(cor_sigA),ind_lowSNR)
            
            paramset_highsnr = np.zeros([len(paramset),len(data_highsnr)])
            paramsetstd_highsnr = np.zeros([len(paramset),len(data_highsnr)])                
            for i in range(len(paramset)):
                paramset_highsnr[i]= np.delete(paramset[i],ind_lowSNR)
                paramsetstd_highsnr[i]= np.delete(paramset_std[i],ind_lowSNR)
                
            
        elif (snr_cut == None) and (zero_ch != []):       
            print "Used no SNR cutoff"          
            """Rename array to be same as when cut-off is used"""
            """If no SNR cutoff is used, remove channels with all zeroes 
            -- these will automatically be removed by any snr_cut > 0"""
            data_highsnr = np.delete(np.array(datas),zero_ch,0)
            model_highsnr = np.delete(np.array(noiselessmodels),zero_ch,0)
            taus_highsnr = np.delete(np.array(obtainedtaus),zero_ch)
            lmfitstds_highsnr = np.delete(np.array(lmfittausstds),zero_ch)
            freqMHz_highsnr = np.delete(np.array(freqmsMHz),zero_ch)
            # New:
            comp_rmss_highsnr = np.delete(np.array(comp_rmss),zero_ch)
            redchis_highsnr = np.delete(np.array(redchis),zero_ch)
            #corsigA_highsnr = np.delete(np.array(cor_sigA),zero_ch)
              
            paramset_highsnr = np.zeros([len(paramset),len(data_highsnr)])
            paramsetstd_highsnr = np.zeros([len(paramset),len(data_highsnr)])     
            for i in range(len(paramset)):
                paramset_highsnr[i]= np.delete(paramset[i],zero_ch)
                paramsetstd_highsnr[i]= np.delete(paramset_std[i],zero_ch)
                
                
        else:
            print "Used no SNR cutoff and there are no empty channels"          
            data_highsnr = np.array(datas)
            model_highsnr = np.array(noiselessmodels)
            taus_highsnr = np.array(obtainedtaus)
            lmfitstds_highsnr = np.array(lmfittausstds)
            freqMHz_highsnr = np.array(freqmsMHz)
            # New:
            comp_rmss_highsnr = np.array(comp_rmss)
            redchis_highsnr = np.array(redchis)
            paramset_highsnr = np.array(paramset)
            paramsetstd_highsnr = np.array(paramset_std)
            
            
            
        
        taussec_highsnr = taus_highsnr*pbs
        lmfitstdssec_highsnr = lmfitstds_highsnr*pbs
        number_of_plotted_channels = len(data_highsnr)
        npch = number_of_plotted_channels
        print "Will plot remaining %d/%d channels" %(npch, nch)
        

        """Plotting starts"""

        #plot onedim in blue dashed
        #else plot in red
        if meth == 'onedim':
            prof = 'b--'
            lcol='b'
        else:
            prof = 'r-'
            lcol ='r'

        """1. PLOT PROFILES"""
        dimx, dimy = 3., 3.
        numsubplots = dimx*dimy
        numplots = int(np.ceil(npch/numsubplots))
        print "Num profile plots:", numplots
        """Compute residuals"""


       #"""Plot 1: Pulse profiles and fits"""
    
        if npch > 0:
            resdata = data_highsnr - model_highsnr
            resnormed = (resdata-resdata.mean())/resdata.std()
            
            if taussec_highsnr[0] > 1:
                taulabel =  taussec_highsnr
                taulabelerr = lmfitstdssec_highsnr
                taustring = 'sec'
            else:
                taulabel = taussec_highsnr*1000
                taulabelerr = lmfitstdssec_highsnr*1000
                taustring = 'ms'

            for k in range(numplots):
                j = int(numsubplots*k)
                figg = plt.figure(k+1,figsize=(int(4*dimx),int(3*dimy)))
                plots_remaining = int(npch - numsubplots*k)
                #print "Plots remaining", plots_remaining 
                for i in range(np.min([int(numsubplots),int(plots_remaining)])):
                    figg.subplots_adjust(left = 0.08, right = 0.98, wspace=0.35,hspace=0.35,bottom=0.15)
                    #plt.rc('text', usetex=True)
                    plt.rc('font', family='serif')              
                    plt.subplot(dimx,dimy,i+1)
                    plt.plot(profilexaxis,data_highsnr[j+i],alpha = 0.20)    
                    plt.plot(profilexaxis,model_highsnr[j+i],lw = 2.0, alpha = 0.85, label=r'$\tau: %.2f \pm%.2f$ %s' %(taulabel[j+i], taulabelerr[j+i], taustring))
                    plt.title('%s at %.1f MHz' %(pulsar, freqMHz_highsnr[j+i]))
                    plt.ylim(ymax=1.3*np.max(data_highsnr[j+i]))
                    plt.xlim(xmax=pulseperiod)
                    plt.xticks(fontsize=11)
                    plt.yticks(fontsize=11)
                    plt.xlabel('time (s)',fontsize=11)
                    plt.legend(fontsize=11,numpoints=1)
                    plt.ylabel('normalized intensity',fontsize=11)
                    plt.tight_layout()
                
                if savefigure == True:
                    figname = '%s_%s_%s_%d.png' %(os.path.basename(filepath),'fitted_profiles', meth, k)
                    plt.savefig(figname, dpi=200)
                    print "Saved figure %s in ./" %figname
                    if noshow == False:
                        plt.show()

            if verboseTag:
                for i in range(npch):
                    print "Channel %d" %i
                    print'Tau (ms): %.2f' %(1000*taussec_highsnr[i])
                    tau1GHz = tauatfreq(freqMHz_highsnr[i]/1000.,taussec_highsnr[i],1.0,4)
                    print 'tau1GHz_alpha_4 (ms) ~ %.4f' %(tau1GHz*1000)

            lmfitstdssec_highsnr = lmfitstdssec_highsnr[np.nonzero(lmfitstdssec_highsnr)]
            taussec_highsnr = taussec_highsnr[np.nonzero(lmfitstdssec_highsnr)]
            freqMHz_highsnr = freqMHz_highsnr[np.nonzero(lmfitstdssec_highsnr)]
            
            """Plot 2: Plot Gaussian fitting parameters and DM if selected"""
        
            if plotparams==True:
                print "\nPlotting Gaussian fit parameters w.r.t frequency\n"
                """Set plotting parameters"""
                alfval = 0.6
                markr= '*'
                msize=12
                plt.figure(numplots+1, figsize=(12,8))
                plt.subplots_adjust(left = 0.055, right=0.98,wspace=0.35,hspace=0.4,bottom=0.08)               
                """Fit models to sigma"""
                powmod = PowerLawModel()
                powpars = powmod.guess(paramset_highsnr[0], x=freqMHz_highsnr)
                powout = powmod.fit(paramset_highsnr[0], powpars, x=freqMHz_highsnr, weights=1/((paramsetstd_highsnr[0])**2))

                linmod = LinearModel()
                
                if len(freqMHz_highsnr) < 3:
                    raise RuntimeError("plotparams == True: Less than three frequency channels. Cannot compute quadratic or exponential fit for width evolution. Consider lowering snr_cut.")
                    
                else:
                    quadmod = QuadraticModel()          
                    quadpars = quadmod.guess(paramset_highsnr[0], x=freqMHz_highsnr)
                    quadout  = quadmod.fit(paramset_highsnr[0], quadpars, x=freqMHz_highsnr, weights=1/((paramsetstd_highsnr[0])**2))

                    expmod = ExponentialModel()
                    exppars = expmod.guess(paramset_highsnr[0], x=freqMHz_highsnr)
                    expout = expmod.fit(paramset_highsnr[0], exppars, x=freqMHz_highsnr, weights=1/((paramsetstd_highsnr[0])**2))


                """Fit a DM model to delta mu"""
                delnuarray = [-(1/freqMHz_highsnr[-1]**2-1/freqMHz_highsnr[i]**2) for i in range(npch)] ##in MHz
                delmuarray = [(paramset_highsnr[1][-1] - paramset_highsnr[1][i])*pbs for i in range(npch)] ##in seconds
                delmu_stdarray = [(paramsetstd_highsnr[1][-1] - paramsetstd_highsnr[1][i])*pbs for i in range(npch)]

                DM_linpars = linmod.guess(delmuarray, x=delnuarray)
                DM_linout  = linmod.fit(delmuarray, DM_linpars, x=delnuarray)

                DM_CCval = DM_linout.best_values['slope']
                DM_CCvalstd = DM_linout.params['slope'].stderr

                DMmodelfit = DM_linout.best_fit ##model gives deltime in seconds (used to shift data)

                DMconstant = 4148.808
                #uncertainty in the constant is 0.003 - only affects the Delta DM value in the 9th decimal
                DMval = (DM_CCval/DMconstant)
                DMvalstd = (DM_CCvalstd/DMconstant)
                #DMcheck = psr.DM_checker(freqmsMHz,bestpT_highSNR[1]*pbs)
                
                
                ## Plot reduced chi square:
                
                plt.subplot(2,3,1)
                plt.plot(freqMHz_highsnr, redchis_highsnr/np.power(comp_rmss_highsnr,2), markr,alpha=alfval,markersize = msize)
                plt.title(r'Reduced $\chi^2$ values', fontsize=12)
                plt.yticks(fontsize=12)
                plt.xticks(fontsize=12)
                plt.xlabel(r'$\nu$ MHz',fontsize=12)
                plt.ylabel(r'$\chi^2$',fontsize=12)
                
                ## Plot sigma:
                
                plt.subplot(2,3,2)
                #plt.errorbar(freqMHz_highsnr,paramset_highsnr[0]*pbs)
                plt.errorbar(freqMHz_highsnr,paramset_highsnr[0]*pbs,yerr =paramsetstd_highsnr[0]*pbs, fmt = markr,markersize=msize,capthick=2,linewidth=1.5,alpha=alfval)
                plt.plot(freqMHz_highsnr,powout.best_fit*pbs,'-', alpha=alfval,label='pow = %.2f' %powout.best_values['exponent'])
                plt.plot(freqMHz_highsnr,quadout.best_fit*pbs,'-',alpha=alfval, label='quad: a,b = %.3f,%.3f' %(quadout.best_values['a'],quadout.best_values['b']))
                plt.ylabel(r'$\sigma$ (sec)')
                plt.title(r'Width evolution', fontsize=12)
                plt.yticks(fontsize=12)
                plt.xticks(fontsize=12)
                plt.xlabel(r'$\nu$ MHz',fontsize=14)
                plt.legend(fontsize = 10, loc='best')
                plt.ylim(np.min(paramset_highsnr[0]*pbs)-1.5*np.median(paramsetstd_highsnr[0]*pbs),np.max(paramset_highsnr[0]*pbs)+1.5*np.median(paramsetstd_highsnr[0]*pbs))
                
                 ## Plot mean:
                
                plt.subplot(2,3,3)
                #plt.errorbar(freqMHz_highsnr,paramset_highsnr[1]*pbs)
                plt.errorbar(freqMHz_highsnr,paramset_highsnr[1]*pbs,yerr =paramsetstd_highsnr[1]*pbs, fmt = markr,markersize=msize,capthick=2,linewidth=1.5,alpha=alfval)
                plt.ylabel(r'$\mu$ (sec)')
                plt.title(r'Centroid evolution', fontsize=12)
                plt.yticks(fontsize=12)
                plt.xticks(fontsize=12)
                plt.xlabel(r'$\nu$ MHz',fontsize=14)
                plt.ylim(np.min(paramset_highsnr[1]*pbs)-1.5*np.median(paramsetstd_highsnr[1]*pbs),np.max(paramset_highsnr[1]*pbs)+1.5*np.median(paramsetstd_highsnr[1]*pbs))
                #plt.legend(fontsize = 9, loc='best')
                
                ## Plot amplitude:
                plt.subplot(2,3,4)
                #plt.errorbar(freqMHz_highsnr,paramset_highsnr[2]*pbs)
                plt.errorbar(freqMHz_highsnr,paramset_highsnr[2]*pbs,yerr =paramsetstd_highsnr[2]*pbs, fmt = markr,markersize=msize,capthick=2,linewidth=1.5,alpha=alfval)
                plt.ylabel(r'$\mu$ (sec)')
                plt.title(r'Amplitude evolution', fontsize=12)
                plt.yticks(fontsize=12)
                plt.xticks(fontsize=12)
                plt.xlabel(r'$\nu$ MHz',fontsize=14)
                plt.ylim(np.min(paramset_highsnr[2]*pbs)-1.5*np.median(paramsetstd_highsnr[2]*pbs),np.max(paramset_highsnr[2]*pbs)+1.5*np.median(paramsetstd_highsnr[2]*pbs))
                #plt.legend(fontsize = 9, loc='best')
                
                ## Plot DC:
                
                plt.subplot(2,3,5)
                #plt.errorbar(freqMHz_highsnr,paramset_highsnr[3]*pbs)
                plt.errorbar(freqMHz_highsnr,paramset_highsnr[3]*pbs,yerr =paramsetstd_highsnr[3]*pbs, fmt = markr,markersize=msize,capthick=2,linewidth=1.5,alpha=alfval)
                plt.ylabel(r'$\mu$ (sec)')
                plt.title(r'DC offset', fontsize=12)
                plt.yticks(fontsize=12)
                plt.xticks(fontsize=12)
                plt.xlabel(r'$\nu$ MHz',fontsize=14)
                plt.ylim(np.min(paramset_highsnr[3]*pbs)-1.5*np.median(paramsetstd_highsnr[3]*pbs),np.max(paramset_highsnr[3]*pbs)+1.5*np.median(paramsetstd_highsnr[3]*pbs))
                #plt.legend(fontsize = 9, loc='best')
                
                 ## Plot DM:
                plt.subplot(2,3,6)
                plt.errorbar(delmuarray,freqMHz_highsnr, fmt = markr, xerr=delmu_stdarray, alpha = alfval, markersize=msize)
                plt.plot(DMmodelfit,freqMHz_highsnr, '-', label=r'DM: $%.3f \pm %.3f$ $\rm{pc.cm}^{-3}$' %(DMval,DMvalstd), alpha = alfval)
                plt.xlabel(r'$\Delta \mu$ (sec)', fontsize =12)
                plt.yticks(fontsize=12)
                plt.xticks(fontsize=12)
                plt.title('Delta DM', fontsize=12)
                plt.ylabel(r'$\nu$ (MHz)',fontsize=14)
                plt.ticklabel_format(style='sci', axis='x',scilimits=(0,0))
                plt.legend(fontsize = 10, loc='best')
                plt.xlim(np.min(delmuarray)-5.0*np.median(delmu_stdarray),np.max(delmuarray)+5.0*np.median(delmu_stdarray))
                plt.tight_layout()

                if savefigure == True:
                    figname2 = '%s_%s.png' %(os.path.basename(filepath),'fitting_parameters')
                    plt.savefig(figname2, dpi=200)
                    print "Saved figure %s in ./" %figname2
                if noshow == False: 
                    plt.show()

            if plotflux == True:  ##Flux section needs debugging
                    ls = 'solid'
                    """Plot flux, and corrected flux spectrum"""
                    """Create unscattered profiles, i.e. Guassians"""   
                    
                    bins, profiles = [],[makeprofile(nbins = nbins, ncomps = 1, amps = paramset_highsnr[2][j], means = paramset_highsnr[1][j], sigmas = paramset_highsnr[0][j]) for j in range(npch)]
                    
                    unscatflux = []
                    for i in range(npch):
                        unscatfl = np.sum(profiles[j])/nbins
                        unscatflux.append(unscatfl)
                        
                     #smootheddata = smooth(data_highsnr[j],int(0.05*nbins))     
                    scatflux = [find_modelflux(model_highsnr[i],nbins) for i in range(npch)]
                    
                    climbvals = []
                    for i in range(npch):
                        climb = returnclimb(np.linspace(1,nbins,nbins),paramset_highsnr[1][i],paramset_highsnr[0][i],paramset_highsnr[2][i],taussec_highsnr[i],paramset_highsnr[3][i],nbins)
                        climbvals.append(climb)
                    
                    correctedflux = np.array([scatflux[i] + climbvals[i] for i in range(npch)])
                    print scatflux
                    print climbvals
                    
                    #per bin
                    meancorflux = np.mean(correctedflux)
                    meancorfluxerr = np.sqrt(np.sum(correctedflux**2))/len(correctedflux)

                    """Calculate error in Flux"""

                    sigmaWIDTH = paramsetstd_highsnr[0]*pbs #in seconds
                    sigmaAMP = paramsetstd_highsnr[2]  #in mJy
                    WIDTHS =paramset_highsnr[0]*pbs #in seconds
                    AMPS =paramset_highsnr[2] #in mJy

                    Expr1 = np.sqrt(2*np.pi)*AMPS
                    Expr2 = np.sqrt(WIDTHS)
                    AreaExpression = Expr1*Expr2

                    sigmaEx1 = np.sqrt(2*np.pi)*sigmaAMP
                    sigmaEx2 = Expr2*0.5*sigmaWIDTH/WIDTHS
                    
                    sigmaFlux =AreaExpression*np.sqrt(np.power(sigmaEx1/Expr1,2)+ np.power(sigmaEx2/Expr2,2))#+ 2*corsigA_highsnr*sigmaEx1*sigmaEx2/(Expr1*Expr2)

                    plt.figure(figsize=(10,6))
                    plt.plot(freqMHz_highsnr, correctedflux,'k-', linewidth=2.0)
                    plt.plot(freqMHz_highsnr, scatflux,'r--', linewidth=2.0)                
                    plt.fill_between(freqMHz_highsnr,scatflux,correctedflux, alpha=alfval,facecolor='r')
                    eb = plt.errorbar(freqMHz_highsnr,correctedflux,yerr=sigmaFlux, fmt=markr,markersize=10.0, alpha=1.0,capthick=2,linewidth=1.5)
                    eb[-1][0].set_linestyle(ls)
                    #plt.errorbar(freqMHz_highsnr, unscatflux,yerr=sigmaFlux, fmt=markr,markersize=10.0, alpha=alfval)
                    plt.title('Flux Spectrum', fontsize=12)
                    plt.xticks(fontsize=12)
                    plt.yticks(fontsize=12)
                    plt.xlabel(r'$\nu$ (MHz)',fontsize=12)
                    plt.ylabel(r'Calibrated flux (mJy)',fontsize=12)
                
        
        return freqMHz_highsnr, taussec_highsnr, lmfitstdssec_highsnr


def produce_tauspectrum(freqMHz,taus,tauserr,log=True, plotspecevo=False,
        savefigure=False, noshow=False):
    
    plt.figure(figsize=(10,6))
    if noshow == False:
        plt.show()
    if len(taus) <= 2:
        print "Can't perform power law fit on single value. Returning spectral index = 0.0"
        alpha=0
        alphaerr=0
        fit=taus

    
    else:  
        freqGHz = freqMHz/1000.
        powmod = PowerLawModel()
        powparstau = powmod.guess(taus,x=freqGHz)

        #remove tauserr = 0 entries

        tauserr = tauserr[np.nonzero(tauserr)]
        taus = taus[np.nonzero(tauserr)]
        freqGHz = freqGHz[np.nonzero(tauserr)]
        freqMHz = freqMHz[np.nonzero(tauserr)]
        
        if len(taus) < 3:
            print "Two or fewer tau-values. Cannot compute tau-spectrum."
            alpha=0
            alphaerr=0
            fit=taus
        else:
            powout = powmod.fit(taus,powparstau,x=freqGHz,weights=1.0/(np.power(tauserr,2)))

            print "\nPlotting fitted tau-spectrum\n"
            print(fit_report(powout.params))
            fit = powout.best_fit
            alpha = -powout.best_values['exponent']
            alphaerr = powout.params['exponent'].stderr

     
            plt.errorbar(freqMHz,taus,yerr=tauserr,fmt='*-', markersize=10.0,capthick=2,linewidth=1.5, label='data')
            plt.plot(freqMHz,fit,ls='dashed',linewidth=1.5,label=r'$\alpha = %.1f \pm %.1f$' %(alpha,alphaerr))
            if log:
                plt.xscale('log')
                plt.yscale('log')
            plt.yticks(fontsize=12)
            #ticksMHz = (freqMHz).astype(np.int)[0:len(freqMHz):2]
            #plt.xticks(ticksMHz,ticksMHz,fontsize=11)
            plt.xticks(fontsize=11)
            leg = plt.legend(fontsize=12,numpoints=None)
            plt.xlabel(r'$\nu$ (MHz)',fontsize=14, labelpad=15.0)
            plt.ylabel(r'$\tau$ (sec)',fontsize=14)
            plt.xlim(xmin = 0.95*freqMHz[0],xmax=1.05*freqMHz[-1])
            plt.gcf().subplots_adjust(bottom=0.15)

            if savefigure == True:
                    figname = '%s_%s.png' %(os.path.basename(filepath),'tau_spectrum')
                    plt.savefig(figname, dpi=200)
                    print "Saved figure %s in ./" %figname
            if noshow == False: 
                plt.show()

            if plotspecevo == True:
                """Calculate spectral index evolution with addition of lowwer freq. channels"""
                npch = len(freqMHz)
                spec_sec = np.zeros(npch-1)
                spec_std_sec = np.zeros(npch-1)

                for i in range(npch-1):
                    bb = npch - (i+2)  ##index begin
                    ee = npch ## index end
                #   print freqms[bb:ee]

                    powparstau_sec = powmod.guess(taus[bb:ee],x=freqMHz[bb:ee])
                    powouttau_sec = powmod.fit(taus[bb:ee],powparstau_sec, x=freqMHz[bb:ee],weights=1/(tauserr[bb:ee]**2))
                    spec_sec[i] = -powouttau_sec.best_values['exponent']
                    spec_std_sec[i] = powouttau_sec.params['exponent'].stderr

                freq_incl = freqMHz[::-1][1:]

                plt.figure(figsize=(10,6))
                plt.errorbar(freq_incl,spec_sec,yerr=spec_std_sec, fmt='*',alpha=0.6,markersize=12.0,capthick=2,linewidth=0.5)                    
                plt.title(r'Change in spectral index', fontsize=12)
                for x,y in zip(freq_incl[0::2], spec_sec[0::2]):
                    plt.annotate('%.2f' %y, xy=(x,1.08*y), xycoords='data',textcoords='data')
                plt.ylim(plt.ylim()[0],1.1*plt.ylim()[1])
                plt.yticks(fontsize=12)
                plt.xticks(fontsize=12)
                plt.xlabel(r'Lowest freq included (MHz)',fontsize=12)
                plt.ylabel(r'$\alpha$',fontsize=12)

                if savefigure == True:
                    figname = '%s_%s.png' %(os.path.basename(filepath),'alpha_evolution')
                    plt.savefig(figname, dpi=200)
                    print "Saved figure %s in ./" %figname
                if noshow == False:
                    plt.show()

    return freqMHz, alpha, alphaerr, fit

        


def produce_tauspectrum_highHBA(freqMHz,taus,tauserr,freqMHzhigh,tauhigh,tauerrhigh):
    freqGHz = freqMHz/1000.
    powmod = PowerLawModel()
    powparstau = powmod.guess(taus,x=freqGHz)
    
    #remove tauserr = 0 entries
    
    tauserr = tauserr[np.nonzero(tauserr)]
    taus = taus[np.nonzero(tauserr)]
    freqGHz = freqGHz[np.nonzero(tauserr)]
    freqMHz = freqMHz[np.nonzero(tauserr)]
    
    powout = powmod.fit(taus,powparstau,x=freqGHz,weights=1/(np.power(tauserr,2)))
    
    print(fit_report(powout.params))
    fit = powout.best_fit
    
    alpha = -powout.best_values['exponent']
    alphaerr = powout.params['exponent'].stderr
    amp = powout.params['amplitude']
    
    allfreq = np.append(freqMHz,freqMHzhigh)
    fithigh = amp*np.power(allfreq/1000, -alpha)

    fig = plt.figure(figsize=(12,6))         
    plt.errorbar(freqMHz,taus,yerr=tauserr,fmt='*-',markersize=10.0,capthick=2,linewidth=1.5,label=r'$\alpha = %.1f \pm %.1f$'%(alpha,alphaerr))
    plt.errorbar(freqMHzhigh,tauhigh,yerr=tauerrhigh,fmt='o')
    plt.plot(freqMHz,fit,'k--',linewidth=1.5)
    plt.xscale('log')
    plt.yscale('log')
    plt.yticks(fontsize=12)
    ticksMHz = (freqMHz).astype(np.int)[0:len(freqMHz):2]
    plt.xticks(ticksMHz,ticksMHz,fontsize=12)
    plt.legend(fontsize=12,numpoints=None)
    plt.xlabel(r'$\nu$ (MHz)',fontsize=14, labelpad=15.0)
    plt.ylabel(r'$\tau$ (sec)',fontsize=14)
    plt.xlim(xmin = 0.95*freqMHz[0],xmax=1.05*freqMHzhigh[-1])
    plt.gcf().subplots_adjust(bottom=0.15)
    
    return freqMHZ, alpha, alphaerr, fit, fithigh


def produce_tauspectrum_noplot(freqMHz,taus,tauserr):
    freqGHz = freqMHz/1000.
    powmod = PowerLawModel()
    powparstau = powmod.guess(taus,x=freqGHz)
    
    tauserr = tauserr[np.nonzero(tauserr)]
    taus = taus[np.nonzero(tauserr)]
    freqGHz = freqGHz[np.nonzero(tauserr)]
    freqMHz = freqMHz[np.nonzero(tauserr)]
    
    powout = powmod.fit(taus,powparstau,x=freqGHz,weights=1/(np.power(tauserr,2)))
    
    print(fit_report(powout.params))
    fit = powout.best_fit
    alpha = -powout.best_values['exponent']
    alphaerr = powout.params['exponent'].stderr

    return freqMHz, alpha, alphaerr, fit




def read_Taufitter(filename):
    f = open(filename)
    lines = f.readlines()
    header0 = lines[0]
    h0_lines = header0.split()
    pulsar = h0_lines[2]
    nch = int(h0_lines[4])
    nbins = int(h0_lines[6])
    data = np.loadtxt(filename,skiprows=2)
    if 'onedim' in filename:
        meth = 'onedim'
    elif 'iso' in filename:
        meth = 'iso'
    else:
        meth = 'unknown'
    return pulsar, nch, meth, data


if __name__ == '__main__':
    import seaborn as sns
    sns.set_palette("muted")
    
    """Define options to the script"""
    parser = argparse.ArgumentParser()
    parser.add_argument('-f','--filename', help="Required. Provide the pathway to the ascii data files (created using pdv -KAt).")
    parser.add_argument('-p','--period',type=float, help="Optional. Provide the period of the pulsar in seconds. Default will call psrcat for P0.")
    parser.add_argument('-m','--method', help="Choosing method to fit data or simulation. Choose between 'onedim' and 'iso'. Default is 'iso'.")
    parser.add_argument('-snrcut','--snr_threshold', help="Optional. S/N threshold cut-off. Average per channel profiles with estimated S/N values lower than this will be excluded. Default value is None.")
    parser.add_argument('-noshow','--no_show', default=False, help="Optional. Default value is False (i.e does display figures). Use -noshow to suppress all plt.show() commands.", action='store_true')
    parser.add_argument('-nosavefig','--no_save_figure', default=False, help="Optional. Default is False (i.e. saves figures). To suppress all plt.savefig() commands use -nosavefig", action='store_true')
    args = parser.parse_args()

    """Allocate variable names to the parsed options"""
    filepath = args.filename
    if filepath == None:
        raise RuntimeError('No filename specified. Please use -f option.')

    pulseperiod = args.period
    meth = args.method
    snrcut = args.snr_threshold
    if snrcut:
        snrcut = float(snrcut)

    noshow = args.no_show
    no_savefig = args.no_save_figure
    savefig = not no_savefig

    """Read ascii file header in full"""
    pulsar, nch, nbins,nsub, lm_rms, tsub = read_headerfull(filepath)
    
    if meth == None:
        print "No fitting-method chosen. Will default to an isotropic fitting model. \n Use option -m with 'onedim' to change."
        meth = 'iso'
        
    if pulseperiod == None:
        print "Calling prscat subprocess to obtain pulsar period"
        command = ['psrcat','-nohead', '-o', 'short', '-c', 'JNAME p0 DM', pulsar]
        pulsar_psrcat = subprocess.Popen(command, shell=False, cwd='.', stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        (stdoutdata, stderrdata) = pulsar_psrcat.communicate()
        return_code = pulsar_psrcat.returncode

        if return_code != 0:
            raise RuntimeError('Problems with psrcat: %s') % stdoutdata
        elif ' not in catalogue' in stdoutdata:
            raise RuntimeError('Pulsar not found in psrcat.')

        pname = stdoutdata.split()[1]
        per= float(stdoutdata.split()[2])
        dm = float(stdoutdata.split()[3])
        
        if pname != pulsar:
            raise RuntimeError("Pulsar name mismatch")
        print "%s. Period: %.3f sec, DM: %.3f pc cm^-3"   %(pname, per,dm)

    else:
        per=pulseperiod

    """Produce scattering fits to pulse profiles"""
    freqMHz, taussec, taustdssec = produce_taufits(filepath, meth=meth,
            pulseperiod=per, snr_cut = snrcut, verbose=False, plotparams=True,
            savefigure=savefig)
    
    """Produce tauspectrum"""
    freqMHz, alpha, alphaerr, fit = produce_tauspectrum(freqMHz, taussec,
            taustdssec, log=True, plotspecevo=True, savefigure=savefig)

