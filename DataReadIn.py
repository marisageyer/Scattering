# -*- coding: utf-8 -*-
"""
Created on Tue Jan 26 09:11:24 2016

@author: marisa
"""

import numpy as np
import os

def read_header(filepath):
    f = open(filepath)
    lines = f.readlines()
    header0 = lines[0]
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
#    return file_name, pulsar_name, nsub, nch, npol, nbins, rms
    return pulsar_name, nch, nbins, rms

def read_headerfull(filepath):
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
    
    
def read_data(filepath, profilenumber, nbins):
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
    

def read_timestringdata(filepath, profilenumber, nbins):
    #single frequency scrunched, time string data
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
                
    return np.array(datalist), freqm


def read_fittxtfiles(filepath,filename):
    pulsar_name = filename[0:10]
    if 'cycle5' in filename:    
        datac = 'cycle'
    elif 'census' in filename:
        datac = 'census'
    elif 'comm' in filename:
        datac = 'comm'
    else:
        datac = 'unknown'
    if 'onedim' in filename:
        meth = 'onedim'
    elif 'iso' in filename: 
        meth = 'iso'
    else:
        meth = 'unknown'
    loadedarray = np.loadtxt(os.path.join(filepath,filename))
    return pulsar_name, loadedarray, datac, meth


def read_raw(filename):
    pulsar_name = filename[0:10]
    loaded = np.loadtxt(filename)
    nch = np.shape(loaded)[0]
    nbins = np.shape(loaded)[1] -1
    datas = loaded[:,0:nbins]
    freq = loaded[:,nbins]
    return pulsar_name, nch, nbins, datas, freq
    
def read_Krishnak(filename):
    #function for reading data downloaded from: http://rac.ncra.tifr.res.in/da/pulsar/pulsar.html
    pulsar_name = filename[-17:-5]
    loaded = np.loadtxt(filename)
    nch = 1
    nbins = np.shape(loaded)[0]
    datas = loaded[:,1]
    freq = 0.327 #(327MHz)
    return pulsar_name, nch, nbins, datas, freq    

def log_writer(logpath,logfilename, my_string):
    logplace = os.path.join(logpath,logfilename)
    with open(logplace,'a') as log_file:
        log_file.write(my_string + '\n')

def log_reader(logpath,logfile):
    logplace = os.path.join(logpath,logfile)
    f = open(logplace)
    lines = f.readlines()
    nbins = lines[2].split()[3]
    nbins = int(nbins)
    return lines, nbins
    

def read_Stefan_data(filepath, profilenumber, nbins):
    d = open(filepath)
    lines = d.readlines()
    
    profile_start = profilenumber*(nbins+1)
    profile_end = profile_start + nbins
    
    lines_block = lines[profile_start:profile_end]
    
    datalist = []
    for i in range(nbins):
        data= float(lines_block[i].split()[3])
        datalist.append(data)
    return np.array(datalist)


def read_Paul(filename):
    #function for reading /Paul/1937+21_mjd_57397.txt at the moment
    #Single column of ascii data
    pulsar_name = "B1937+21"
    datas = np.loadtxt(filename)
    nch = 1
    nbins = np.shape(datas)[0]
    mainpulse = datas[nbins/2:nbins]
    nbinsmain = nbins/2
    freq = 1.4 #(1.4GHz)
    return pulsar_name, nch, nbinsmain, mainpulse, freq  
    

def skip_first_col(fname, delimiter=None):
    with open(fname, 'r') as fin:
        for line in fin:
            try:
               yield line.split(delimiter, 1)[1]
            except IndexError:
               continue

           
#def read_Stefan(filepath):
#    datamatrix = np.loadtxt(filepath)
#    ##determine number of blocks
#    number_blocks = datamatrix[-1][0]
#    blocks = []
#    for i in range(number_blocks):
#      block = datamatrix[datamatrix[:, 0] == i, :]
#        blocks.append(block)
#    return blocks[0:10]
#
    ##divide observation up into sub_int blocks


    ##divide each block into channels

    
