# -*- coding: utf-8 -*-
"""
Created on Wed Sep  9 13:39:55 2015

@author: diego
"""

import numpy as np
import pylab as plt
import glob

ell = 4
colors = ["red","blue","green"]

loc = "/home/diego/Documents/ROTORCmodels/visibilities/2p5Msun/V0"
fmag = glob.glob(loc+"/MODE_EVEN_"+str(51)+"/magnitudes*")
bfile = np.genfromtxt(fmag[0])
base = bfile[0,2]

if ell==0:
    modes = [51, 69, 78, 85] #l=0 M2p5 V0
    freqs = [1.58717,2.05922,2.49359,2.95717] #l=0 M2p5 V0

if ell==2:
    modes = [8,34,52,66,76] #l=2 M2p5 V0
    freqs = [0.83464,1.17430,1.60199,1.94807,2.36756]

if ell==4:
    modes = [12,26,40,49,61,72,80]
    freqs = [0.85613,1.03514,1.28624,1.54099,1.79625,2.12890,2.64669] # l=4 M2p5 V0

mag = []
for i in modes:
    fmag = glob.glob(loc+"/MODE_EVEN_"+str(i)+"/magnitudes*")
    print fmag
    mag.append(np.genfromtxt(fmag[0]))
    

for i in range(len(modes)):
    #plt.plot([freqs[i],0],[freqs[i],mag[i][0,2]/base])
    #plt.plot([[freqs[i],0],[freqs[i],1]])
    plot3 = plt.vlines(freqs[i], 0, mag[i][-1,2]/base,color=colors[ell/2])
    
#plt.legend([plot1,plot2,plot3],["$\ell = 0$","$\ell = 2$","$\ell = 4$"],loc="best")
plt.ylim(0,1.5)
plt.title("Mode luminosities M=2.5 V=0")
plt.ylabel("Normalized L")
plt.xlabel("Frequency")
plt.show()