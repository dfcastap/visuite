# -*- coding: utf-8 -*-
"""
Created on Wed Sep  9 13:39:55 2015

@author: diego
"""

import numpy as np
import pylab as plt
import glob

ell = 0
colors = ["red","blue","green","k"]
if ell%2==1:
    par = "ODD"
else:
    par ="EVEN"


loc = "/home/diego/Documents/ROTORCmodels/visibilities/2p5Msun/V0"
locbase = "/home/diego/Documents/ROTORCmodels/visibilities/2p5Msun/"
#fmag = glob.glob(loc+"/MODE_EVEN_"+str(51)+"/magnitudes*")

bfile = np.genfromtxt(locbase+"2p5Msun_V0_magnitudes_static")
base = bfile[2]
b_v_base = -2.5*np.log10(bfile[4]/bfile[5])

if ell==0:
    modes = [51, 69, 78, 85, 92] #l=0 M2p5 V0
    freqs = [1.58717,2.05922,2.49359,2.95717,3.46299] #l=0 M2p5 V0
    
if ell==1:
    modes = [5, 61] #l=1 M2p5 V0
    freqs = [0.78296,1.6348,2.16027,2.66899,3.18255,3.70941] #l=1 M2p5 V0
    
if ell==2:
    modes = [8,34,52,66,76,84,90] #l=2 M2p5 V0
    freqs = [0.83464,1.17430,1.60199,1.94807,2.36756,2.86427,3.38705,3.92415]

if ell==4:
    modes = [12,26,40,49,61,72,80,88]
    freqs = [0.85613,1.03514,1.28624,1.54099,1.79625,2.12890,2.64669,3.17963] # l=4 M2p5 V0
    
if ell==6:
    modes = [67,75,83,91]
    freqs = [1.97688,2.25056,2.82636,3.39865]

mag = []
b_v = []
for i in modes:
    fmag = glob.glob(loc+"/MODE_"+par+"_"+str(i)+"/magnitudes*")
    print fmag
    mag.append(np.genfromtxt(fmag[0]))

for i in range(len(modes)):
    #plt.plot([freqs[i],0],[freqs[i],mag[i][0,2]/base])
    #plt.plot([[freqs[i],0],[freqs[i],1]])
    c = colors[ell/2]
    lstyle = "-"
    if par == "ODD":
        c = colors[ell-1]
        lstyle = "--"
        
    b_v = -2.5*np.log10(mag[i][4]/mag[i][5])
    #plot1 = plt.vlines(freqs[i], 0, np.abs(mag[i][2]-base)/base,color=c,linestyle=lstyle)
    plot4 = plt.vlines(freqs[i], 0, (b_v-b_v_base)/b_v_base,color=c,linestyle=lstyle)
    
#plt.legend([plot1,plot2,plot3],["$\ell = 0$","$\ell = 2$","$\ell = 4$"],loc="best")
#plt.ylim(0,1.5)
plt.title("Mode luminosities M=2.5 V=0")
plt.ylabel("Normalized L")
plt.xlabel("Frequency")
plt.show()