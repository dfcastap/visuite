# -*- coding: utf-8 -*-
"""
Created on Wed Sep  9 13:39:55 2015

@author: diego
"""

import numpy as np
import pylab as plt
import glob
import seaborn as sns
sns.set(style="white",rc={"figure.figsize": (8, 8),'axes.labelsize': 16,
                              'ytick.labelsize': 12,'xtick.labelsize': 12,
                              'legend.fontsize': 16,'axes.titlesize':18,'font.size':14})
wavebands = [3650.,4450.,5510.,6580.,8060.]

wal_wave = [3250,3630,3840,4320,5470]                        

ell = 0
colors = ["red","blue","green","magenta","purple"]
if ell%2==1:
    par = "ODD"
else:
    par ="EVEN"


incl = [0,15,30,45,60,75,90.0]

loc = "/home/diego/Documents/ROTORCmodels/visibilities/2p5Msun/V0"
locbase = "/home/diego/Documents/ROTORCmodels/visibilities/2p5Msun/"
#fmag = glob.glob(loc+"/MODE_EVEN_"+str(51)+"/magnitudes*")

bfile = np.genfromtxt(locbase+"mags_walraven_2p5Msun_V0_static")
base = bfile[2]
b_v_base = -2.5*np.log10(bfile[4]/bfile[5])

if ell==0:
    modes = [51, 69, 78, 85, 92, 100, 106, 114, 122, 129] #l=0 M2p5 V0
    freqs = [1.58717,2.05922,2.49359,2.95717,3.46299,3.99529,4.54267,5.09092,5.64618] #l=0 M2p5 V0
    
if ell==1:
    modes = [5, 61, 81, 89, 96, 104, 111, 133] #l=1 M2p5 V0
    freqs = [0.78296,1.63480,2.16027,2.66899,3.18255,3.70941] #l=1 M2p5 V0
    
if ell==2:
    modes = [8, 34, 52, 66, 76, 84, 90, 97, 105, 113, 128] #l=2 M2p5 V0
    freqs = [0.83464,1.17430,1.60199,1.94807,2.36756,2.86427,3.38705,3.92415,4.47259,5.01979]
    
if ell==3:
    modes = [34, 52, 64, 78, 87, 94, 102, 115, 139]
    freqs = [1.08830,1.42856,1.68216,2.06952,2.54134,3.04737,3.57597] #l=3 M2p5 V=0

if ell==4:
    modes = [12,26,40,49,61,72,80,88]
    freqs = [0.85613,1.03514,1.28624,1.54099,1.79625,2.12890,2.64669,3.17963] # l=4 M2p5 V0
    
if ell==6:
    modes = [67,75,83,91]
    freqs = [1.97688,2.25056,2.82636,3.39865]

modes = [129]
freqs = [1.63480]
mag = []
b_v = []
for i in modes:
    for j in incl:
        fmag = glob.glob(loc+"/MODE_"+par+"_"+str(i)+"/walraven_magnitudes"+"_i"+str(j)+"_MODE_"+str(i))
        #print fmag
        mag.append(np.genfromtxt(fmag[0]))

def plot_amplitude_ratios():
    global ell,modes,colors,bfile,base,par
    for i in range(len(modes)):
        #plt.plot([freqs[i],0],[freqs[i],mag[i][0,2]/base])
        #plt.plot([[freqs[i],0],[freqs[i],1]])
        c = colors[ell/2]
        lstyle = "-"
        if par == "ODD":
            c = colors[-ell]
            lstyle = "--"
            
        #b_v = -2.5*np.log10(mag[i][4]/mag[i][5])
        #plot5 = plt.vlines(freqs[i], 0, np.abs(mag[i][2]-base)/base,color=c,linestyle=lstyle)
        #print mag[i][2]
        #plot4 = plt.vlines(freqs[i], 0, (b_v-b_v_base)/b_v_base,color=c,linestyle=lstyle)
        #mag_ratios = mag[i][3:]/mag[i][3]
        mag_ratios = mag[i][3:]/bfile[3:]
        mag_ratios /= mag_ratios[0]
        plt.plot(wavebands,mag_ratios)
        plt.text(wavebands[-1]+140,mag_ratios[-1],r"$\ell$ = "+str(ell) + " - "+ str(i))
        plt.xlabel("Wavelength")
        plt.ylabel("Apmplitude ratio")


#plot_amplitude_ratios()

def lum_plot():
    global ell,modes,colors,bfile,base,par,mag,freqs
    for i in range(len(modes)):
        #plt.plot([freqs[i],0],[freqs[i],mag[i][0,2]/base])
        #plt.plot([[freqs[i],0],[freqs[i],1]])
        c = colors[ell/2]
        lstyle = "-"
        if par == "ODD":
            c = colors[-ell]
            lstyle = "--"
            
        #b_v = -2.5*np.log10(mag[i][4]/mag[i][5])
        #plot5 = plt.vlines(freqs[i], 0, np.abs(mag[i][2]-base)/base,color=c,linestyle=lstyle)
        #print mag[i][2]
        plt.vlines(freqs[i], 0, (mag[i][2]-base)/base,color=c,linestyle=lstyle)
        
#lum_plot()
#plt.legend([plot1,plot2,plot3],["$\ell = 0$","$\ell = 2$","$\ell = 4$"],loc="best")
#plt.ylim(0,1.5)
#plt.xlim(1,4)
#plt.title("Mode luminosities M=2.5 V=0")
#plt.ylabel("Normalized L")
#plt.xlabel("Frequency")
#plt.show()
        
def lum_incl_plot():
    global ell,modes,colors,bfile,base,par,mag,freqs,incl
    for i in range(len(modes)):
        #plt.plot([freqs[i],0],[freqs[i],mag[i][0,2]/base])
        #plt.plot([[freqs[i],0],[freqs[i],1]])
        tmp = []
        for j in range(len(incl)):
            tmp.append((mag[j][2]-base)/base)
        
        plt.plot(incl, np.abs(np.array(tmp)))
        
            
        #b_v = -2.5*np.log10(mag[i][4]/mag[i][5])
        #plot5 = plt.vlines(freqs[i], 0, np.abs(mag[i][2]-base)/base,color=c,linestyle=lstyle)
        #print mag[i][2]
        #plt.vlines(freqs[i], 0, (mag[i][2]-base)/base,color=c,linestyle=lstyle)

lum_incl_plot()
plt.grid()