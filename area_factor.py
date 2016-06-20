# -*- coding: utf-8 -*-
"""
Created on Sat May  7 19:57:09 2016

@author: castaned
"""
import main_modules as mmod
import pylab as plt
import seaborn as sns
sns.set(style="white",rc={"figure.figsize": (8, 8),'axes.labelsize': 16,
                              'ytick.labelsize': 12,'xtick.labelsize': 12,
                              'legend.fontsize': 16,'axes.titlesize':16,'font.size':14})

r_ord = 4
kind = "p"
ells = [0,1,2,3,4,5,6]
mds = mmod.list_gen(ells,r_ord,r_ord,kind)
mss = "1p875"

incl = [0,10,20,30,40,50,60,70,80,90]
for i in mds:
    mde = i
    ll = int(mde.split()[0])
    clr = sns.color_palette("Set2", 12)[ll]
    plt.plot(incl,map(lambda x: x,[mmod.emode(mss,0,mde).calc_area(incl=i) for i in incl]),label=r"$\ell$="+str(ll),color=clr)
    

#plt.ylim(0,1.2)
plt.xlabel("Inclination")
plt.ylabel("(1.0 - Area_Factor)")
plt.legend(loc="best")
plt.title(r"M="+mss+"M$\odot$ - V = 0 km s$^{-1}$")
plt.grid()