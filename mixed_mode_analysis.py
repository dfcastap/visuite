# -*- coding: utf-8 -*-
"""
Created on Wed Apr 27 12:47:33 2016

@author: castaned
"""

import main_modules as mmod
import numpy as np
import pylab as plt
import seaborn as sns

sns.set(style="white",rc={"figure.figsize": (8, 8),'axes.labelsize': 16,
                              'ytick.labelsize': 12,'xtick.labelsize': 12,
                              'legend.fontsize': 16,'axes.titlesize':16,'font.size':14})
                              

mname = "0 p4"
#mname = "6 g3"

lst = ["-","--",":"]

inc = [7]
vels = [0,2,8,9]

clr = sns.color_palette("husl", len(vels))

mode = []                           
[mode.append(mmod.emode("1p875",i,mname,lpert=True)) for i in vels]

#[plt.plot(mode[i].xi_r[:,inc[j]],ls=lst[j],color=clr[i],label=r"V="+mode[i].eq_vel+"km s$^{-1}$") for j in range(len(inc)) for i in range(len(vels))]
#plt.grid()
#plt.ylabel(r"$\xi_{r}$")
#plt.xlabel("R [zone #]")
#plt.legend(loc="best")

for i in range(len(vels)):
    if i==9 or i==9:
        mode[i].plot_dr_r(invert=True)
    else:
        mode[i].plot_dr_r()
