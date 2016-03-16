# -*- coding: utf-8 -*-
"""
Created on Wed Mar 16 01:10:49 2016

@author: castaned
"""

import main_modules as mmod
import numpy as np
import seaborn
import pylab as plt


mass = "1p875"
vel = 0
m_list = mmod.list_gen([0,1,2,3,4,5,6,7,8,9,10],1,1,"p")
modes = []
for i in range(len(m_list)):
    modes.append(mmod.emode(mass,vel,m_list[i]))
    
areas = []
mag_ampl = []
for i in range(len(modes)):
    areas.append(modes[i].area)
    mag_ampl.append(np.abs(-2.5*np.log10(modes[i].w_mags_min[:,4]/modes[i].w_mags[:,4])))

mag_ampl_incl = []
order_by_mag_ampl = []
for i in range(len(mag_ampl[0])):
    tmp = []
    for j in range(len(modes)):
        tmp.append(mag_ampl[j][i])
    mag_ampl_incl.append(tmp)
    order_by_mag_ampl.append(sorted(range(len(tmp)), key=lambda k: tmp[k])[::-1])
    
    
areas = np.array(areas)
order_by_area = sorted(range(len(areas)), key=lambda k: areas[k])