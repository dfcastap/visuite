# -*- coding: utf-8 -*-
"""
Created on Wed Mar 16 01:10:49 2016

@author: castaned
"""

import main_modules as mmod
import numpy as np
import seaborn as sns
import pylab as plt

sns.set(style="white",rc={"figure.figsize": (10, 8),'axes.labelsize': 16,
                              'ytick.labelsize': 12,'xtick.labelsize': 12,
                              'legend.fontsize': 16,'axes.titlesize':16,'font.size':14})

mass = "1p875"
vel = 4
ells = [0,1,2,3,4,5,6]
kind = "p"
r_ord = 4
m_list = mmod.list_gen(ells,r_ord,r_ord,kind)
modes = []
phase = False

incl = np.linspace(0,90,10)
if not phase:
    lblphase = "0"
else:
    lblphase = "123"

for i in range(len(m_list)):
    modes.append(mmod.emode(mass,vel,m_list[i],lpert=True))

areas = []
mag_ampl = []
for i in range(len(modes)):
    areas.append(modes[i].area)
    print modes[i].name
    if phase:
        mag_ampl.append(np.abs(-2.5*np.log10(modes[i].w_phi_mags_min[:,4]/modes[i].w_phi_mags[:,4])))
    else:
        mag_ampl.append(np.abs(-2.5*np.log10(modes[i].w_mags_min[:,4]/modes[i].w_mags[:,4])))

mag_ampl_incl = []
order_by_mag_ampl = []
for i in range(len(mag_ampl[0])):
    tmp = []
    for j in range(len(modes)):
        tmp.append(mag_ampl[j][i])
    mag_ampl_incl.append(tmp)
    order_by_mag_ampl.append(sorted(range(len(tmp)), key=lambda k: tmp[k])[::-1])

mag_ampl_incl = np.transpose(np.array(mag_ampl_incl))
order_by_mag_ampl = np.transpose(np.array(order_by_mag_ampl))
areas = np.array(areas)
order_by_area = sorted(range(len(areas)), key=lambda k: areas[k])

#orders = np.zeros((len(modes),len(ells)))
#orders[:,0] = order_by_area
#orders[:,1::] = order_by_mag_ampl

ranked_mag_ampl = np.zeros(mag_ampl_incl.shape)
for i in range(len(ells)):
    ranked_mag_ampl[i,:] = [mag_ampl_incl[order_by_mag_ampl[i,j],j] for j in range(len(incl))]    

fig, ax = plt.subplots()
x = np.linspace(0,90,10)
y = np.linspace(0,ells[-1],len(ells))
z = 1.*(ranked_mag_ampl)
cax = plt.pcolor(x,y,z,cmap="cool")
cbar = fig.colorbar(cax,label="Mag. Amplitude [mag]")
plt.gca().invert_yaxis()
plt.ylabel("Ranking")
plt.xlabel("Inclination")
#plt.gca().get_xaxis().set_visible(False)
#plt.gca().get_yaxis().set_visible(False)
#plt.text(0,len(ells)-len(ells)*0.1,"Inclination ->")
#plt.text(-5,len(ells)-3,"Ranking ->",rotation="vertical")

clr = sns.light_palette("red", input="xkcd", n_colors=len(ells)+4,reverse=True)
clr = sns.light_palette("gray", reverse=True,n_colors=len(ells))
nells = len(ells)
for i in range(len(incl)):
    for j in range(len(y)):
        yloc = float(j)*(float(ells[-1])/len(y))+0.5
        plt.text(i*(9.)+3.5,yloc,str(ells[order_by_mag_ampl[j,i]]),color="k")
        
plt.title(r"M="+mass.replace("p",".")+"M$\odot$, V="+modes[0].eq_vel+"km s$^{-1}$, "+kind+" = "+str(r_ord)+" - L_phase = "+lblphase+"$^{\circ}$")

#pfol = "/home/castaned/Documents/ROTORCmodels/visibilities/plots/toBob_apr22/mag_ampl/"
#order_stack = (np.c_[order_by_mag_ampl,np.array(order_by_area)])
#hdr = "Mode Rankings. Cols: i=0,10,20,30,40,50,60,70,80,90,Area Ranking "
#np.savetxt(pfol+"/M"+mass+"/rankings_V"+modes[0].eq_vel+"_ell-"+str(ells[0])+"-"+str(ells[-1])+"_"+kind+str(r_ord)+"_Lphase_"+lblphase+".txt", order_stack,fmt="%3i", header=hdr)