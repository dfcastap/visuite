# -*- coding: utf-8 -*-
"""
Spyder Editor

This temporary script file is located here:
/home/diego/.spyder2/.temp.py
"""
import numpy as np
import pylab as pl
import scipy.interpolate as interp
from scipy.misc import derivative

"""
def legendre(data1):
    nn,mm = data1.shape;
    #choice = input(' Enter 0 for horizontal variation of radial perturbation, 1 of pressure perturbation  ')  ;
    choice = 0;
    if choice == 0:
        data = data1[:,0] ;
    else:
        data = data1[:,2] ;
    
   # step = 90. / (nn + 1) ;
    ang = np.linspace(10,80,8) ;
    ang_rad = np.deg2rad(np.transpose(ang));
    c2 = np.cos(ang_rad)**2 ;
    c4 = c2 * c2 ;
    c6 = c4 * c2 ;
    c8 = c6 * c2 ;
    c10 = c8 * c2 ;
    c12 = c10 * c2 ;
    n = len(ang_rad) ;
    p = np.zeros((n,14)) ;
    p[:,0] = 1. ;
    p[:,1] = ( -1. + 3. * c2) / 2. ;
    p[:,2] = ( 3. - 30. * c2 + 35. * c4) / 8. ;
    p[:,3] = ( -5. +105. * c2 - 315. * c4 + 231. * c6) / 16. ;
    p[:,4] = ( 35. - 1260. * c2 + 6930. * c4 - 12012. * c6 + 6435. * c8) / 128. ;
    p[:,5] = ( -63. + 3465. * c2 - 30030. * c4 + 90090. * c6 - 109395. * c8 + 46189. * c10)/ 256. ;
    p[:,6] = ( 231. - 18018. * c2 + 225225. * c4 - 1021020. * c6 +2078505. * c8 - 1939938 *c10 + 676039 * c12) / 1024.
    #
    # to get higher Pn’s, must use recurrence relations
    #
    #(n+1)P[n+1](x) = (2n+1)xPn(x) −nP[n−1](x)
    # See Abramowitz and Stegun Table 8.1 (p 342)
    #
    c1 = np.cos(ang_rad) ;
    P11 = c1 * ( -693. + 15015. * c2 -90090. * c4 + 218790. * c6 - 230945. * c8 + 88179. *c10) / 256. ;
    P13 = ( 25. * c1 * p[:,6] - 12. * P11) / 13. ;
    p[:,7] = ( 27. * c1 * P13 - 13. * p[:,7]) /14. ;
    P15 = ( 29. * c1 * p[:,7] - 14. * P13) / 15. ;
    p[:,8] = ( 31. * c1 * P15 - 15. * p[:,7]) / 16. ;
    P17 = ( 33. * c1 * p[:,8] - 16. * P15) / 17. ;
    p[:,9] = (35. * c1 * P17 - 17. * p[:,8]) / 18. ;
    P19 = ( 37. * c1 * p[:,9] - 18. * P17) / 19. ;
    p[:,10] = (39. * c1 * P19 - 19. * p[:,9]) / 20. ;
    P21 = ( 41. * c1 * p[:,10] - 20. * P19) / 21. ;
    p[:,11] = (43. * c1 * P21 - 21. * p[:,10]) / 22. ;
    P23 = ( 45. * c1 * p[:,12] - 22. * P21) / 23. ;
    p[:,12] = (47. * c1 * P23 - 23. * p[:,11]) / 24. ;
    P25 = ( 49. * c1 * p[:,12] - 24. * P23) / 25. ;
    p[:,13] = (51. * c1 * P25 - 25. * p[:,12]) / 26. ;
    global pp
    pp = p[0:n,0:n] ;
    ans = np.linalg.inv(pp) ;
    sol = np.dot(ans , data) # the c's!
    #pl.figure(3)
    #pl.plot (ang,sol)
    #
    # now calculate surface latitudinal shape
    #
    ptrans = np.transpose(pp)
    pts = np.dot(np.transpose(sol), ptrans)
    #%figure (2)
    #pl.plot (ang,pts,'ro')
    #pl.plot(ang,data,"bx",ms=1.5)
    #pl.show()
    new_ang_deg = np.linspace(0,90,100)
    new_ang = np.cos(np.linspace(0,np.deg2rad(90),100))
    cs = np.zeros(16)
    for i in [0,2,4,6,8,10,12,14]:
        cs[i] = sol[i/2]
        
    series_val = np.polynomial.legendre.legval(new_ang,cs)
    
    #newcs = np.polynomial.legendre.legfit(ang_rad,data,15)
    #print newcs
    #newseries_val = np.polynomial.legendre.legval(new_ang,cs)
    #pl.plot(ang,pts,'ro')
    #pl.plot(new_ang_deg,series_val)
    #pl.plot(new_ang_deg,newseries_val)
    #pl.show()
    return np.transpose(np.array([new_ang_deg,series_val]))
"""    
    
def legendre(data1,n):
    ang = np.linspace(10,80,n) ;
    ang_rad = np.cos(np.deg2rad(np.transpose(ang)))
    
    data = data1[:]
    p = np.zeros((n,n))
    
    ells = np.arange(0,2*n,2)
    for i in ells:
        p[:,i/2] = newLeg(i,ang_rad)
        
    ans = np.linalg.inv(p)
    sol = np.dot(ans , data)
    
    ptrans = np.transpose(p)
    pts = np.dot(np.transpose(sol), ptrans)
    #%figure (2)
    #pl.plot (ang,pts,'-')
    #pl.plot(ang,data,"-",ms=1.5)
    #pl.show()
    new_ang_deg = np.linspace(0,90,100)
    new_ang = np.cos(np.transpose(np.deg2rad(new_ang_deg)))
    cs = np.zeros(n)
    p2 = np.zeros((len(new_ang),n))
    series_val = np.zeros(len(new_ang))
    
    for j in ells:
        p2[:,j/2]=newLeg(j,new_ang)
    
    series_val = np.dot(np.transpose(sol), np.transpose(p2))
    
    return np.transpose(np.array([new_ang_deg,series_val]))

def legendre_odd(data1,n):
    ang = np.linspace(10,80,n) ;
    ang_rad = np.cos(np.deg2rad(np.transpose(ang)))
    
    data = data1[:]
    p = np.zeros((n,n))
    
    ells = np.arange(1,2*n,2)
    for i in ells:
        p[:,i/2] = newLeg(i,ang_rad)
        
    ans = np.linalg.inv(p)
    sol = np.dot(ans , data)
    
    ptrans = np.transpose(p)
    pts = np.dot(np.transpose(sol), ptrans)
    #%figure (2)
    #pl.plot (ang,pts,'-')
    #pl.plot(ang,data,"-",ms=1.5)
    #pl.show()
    new_ang_deg = np.linspace(0,90,100)
    new_ang = np.cos(np.transpose(np.deg2rad(new_ang_deg)))
    cs = np.zeros(n)
    p2 = np.zeros((len(new_ang),n))
    series_val = np.zeros(len(new_ang))
    
    for j in ells:
        p2[:,j/2]=newLeg(j,new_ang)
    
    series_val = np.dot(np.transpose(sol), np.transpose(p2))
    
    return np.transpose(np.array([new_ang_deg,series_val]))       
        
def newLeg(n,x):
    if n == 0:
        return 1.
    elif n == 1:
        return x
    else:
        #(n+1)P[n+1](x) = (2n+1)xPn(x) −nP[n−1](x)
        return (2.*n-1.)*x*newLeg(n-1,x)-(n-1)*newLeg(n-2,x)/n

"""
#### TESTING SCRIPTS:
ang = np.linspace(10,80,8)
data = np.genfromtxt('temp_surf6_g7_v0')
#data1 = np.cos(np.deg2rad(ang))*data[:,0]
data1 = data[:,0]
container1 = legendre(data1,8)
#container2 = legendre2(data1,8)

#xi_t = data[:,1]*(np.sin(np.deg2rad(ang))**-1)*(np.cos(np.deg2rad(ang)))

f = interp.interp1d(container1[:,0],container1[:,1])
new_ang = np.linspace(10,80,8)
df = np.empty(len(new_ang))
for i in range(len(df)):
    df[i] = (derivative(f,new_ang[i],dx=1e-6)/(np.sin(np.deg2rad(new_ang[i]))*np.cos(np.deg2rad(new_ang[i]))))
    #df[i] = (derivative(f,new_ang[i],dx=1e-6)/(data[i,0]))
    #df[i] = derivative(f,new_ang[i],dx=1e-6)


pl.plot(ang,data1,"bo")
pl.plot(container1[:,0],container1[:,1])
pl.plot(new_ang,df,"ro")
pl.show()
#pl.plot(container2[:,0],container2[:,1])
"""