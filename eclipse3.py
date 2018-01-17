#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 13:53:02 2017

@author: bonar
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 31 17:56:37 2017

@author: bonar
"""

import matplotlib.pyplot as plt
import numpy as np
import math


from astropy import units as u
from poliastro.twobody import Orbit
from poliastro.plotting import plot
from poliastro.plotting import OrbitPlotter
from poliastro.bodies import Earth, Mars, Sun, Moon
import poliastro
plt.style.use("seaborn")



"""
get the satellite position at time after perigee passage
"""


def sun_angle(step, T, r_p, i, phi_0, mars_0_angle):
    T_mars = (24*60+40)*60           #seconds
    G = 6.67408e-11                  #m^3 kg^-1 s^-2
    M_mars = 6.39e23                 #kg
    R_mars = 3390                    #km
    station_lat = 55                 #deg
    mu = 4.282837e13                 #m^3 s^-2 of mars
    k = mu*1e-9                      #km^3 s^-2
                   #km


    r_asc = 0                                       #deg
    arg_p = 90                                      #deg
    true_anom = 0                                  #deg at time t = 0
    a = ((((T**2)/(4*math.pi**2))  * mu  )**(1/3)) * 1e-3  #km
    e = 1-(r_p)/a
    
    orbit_classical = [mu, a, e, i, r_asc, arg_p, true_anom]
    #print("a = " +str(a)+ "\ne = "+str(e)+"\ni = " +str(i)+ "\nr_asc = " +str(r_asc)+"\n arg_p = " +str(arg_p)+"\ntrue_nom = " +str(true_anom))
    
    T_calc = (2*math.pi*math.sqrt(((a*10**3)**3)/mu))
    v_p = math.sqrt(k*((2/(a*(1-e)))-(1/a)))
    r0 = [(r_p * math.cos(math.radians(i))), 0, (r_p * math.sin(math.radians(i)))] * u.km
    v0 = [0, v_p, 0]  * u.km
    
    
    tof = step
    try:
        prop = poliastro.twobody.propagation.kepler(k, r0, v0, tof, rtol=1e-15, numiter=35)
    except:
        prop = None
        while prop is None:
            try:
                tof -= T
                prop = poliastro.twobody.propagation.kepler(k, r0, v0, tof, rtol=1e-15, numiter=35)
                print("trying")
            except:
                 pass
    r_o = prop[0]           #satellite position w.r.t mars at t
    v_o = prop[1]


    r_mars = 2.28e8   
    W_mars_sun = 2*math.pi/686.98 #rad per day
    #position of mars in sun frame

    x = r_mars * math.cos(W_mars_sun*step/60/60 + mars_0_angle)
    y = r_mars * math.sin(W_mars_sun*step/60/60+ mars_0_angle)
    z = 0
    r_mars_t = [x, y, z]
    # print(math.degrees(math.atan(y/x)), r_o, r_mars_t)

    
    phi = math.acos((np.dot(r_o, r_mars_t))/(np.linalg.norm(r_o)*np.linalg.norm(r_mars_t)))
    phi_e = math.asin(3390/np.linalg.norm(r_o))
    if phi < phi_e:
        print("eclipse!")
    
    return(math.degrees(phi), math.degrees(phi_e))
    
"""find angle between landing site vector and vector from landing site to orbiter"""

    
    
T_mars = (24*60+40)*60             
T =  2*T_mars       #period of orbit seconds
r_p = 20331               #km radius of pericentr

i = -63.435     
mars_0_angle = math.radians(0)
sim_time = 1000 #simulation running time in hours
data_step = 0.5   #time between measurements in hours
phis = []
phies = []
ts = []

 
for j in range(0, int(sim_time/data_step)):
    print(str(j)+"/"+str(sim_time/data_step))
    time = j*data_step
    ts.append(time)
    phis.append(sun_angle(time*60*60, T, r_p, i, 0, mars_0_angle)[0])
    phies.append(sun_angle(time*60*60, T, r_p, i, 0, mars_0_angle)[1])
   
plt.xlabel("Time (hours)")
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

plt.axhline(y=90, color='k', linestyle='--')
plt.ylabel("Angles between position vectors (degrees)")
plt.title(r"$\theta_s $ (blue) and $\theta_e$ (red) for Orbit Perpendicular to Sun", fontsize = 20)  
plt.fill_between(ts, 0, phies, facecolor='red', alpha=0.9)  
plt.plot(ts, phis, alpha=0.6)


plt.ylim(0, 90)



tes = []
eclipse = False
starts = []
ends = []
for i in range(len(phis)):
    if eclipse is False:     
        if phis[i] < 90 and phis[i] < phies[i]:
            start = ts[i]
            starts.append(start)
            print("start eclipse")
            eclipse = True
        else:
            pass
    else:
        if phis[i] < 90 and phis[i] < phies[i]:
            pass
        else:
            end = ts[i]
            ends.append(end)
            print("end eclipse")
            eclipse_time = end-start
            tes.append(eclipse_time)
            eclipse = False
max_index = tes.index(max(tes))
maxwindow = max(tes)
print("maximum eclipse time is "+str(max(tes))+" hours \nThis occurs beginning at time "+str(starts[max_index]))    
plt.annotate("Maximum Eclipse Time = "+str(maxwindow)+" hours\n on day "+str(starts[max_index]/24), xy=(0.05, 0.85), xycoords = 'axes fraction', fontsize=16)
plt.annotate("There are "+str(len(tes))+" eclipses",  xy=(0.05, 0.75), xycoords = 'axes fraction', fontsize=16)
plt.show()