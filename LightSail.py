# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 01:09:23 2021

@author: juliu
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import matplotlib.cm as cm
from matplotlib.animation import FuncAnimation
import matplotlib as mpl

# Fundamental Constants
Msun = 2e30      # Solar Mass (k)
G = 6.67e-11     # Gravitational Constant
AU = 1.5e11      # Astronmical Unit (m)
rE = 1.*AU       # Orbital radius of Earth (m)
Rsun = 7e8       # Radius of the sun (m)
yr2s = 86400*365 # Conversion of year to second

#%%
# Adjustable parameters

# Simulation properties
tmax = 2e0*yr2s # Simulation time

# Sail properties
aE = 1e-3 # Maxium acceleration of the sail at 1AU.
delay = 0.5*yr2s # How long does the sail wait 
                 # on the earth orbit before starting manuver

# Visualization properties
Box_size = 3e11 # Size of the plot
frames = int(1e3) # Output frames
Tracing = False # Viewing the sail with tracing mode.
SAVE_VIDEO = False # Whether you want to save the video

#%%
def initial_condition():
    def cicular(k):
        v_sun_circular = np.sqrt(G*Msun/k)
        return [0, -k, v_sun_circular, 0]
    def eliptical():
        a = AU
        R_init = a*2-1*Rsun
        v = np.sqrt(G*Msun*(2/R_init-1/a))
        return [R_init, 0, 0, v]
    return cicular(AU)
y_0 = initial_condition()
#%%
def Etot(x, y, vx, vy):
    r = np.sqrt(x**2 + y**2)
    K = 0.5*(vx**2 + vy**2)
    return K - G*Msun/r

#%%
def Decide_Pointing(t, x, y, vx, vy):
    r = np.array([x, y])
    v = np.array([vx, vy])
    rhat = r/np.linalg.norm(r)
    vhat = v/np.linalg.norm(v)
    Acc = False
    if t > delay:
        if np.linalg.norm(r)>10*Rsun and Acc == False:
            phat = (rhat - vhat)/np.sqrt(2)
        else:
            Acc = True
            phat = (rhat + vhat)/np.sqrt(2)
    else:
        phat = np.array([-rhat[1], rhat[0]])
    return phat
#%%
def function(t, y):
    r_vec = y[:2]
    r = np.linalg.norm(r_vec)
    v_vec = y[2:]
    phat = Decide_Pointing(t, y[0], y[1], y[2], y[3])
    dxdt = v_vec[0]
    dydt = v_vec[1]
    a_rp = aE*rE**2/r**2*np.dot(r_vec, phat)/np.linalg.norm(r_vec)*phat
    a_g = -G*Msun/r**3*r_vec
    a = a_rp + a_g
    dvxdt = a[0]
    dvydt = a[1]
    return np.array([dxdt, dydt, dvxdt, dvydt])
#%%
sol = solve_ivp(fun=function,
                t_span=[0, tmax],
                y0=y_0,
                t_eval=np.linspace(0,tmax,frames),
                method='DOP853')

t = sol.t
Data = sol.y
x = Data[0,:]
y = Data[1,:]
vx = Data[2,:]
vy = Data[3,:]

#%%
# Visualization Setup
COLOR = '#303030'
fig, ax = plt.subplots(facecolor=COLOR)
ax.set_facecolor(COLOR)
ax.xaxis.set_visible(False)
ax.yaxis.set_visible(False)
ax.spines['bottom'].set_color(COLOR)
ax.spines['top'].set_color(COLOR) 
ax.spines['right'].set_color(COLOR)
ax.spines['left'].set_color(COLOR)

# Solar system bodies
sun = plt.Circle((0, 0), Rsun, color='y')
mercury = plt.Circle((0, 0), 0.387*AU, edgecolor='cyan', fill=False)
venus = plt.Circle((0, 0), 0.723*AU, edgecolor='y', fill=False)
earth = plt.Circle((0, 0), 1.*AU, edgecolor='skyblue', fill=False)
mars = plt.Circle((0, 0), 1.524*AU, edgecolor='r', fill=False)

ax.add_patch(sun)
ax.add_patch(mercury)
ax.add_patch(venus)
ax.add_patch(earth)
ax.add_patch(mars)

ax.set_aspect('equal', 'box')

line, = ax.plot(x[0], y[0], color='silver', linestyle='-', linewidth=1)
dot, = ax.plot([], [], color='silver', marker='o', markersize=1, markeredgecolor='w', linestyle='')
Vel = ax.text(0.05, 0.9, 'Velocity: {:.2e} m/s'.format(np.sqrt(vx[0]**2 + vy[0]**2)), horizontalalignment='left',
              verticalalignment='top', transform=ax.transAxes, color='w')
E_tot = ax.text(0.05, 0.85, 'Specific Total Energy: {:.2e} J/kg'.format(Etot(x[0], y[0], vx[0], vy[0])), horizontalalignment='left', 
                verticalalignment='top', transform=ax.transAxes, color='w')
Time = ax.text(0.05, 0.95, 'Time: {:.2f} yr'.format(t[0]/86400/365), horizontalalignment='left', 
                verticalalignment='top', transform=ax.transAxes, color='w')

ax.set_xlim([-Box_size,Box_size])
ax.set_ylim([-Box_size,Box_size])
plt.tight_layout()
#Sail = mpl.patches.Rectangle((x[0], y[0]), width=5, height=50, angle=0, color='b')
#ax.add_patch(Sail)
#%%
ms2AUyr = 86400*365/1.5e11
def update(i):
    dot.set_data(x[i], y[i])
    line.set_data(x[:i], y[:i])
    r = np.sqrt(x[i]**2 + y[i]**2)
    if Tracing:
        ax.set_xlim([-1.5*r,1.5*r])
        ax.set_ylim([-1.5*r,1.5*r])
    #Vel.set_text('Velocity: {:.2e} m/s'.format(np.sqrt(vx[i]**2 + vy[i]**2)))
    Vel.set_text('Velocity: {:.2e} AU/yr'.format(np.sqrt(vx[i]**2 + vy[i]**2)*ms2AUyr))
    E_tot.set_text('Total Energy: {:.2e} J/kg'.format(Etot(x[i], y[i], vx[i], vy[i])))
    Time.set_text('Time: {:.2f} yr'.format(t[i]/86400/365))
    return [dot, line, Vel, E_tot, Time]


ani = FuncAnimation(fig=fig, 
                    func=update,
                    frames=frames, 
                    interval=10000/frames, 
                    blit=True, 
                    repeat=False)

if SAVE_VIDEO:
    ani.save("sail.mp4", dpi=300)
plt.show()