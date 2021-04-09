from constants import *
import numpy as np
# Adjustable parameters

# Simulation properties
tmax = 1e0*yr2s # Simulation time

# Sail properties
aE = 1e-3 # Maxium acceleration of the sail at 1AU.
delay = 0.*yr2s # How long does the sail wait 
                 # on the earth orbit before starting manuver
perihelion = 1e1*Rsun

# Visualization properties
Box_size = 3e11 # Size of the plot
frames = int(1e3) # Output frames
Tracing = False # Viewing the sail with tracing mode.
SAVE_VIDEO = True # Whether you want to save the video
VIDEO_LEN = 10 # s
VIDEO_FPS = 40 # frames/s

COLOR = '#303030'
LineColor = 'silver'

# Set initial condition
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

n_process = 10
