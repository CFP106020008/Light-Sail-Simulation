import numpy as np
import matplotlib.pyplot as plt
import os
import time

import physics as p
import solver as s
import Visualization as v
from constants import *
from initial import *

t0 = time.time()

y_0 = initial_condition()
print("Initial condition set up.")

t, x, y, vx, vy = s.IVP(p.core, y_0, tmax, frames, delay)
print("Orbit integration complete.")

ani = v.Draw_ani(t, x, y, vx, vy)
print("Initialize the plot")

if SAVE_VIDEO:
    print("Start saving movie.")
    #ani.save("sail.mp4", dpi=300, savefig_kwargs={'facecolor':COLOR})
    print("Movie saved.")
    os.system('ffmpeg -y -r 40 -pattern_type glob -i "./Temp/*.png" -vcodec libx264 -s 2048x1152 -pix_fmt yuv420p sail.mp4')
    os.system("rm ./Temp/*.png")
else:
    plt.show()

t1 = time.time()
print("The program took {} s to run".format(t1-t0))
