import numpy as np
import matplotlib.pyplot as plt

import physics as p
import solver as s
import Visualization as v
from constants import *
from initial import *

y_0 = initial_condition()

t, x, y, vx, vy = s.IVP(p.core, y_0, tmax, frames, delay)

ani = v.Draw_ani(t, x, y, vx, vy)

if SAVE_VIDEO:
    ani.save("sail.mp4", dpi=300, savefig_kwargs={'facecolor':COLOR})
plt.show()

