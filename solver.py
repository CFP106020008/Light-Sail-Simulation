from scipy.integrate import solve_ivp
import numpy as np

def IVP(function, y_0, tmax, frames, delay=0):
    sol = solve_ivp(fun=function,
                    t_span=[0, tmax],
                    y0=y_0,
                    t_eval=np.linspace(0,tmax,frames),
                    method='LSODA',
                    args=[delay])

    t = sol.t
    Data = sol.y
    x = Data[0,:]
    y = Data[1,:]
    vx = Data[2,:]
    vy = Data[3,:]
    return [t, x, y, vx, vy]

def BVP(IVP):
    delay = 1
    return t0
