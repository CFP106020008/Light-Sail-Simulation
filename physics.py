import numpy as np
from constants import *
from initial import *

def Etot(x, y, vx, vy):
    r = np.sqrt(x**2 + y**2)
    K = 0.5*(vx**2 + vy**2)
    return K - G*Msun/r

def predict_peri(x, y, vx, vy):
    E = Etot(x, y, vx, vy)
    vec_r = np.array([x, y]) 
    vec_v = np.array([vx, vy]) 
    a = -G*Msun/2/E
    vec_e = ((np.linalg.norm(vec_v)**2 - G*Msun/np.linalg.norm(vec_r))*vec_r - np.dot(vec_r, vec_v)*vec_v)/(G*Msun)
    e = np.linalg.norm(vec_e)
    p = a*(1-e*e)
    return p

def Decide_Pointing(t, x, y, vx, vy, delay=0):
    r = np.array([x, y])
    v = np.array([vx, vy])
    rhat = r/np.linalg.norm(r)
    vhat = v/np.linalg.norm(v)
    Acc = False
    #limit = predict_peri(x, y, vx, vy)
    limit = np.linalg.norm(r)
    if t > delay:
        if limit > perihelion and Acc == False:
            phat = (rhat - vhat)/np.sqrt(2)
        else:
            Acc = True
            phat = (rhat + vhat)/np.sqrt(2)
            #phat = np.array([-rhat[1], rhat[0]])
    else:
        phat = np.array([-rhat[1], rhat[0]])
    return phat
#%%
def core(t, y, delay=0):
    r_vec = y[:2]
    r = np.linalg.norm(r_vec)
    v_vec = y[2:]
    phat = Decide_Pointing(t, y[0], y[1], y[2], y[3], delay)
    dxdt = v_vec[0]
    dydt = v_vec[1]
    a_rp = aE*rE**2/r**2*np.dot(r_vec, phat)/np.linalg.norm(r_vec)*phat
    a_g = -G*Msun/r**3*r_vec
    a = a_rp + a_g
    dvxdt = a[0]
    dvydt = a[1]
    return np.array([dxdt, dydt, dvxdt, dvydt])

def predicted_orbit(x, y, vx, vy):
    E = Etot(x, y, vx, vy)
    if E < 0:
        vec_r = np.array([x, y]) 
        vec_v = np.array([vx, vy]) 
        a = -G*Msun/2/E
        vec_e = ((np.linalg.norm(vec_v)**2 - G*Msun/np.linalg.norm(vec_r))*vec_r - np.dot(vec_r, vec_v)*vec_v)/(G*Msun)
        e = np.linalg.norm(vec_e)
        h = np.linalg.norm(np.cross(np.array([vec_r[0], vec_r[1], 0]), np.array([vec_v[0], vec_v[1], 0])))
        theta = np.linspace(0, 2*np.pi, 1000)
        r_predict = h**2/G/Msun/(1+np.linalg.norm(vec_e)*np.cos(theta))
        P = np.array([r_predict*np.cos(theta), -r_predict*np.sin(theta)])
        if e == 0:
            R = np.array([[vec_e[0], -vec_e[1]], [vec_e[1], vec_e[0]]])
        else:
            R = np.array([[vec_e[0], -vec_e[1]], [vec_e[1], vec_e[0]]])/e
        OUT = np.matmul(R, P)
    else: 
        OUT = [[], []]
    return OUT
