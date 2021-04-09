import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.animation import FuncAnimation
import matplotlib as mpl
from matplotlib.gridspec import GridSpec
from constants import *
from initial import *
import physics as p
from physics import Etot, predicted_orbit
import multiprocessing as mp
from tqdm import tqdm

def Draw_ani(t, x, y, vx, vy):

    # Visualization Setup
    fig = plt.figure(figsize = (8, 4.5), facecolor=COLOR)
    gs = GridSpec(2, 4, figure=fig)

    # Picture
    ax = fig.add_subplot(gs[:, :2])
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
    predx, predy = predicted_orbit(x[0], y[0], vx[0], vy[0])
    pred_traj, = ax.plot(predx, predy, color='silver', linestyle=':', linewidth=1)
    dot, = ax.plot([], [], color='silver', marker='o', markersize=1, markeredgecolor='w', linestyle='')
    #Vel = ax.text(0.05, 0.9, 'Velocity: {:.2e} m/s'.format(np.sqrt(vx[0]**2 + vy[0]**2)), horizontalalignment='left',
    #              verticalalignment='top', transform=ax.transAxes, color='w')
    #E_tot = ax.text(0.05, 0.85, 'Specific Total Energy: {:.2e} J/kg'.format(p.Etot(x[0], y[0], vx[0], vy[0])), horizontalalignment='left', 
    #                verticalalignment='top', transform=ax.transAxes, color='w')
    #Time = ax.text(0.05, 0.95, 'Time: {:.2f} yr'.format(t[0]/86400/365), horizontalalignment='left', 
    #                verticalalignment='top', transform=ax.transAxes, color='w')

    ax.set_xlim([-Box_size,Box_size])
    ax.set_ylim([-Box_size,Box_size])
    #%%
    # Velocity Plot
    ax1 = fig.add_subplot(gs[0, 2:])
    ax1.set_facecolor(COLOR)
    velline, = ax1.plot(t[0]/yr2s, np.sqrt(vx[0]**2+vy[0]**2), color='silver')
    ax1.spines['bottom'].set_color(LineColor)
    ax1.spines['top'].set_color(LineColor) 
    ax1.spines['right'].set_color(LineColor)
    ax1.spines['left'].set_color(LineColor)
    ax1.set_xlim([0,tmax/yr2s])
    ax1.set_ylim([0,np.max(np.sqrt(vx**2+vy**2))*1.2])
    ax1.tick_params(labelcolor=LineColor, labelsize='medium', width=3, colors=LineColor)
    ax1.ticklabel_format(axis='y', style='sci', useMathText=True, scilimits=(4,5))
    ax1.set_xlabel('Time (yr)')
    ax1.set_ylabel('Velocity (m/s)')
    ax1.xaxis.label.set_color(LineColor)
    ax1.yaxis.label.set_color(LineColor)

    #%%
    # Energy Plot
    ax2 = fig.add_subplot(gs[1, 2:])
    ax2.set_facecolor(COLOR)
    Etotline, = ax2.plot(t[0]/yr2s, Etot(x[0], y[0], vx[0], vy[0]), color='silver')
    ax2.spines['bottom'].set_color(LineColor)
    ax2.spines['top'].set_color(LineColor) 
    ax2.spines['right'].set_color(LineColor)
    ax2.spines['left'].set_color(LineColor)
    ax2.set_xlim([0, tmax/yr2s])
    ax2.set_ylim([np.min(Etot(x, y, vx, vy))*1.2, np.max(Etot(x, y, vx, vy))*1.2])
    ax2.tick_params(labelcolor=LineColor, labelsize='medium', width=3, colors=LineColor)
    ax2.ticklabel_format(style='sci', useMathText=True)
    ax2.set_xlabel('Time (yr)')
    ax2.set_ylabel('Specific total energy (J/kg)')
    ax2.xaxis.label.set_color(LineColor)
    ax2.yaxis.label.set_color(LineColor)

    plt.tight_layout()
    #%%
    ms2AUyr = 86400*365/1.5e11
    def update(i):
        dot.set_data(x[i], y[i])
        line.set_data(x[:i], y[:i])
        velline.set_data(t[:i]/yr2s, np.sqrt(vx[:i]**2+vy[:i]**2))
        Etotline.set_data(t[:i]/yr2s, Etot(x[:i], y[:i], vx[:i], vy[:i]))
        predx, predy = predicted_orbit(x[i], y[i], vx[i], vy[i])
        pred_traj.set_data(predx, predy)
        r = np.sqrt(x[i]**2 + y[i]**2)
        if Tracing:
            ax.set_xlim([-1.5*r,1.5*r])
            ax.set_ylim([-1.5*r,1.5*r])
        O1 = ax.add_patch(sun)
        O2 = ax.add_patch(mercury)
        O3 = ax.add_patch(venus)
        O4 = ax.add_patch(earth)
        O5 = ax.add_patch(mars)
        #Vel.set_text('Velocity: {:.2e} m/s'.format(np.sqrt(vx[i]**2 + vy[i]**2)))
        #Vel.set_text('Velocity: {:.2e} AU/yr'.format(np.sqrt(vx[i]**2 + vy[i]**2)*ms2AUyr))
        #E_tot.set_text('Total Energy: {:.2e} J/kg'.format(Etot(x[i], y[i], vx[i], vy[i])))
        #Time.set_text('Time: {:.2f} yr'.format(t[i]/86400/365))
        return [dot, line, velline, Etotline, pred_traj, O1, O2, O3, O4, O5]


    ani = FuncAnimation(fig=fig, 
                        func=update,
                        frames=frames, 
                        interval=10000/frames, 
                        blit=True, 
                        repeat=False)

    if SAVE_VIDEO:
        Total_frames = VIDEO_FPS*VIDEO_LEN
        Frames = np.linspace(0, frames, Total_frames).astype(int)
        def save(start, end):
            for i in tqdm(range(start, end)):
                update(Frames[i])
                plt.savefig("./Temp/Frame_{:04d}.png".format(i), dpi=300, facecolor=COLOR)
        Pro_List = []
        FPT = int(len(Frames)/(n_process-1))
        resid = len(Frames)%FPT
        for i in range(n_process-1): # Append the processes into the list
            start = int(i*FPT)
            end = int((i+1)*FPT)
            print(start, end)
            Pro_List.append(mp.Process(target=save, args=(start, end)))
        Pro_List.append(mp.Process(target=save, args=((n_process-1)*FPT, len(Frames))))
        for p in Pro_List: # Start running
            p.start()
        for p in Pro_List: # Wait for all the processes to finish before moving on
            p.join()
    return ani
