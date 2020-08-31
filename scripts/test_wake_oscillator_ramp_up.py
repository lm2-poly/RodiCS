import numpy as np

from fenics import Mesh

from wake_oscillator import initialise_results, run_wake_oscillator
from visualisation   import *
from axes_world      import one_by_one
# =============================================================================
# Mesh
#name = 'straight_10'
name = 'straight_30'

mesh_file = '../xml_files/%s.xml' % name

# Define speed direction
u0_th = (0.50)*np.pi
u0_ph = (0.00)*np.pi

# u0 is unitary by definition.
u0 = np.array([np.sin(u0_th)*np.cos(u0_ph),
               np.sin(u0_th)*np.sin(u0_ph),
               np.cos(u0_th)])

def speed_func(t_n):
    return min(1.,t_n)

Gamma = 25
St    = 0.16

Cy    = 5
Ur    = St * np.sqrt(np.pi * Cy / 4.)

dt = 5e-4
Tf = 10
Nt = int(Tf/dt)

rlx = 0.8

meshy   = Mesh(mesh_file)
results = initialise_results()

run_wake_oscillator(results, meshy, u0, Cy, dt, Nt, Gamma, Ur, speed_func, rlx)

# Image sequence saving.
# You can animate your sequence with the software ffmpeg or imageJ.
nb_frames = 100
save_sequence(Nt//nb_frames, dt, meshy, results[:-1], 'output/test_ramp_up_wake_oscillator')
# =============================================================================
def plot_at_n(n):
    """
    If you want to know how to hide any vector, see the script static_exec.py
    """
    results_T_n = results[0][n]
    results_N_n = results[1][n]
    results_B_n = results[2][n]
    results_W_n = results[3][n]
    
    results_Fext_n  = results[4][n]
    results_Speed_n = results[5][n]

    results_n = results_T_n, results_N_n, results_B_n, results_W_n, \
                results_Fext_n, results_Speed_n

    ax = Axes3D(pp.figure())
    prepare_frame(ax)
    plotIt(ax, meshy, results_n, to_scale=1, freq=1)

#n = int(0.75*Nt)    
#plot_at_n(n)
# =============================================================================
# Analysis
#wx_tip = results[3][:,-1,0]
#wz_tip = results[3][:,-1,2]

#wx_all, wy_all, wz_all = results[3][:,:,0], results[3][:,:,1], results[3][:,:,2]

# Lift fluctuating variable (at the tip)
#q = results[5][:,-1]

# Speed function U(t)/U0, for checking
#sp = results[6][:,-1]

#ax = one_by_one()
#plot_vs_time(ax, wx_tip, dt, r'$w_{0}$', color='forestgreen')
#plot_vs_time(ax, wz_tip, dt, r'$w_{\perp}$', color='firebrick')

#plot_profiles(ax, wx_all, wy_all) # lateral profile
#plot_profiles(ax, wz_all, wy_all) # frontal profile


