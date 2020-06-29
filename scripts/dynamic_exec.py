import numpy as np

from fenics import Mesh

from dynamic       import initialise_results, run_dynamic
from visualisation import *
from axes_world    import one_by_one
# =============================================================================
# Mesh
name = 'straight_10'
#name = 'straight_30'

mesh_file = '../xml_files/%s.xml' % name

# Define speed direction
u0_th = (0.50)*np.pi
u0_ph = (0.00)*np.pi

# u0 is unitary by definition.
u0 = np.array([np.sin(u0_th)*np.cos(u0_ph),
               np.sin(u0_th)*np.sin(u0_ph),
               np.cos(u0_th)])

omeg  = 2*np.pi
def speed_func(t_n):
    """
    Returns the normalised velocity U(t)/U0.
    """
    return min(1,t_n)

Cd    = 1.2
St    = 0.16
Gamma = 10
Cy    = 100
Ur    = St * np.sqrt(np.pi * Cy / 4.)

dt = 5e-4
Tf = 1
Nt = int(Tf/dt)

rlx = 0.8

which_force = 'distributed'
force_mag   = 10
   
meshy   = Mesh(mesh_file)
results = initialise_results()

run_dynamic(results, meshy, u0, which_force, force_mag, dt, Nt, Ur, Gamma,
            speed_func, rlx)

def plot_at_n(n):
    """
    If you want to hide the material frame, replace the arrays of results_T,
    results_N, and results_B (i.e. results[0][n], results[1][n], and results[2][n])
    with an array full of None.
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
    
n = int(0.5*Nt)
plot_at_n(n)
# =============================================================================
# Image sequence saving.
# You can animate your sequence with the software ffmpeg or imageJ.
#nb_frames = 100
#save_sequence(Nt//nb_frames, dt, meshy, results, 'test_dynamic')
# =============================================================================
# Analysis
#wx_tip = results[3][:,-1,0]

#ax = one_by_one()
#plot_vs_time(ax, wx_tip, dt, r'$w_{0}$', color='black')

