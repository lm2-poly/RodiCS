import numpy as np

from fenics import Mesh

from dynamic       import initialise_results, run_dynamic
from visualisation import *
from axes_world    import one_by_one
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

omeg  = 6.4*np.sqrt(0.028)
def sine_func(t_n):
    """
    Returns the normalised velocity U(t)/U0.
    """
    return np.sin(omeg*t_n)

Cd    = 1.2
St    = 0.16
Gamma = 10

# These are formulae used in Leclercq and de Langre (2018). In my thesis the
# pulsation is normalised by another constant, whence the factor sqrt(0.028).
alpha = 0.65
Cy    = 12.7 * ( (omeg/np.sqrt(0.028)) * alpha)**2 / 2.0
Ur    = St * alpha*omeg*Gamma

dt = 5e-4
Nt = 4*int(2*np.pi/omeg/dt)

rlx = 0.8

meshy   = Mesh(mesh_file)
results = initialise_results()

run_dynamic(results, meshy, u0, 'drag', Cd*Cy, dt, Nt, Ur, Gamma, sine_func, rlx)

# Image sequence saving.
# You can animate your sequence with the software ffmpeg or imageJ.
nb_frames = 100
save_sequence(Nt//nb_frames, dt, meshy, results, 'output/test_oscillatory_dynamic')
# =============================================================================
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
    
#n = int(0.75*Nt)
#plot_at_n(n)
# =============================================================================
# Analysis
#wx_tip = results[3][:,-1,0]
#wx_all, wx_all = results[3][:,:,0], results[3][:,:,1]

#ax = one_by_one()
#plot_vs_time(ax, wx_tip, dt, r'$w_{0}$', color='black')
#plot_profiles(ax, wx, wy)

