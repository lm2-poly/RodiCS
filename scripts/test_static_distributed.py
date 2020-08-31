import numpy as np
import matplotlib.pyplot as pp

from mpl_toolkits.mplot3d import Axes3D

from fenics import Mesh

from static         import run_static
from visualisation  import prepare_frame, plotIt
from postprocessing import *
# =============================================================================
# Mesh
#name = 'straight_10'
name = 'straight_100'
#name = 'straight_300'

mesh_file = '../xml_files/%s.xml' % name

# Define speed direction
u0_th = (0.50)*np.pi
u0_ph = (0.00)*np.pi

# u0 is unitary by definition.
u0 = np.array([np.sin(u0_th)*np.cos(u0_ph),
               np.sin(u0_th)*np.sin(u0_ph),
               np.cos(u0_th)])

which_force = 'distributed'
force_mag   = 10

meshy   = Mesh(mesh_file)
results = run_static(meshy, u0, which_force, force_mag, relaxation=0.5)

def view_3d():
    """
    If you want to hide the material frame, replace the arrays of results_T,
    results_N, and results_B (i.e. results[0], results[1], and results[2]) with
    an array full of None.
    """
    n_vertices = len(results[0])
    
    hidden_results = np.copy(results)

    hidden_results[0] = np.array([None]*n_vertices*3).reshape(n_vertices,3)
    hidden_results[1] = np.array([None]*n_vertices*3).reshape(n_vertices,3)
    hidden_results[2] = np.array([None]*n_vertices*3).reshape(n_vertices,3)
    
    ax = Axes3D(pp.figure(figsize=(5,3.5)))
    prepare_frame(ax, 50, 30)
    plotIt(ax, meshy, hidden_results, static=True)

view_3d()
pp.show()
# =============================================================================
# Quick postprocessing
#total_L = total_length(meshy, results)

#max_deflection = extract_deflection(meshy, u0, results)
#reconf_number  = extract_reconf_number(meshy, u0, results)
