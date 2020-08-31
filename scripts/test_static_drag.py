import numpy as np
import matplotlib.pyplot as pp

from mpl_toolkits.mplot3d import Axes3D

from fenics import Mesh

from static         import run_static
from visualisation  import *

from axes_world import one_by_one
# =============================================================================
# Mesh
#name = 'straight_50'
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

# Drag coefficient and Cauchy number
Cd = 1.2
Cy = 50

# Relaxation parameter
rlx = 0.8

# Import mesh
meshy   = Mesh(mesh_file)

# Successive simulations
def run_succession(list_Cy):
    list_results = []
    
    meshy = Mesh(mesh_file)
    for Cy in list_Cy:
        print('Cy = %s' % Cy)
        results = run_static(meshy, u0, 'drag', Cd*Cy, rlx)
        list_results.append(results)
        print('')
    
    return list_results

def draw_succession(ax, list_Cy, list_results, colors):
    for Cy, results, color in zip(list_Cy, list_results, colors):
        pos_x, pos_y = results[3,:,0], results[3,:,1]
        
        ax.plot(pos_y, -pos_x, linewidth=1.5, linestyle='-', color=color,
                label=r'$C_{\mathrm{Y}} = %.1f$' % Cy)
        ax.plot(-pos_y, -pos_x, linewidth=1.5, linestyle='-', color=color)
        
        ax.set_aspect('equal')
        ax.axis('off')
        
        ax.legend(loc='best',
                  fontsize=10,
                  frameon=False,
                  ncol=1,
                  labelspacing=0.3,
                  handlelength=1.5)

# Parameters used to plot Figure 1(a) in the JOSS paper
list_Cy_tilde = np.array([0.001, 2, 4.7, 9.0, 27, 72, 99])
list_Cy       = list_Cy_tilde/Cd
colors  = ['lightgray', 'yellow', 'khaki', 'goldenrod', 'olive', 'orange', 'red']

list_results = run_succession(list_Cy)

ax = one_by_one()
draw_succession(ax, list_Cy, list_results, colors)
pp.show()