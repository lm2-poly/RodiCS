"""
Code verification for the case of a rod under uniformly distributed load.
"""

import numpy as np

from fenics import Mesh

from static import run_static

from axes_world import one_by_two
# =============================================================================
u0 = np.array([1., 0., 0.])

# Verification loop
def verification_distributed_L2(list_Nelem, g0):
    profiles = []

    for n_elems in list_Nelem:
        print('number of elements = %d' % n_elems)
        
        mesh_file = '../xml_files/straight_%s.xml' % n_elems
        
        meshy   = Mesh(mesh_file)
        results = run_static(meshy, u0, 'distributed', g0, relaxation=0.5)
        
        # x-displacement for all nodes
        displ_x = results[3,:,0]
        profiles.append(displ_x)
        
    l2s = np.empty(0)
    for i, prof in enumerate(profiles[:-1]):
        diff = profiles[i] - profiles[i+1][::2]
        l2 = np.sqrt(np.sum(diff**2)/len(profiles[i]))
        
        l2s = np.append(l2s, l2)
      
    return l2s

# =============================================================================
# Graphics functions
def plot_discretisation_error_L2(ax, list_Nelem, errorsL2, color, marker):
    ax.plot(list_Nelem[1:], errorsL2, linestyle='-',
            linewidth=1, color=color, marker=marker, markeredgecolor='black',
            markeredgewidth=0.5, alpha=0.75)

    ax.set(xscale='log', yscale='log')
    
    ax.set_ylabel(r'$\varepsilon_{N}$', fontsize=12)
    ax.set_xlabel(r'$N$', fontsize=12)
    
def plot_estimated_order_L2(ax, list_Nelem, errorsL2, color, marker):
    obs_p = np.log(errorsL2[:-1]/errorsL2[1:])/np.log(2)
    print('observed orders = %s' % obs_p)
    
    ax.plot(list_Nelem[2:], obs_p, linestyle='',
            linewidth=1, markeredgewidth=0.5, color=color, marker=marker,
            markeredgecolor='black', alpha=0.75)
    
    ax.set_xscale('log')
    
    ax.set_ylim([1.89, 2.11])
    ax.set_yticks([1.9, 1.95, 2, 2.05, 2.1])

    ax.set_xlabel(r'$N$', fontsize=12)
    ax.set_ylabel(r'$\hat{p}_{N}$', fontsize=12)
# =============================================================================
list_Nelem = np.array([50, 100, 200, 400, 800, 1600, 3200, 6400, 12800])

light_load = verification_distributed_L2(list_Nelem, g0=0.1)
heavy_load = verification_distributed_L2(list_Nelem, g0=10.)

ax_a, ax_b = one_by_two()

ax_a.plot(list_Nelem[1:5], 5e-5/list_Nelem[:4]**2, linestyle='--',
          color='black')
ax_b.plot(list_Nelem[2:], 2 + 0*list_Nelem[2:], linestyle='--',
          color='black')

plot_discretisation_error_L2(ax_a, list_Nelem, heavy_load, 'red', 's')
plot_discretisation_error_L2(ax_a, list_Nelem, light_load, 'blue', '^')

plot_estimated_order_L2(ax_b, list_Nelem, heavy_load, 'red', 's')
plot_estimated_order_L2(ax_b, list_Nelem, light_load, 'blue', '^')
