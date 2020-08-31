"""
Code verification using the benchmark of a rod under uniformly distributed load.
"""

import numpy as np

from fenics import Mesh

from dynamic    import initialise_results, run_dynamic
from axes_world import one_by_two, fontsize
# =============================================================================
# Verification
u0 = np.array([1., 0., 0.])
Tf = 1

Ur, Gamma = 5, 50 # not necessary to simulate distributed load

def verification_distributed_L2(list_dt, mesh_file, g0):
    
    def speed_function(t_n):
        return 1.

    meshy   = Mesh(mesh_file)

    tips_evolution = []
    for dt in list_dt:
        print('dt = %e' % dt)
        
        Nt  = int(Tf/dt)
        rlx = 0.8
        
        results = initialise_results()

        run_dynamic(results, meshy, u0, 'distributed', g0, dt, Nt, Ur, Gamma, 
                    speed_function, rlx)
                    
        wx_tip = results[3][:,-1,0]
        tips_evolution.append(wx_tip)
        
    l2s = np.empty(0)
    for i, tip in enumerate(tips_evolution[:-1]):
        diff = tip[i] - tip[i+1][::2]
        l2 = np.sqrt(np.sum(diff**2)/len(tips_evolution[i]))
        
        l2s = np.append(l2s, l2)
      
    return l2s

# =============================================================================
# Graphics
def plot_discretisation_error_L2(ax, list_dt, errorsL2, color, marker):
    ax.plot(list_dt[1:], errorsL2, linestyle='-',
            linewidth=1, color=color, marker=marker, markeredgecolor='black',
            markeredgewidth=0.5, alpha=0.75)

    ax.set(xscale='log', yscale='log')
    
    ax.set_ylabel(r'$\varepsilon_{\Delta t}$', fontsize=fontsize)
    ax.set_xlabel(r'$\Delta t$', fontsize=fontsize)
    
def plot_estimated_order_L2(ax, list_dt, errorsL2, color, marker):
    obs_p = np.log(errorsL2[:-1]/errorsL2[1:])/np.log(2)
    print('observed orders = %s' % obs_p)
    
    ax.plot(list_dt[2:], obs_p, linestyle='', linewidth=1,
            color=color, marker=marker, markeredgecolor='black',
            markeredgewidth=0.5, alpha=0.75)
    
    ax.set_xscale('log')
    
    ax.set_xlabel(r'$\Delta t$', fontsize=fontsize)
    ax.set_ylabel(r'$\hat{p}_{\Delta t}$', fontsize=fontsize)

# =============================================================================
list_dt = np.array([2e-3, 1e-3, 5e-4, 2.5e-4, 1.25e-4, 6.25e-5, 3.125e-5])

light_load = verification_distributed_L2(list_dt, g0=0.1)
heavy_load = verification_distributed_L2(list_dt, g0=10.)

ax_a, ax_b = one_by_two()

ax_a.plot(list_dt[1:5], 5e-5/list_dt[:4]**2, linestyle='--',
          color='black')
ax_b.plot(list_dt[2:], 2 + 0*list_dt[2:], linestyle='--',
          color='black')

plot_discretisation_error_L2(ax_a, list_dt, heavy_load, 'red', 's')
plot_discretisation_error_L2(ax_a, list_dt, light_load, 'blue', '^')

plot_estimated_order_L2(ax_b, list_dt, heavy_load, 'red', 's')
plot_estimated_order_L2(ax_b, list_dt, light_load, 'blue', '^')

