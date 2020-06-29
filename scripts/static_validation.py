"""
Code validation for the benchmark of a rod under uniformly distributed load.
"""

import numpy as np

from fenics import Mesh

from static         import run_static
from postprocessing import extract_deflection

from axes_world import one_by_one, fontsize
# =============================================================================
# Theoretical deflection by Rohde (1952)
def g0(alpha):
    discriminant = np.cos(alpha)**2/36. + alpha*np.sin(alpha)*np.cos(alpha)/45.
    return 90*(np.sqrt(discriminant) - np.cos(alpha)/6.)/(np.sin(alpha)*np.cos(alpha))

def a3(alpha):
    G0 = g0(alpha)
    return -G0*np.cos(alpha)/6.

def a6(alpha):
    G0 = g0(alpha)
    return -G0**2*np.cos(alpha)*np.sin(alpha)/180.

def T1(s, alpha):
    A3 = a3(alpha)
    A6 = a6(alpha)
    return s - (A3**2/14.)*s**7 - (A3*A6/10.)*s**10

def T2(s, alpha):
    A3 = a3(alpha)
    A6 = a6(alpha)
    return (A3/4.)*s**4 + (A6/7.)*s**7

def th_y_dist(s, alpha):
    return T1(s, alpha)*np.sin(alpha) + T2(s, alpha)*np.cos(alpha)

# Script for validation
u0 = np.array([1., 0., 0.])

def validation_distributed(mesh_file):
    alphas = np.linspace(0.015, 1.25, 60)
    g0s    = g0(alphas)
    
    delta_ROHDE  = th_y_dist(1, alphas)
    delta_FEniCS = np.empty(0)
    
    relaxation = 0.8
    for G0 in g0s:
        print('g0 = %f' % G0)
        
        if G0 > 10:
            relaxation = 0.5
        
        meshy   = Mesh(mesh_file)
        results = run_static(meshy, u0, 'distributed', G0, relaxation=relaxation)
        
        delta = extract_deflection(meshy, u0, results)
        delta_FEniCS = np.append(delta_FEniCS, delta)
        
        print('')
        
    return g0s, delta_ROHDE, delta_FEniCS

# =============================================================================
# Graphics function
def plot_validation_distributed(ax, g0s, delta_ROHDE, delta_FEniCS):
    ax.plot(g0s, delta_FEniCS, color='goldenrod', linewidth=2,
            linestyle='-', label='FEniCS')
    
    ax.plot(g0s, delta_ROHDE, color='black', linewidth=1,
            linestyle='--', label='Rohde (1952)')

    ax.legend(loc='lower right',
              fontsize=fontsize,
              frameon=False,
              ncol=1,
              labelspacing=0.3,
              handlelength=1.2)
    
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8])
    ax.set_xticks(np.arange(0, 16, 5))
    
    ax.set_xlabel(r'$g_{0}$', fontsize=fontsize)
    ax.set_ylabel(r'$\delta$', fontsize=fontsize)
    
# =============================================================================
name = 'straight_300'
mesh_file = '../xml_files/%s.xml' % name

g0s, delta_ROHDE, delta_FEniCS = validation_distributed(mesh_file)

ax = one_by_one()
plot_validation_distributed(ax, g0s, delta_ROHDE, delta_FEniCS)