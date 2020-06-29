"""
Some useful functions for postprocessing.
"""

import numpy as np

from scipy.integrate import trapz
#==============================================================================
def total_length(meshy, results):
    results_W = results[3]
    
    cells = meshy.cells()

    L = 0
    for c in cells:
        dl = results_W[c[1]] - results_W[c[0]]
        L += np.linalg.norm(dl)

    print('Total length = %s' % L)
    
    return L
        
def extract_deflection(meshy, u0, results):
    coords = meshy.coordinates()

    results_W = results[3]

    displacement = results_W - coords
    
    max_deflection = np.max(u0*displacement)

    print('Numerical maximum deflection = %e' % max_deflection)

    return max_deflection

def extract_reconf_number(meshy, u0, results):
    coords = meshy.coordinates()
    
    ds = np.linalg.norm(coords[1] - coords[0])
    
    results_B = results[2]
    
    ub = np.array([np.dot(b,u0) for b in results_B])
    
    dR = np.abs(ub)*ub*ub
    
    R = trapz(dR, dx=ds)
    
    print('trapz R = %f' % R)
    
    return R

def extract_reconf_number_dynamic(meshy, u0, results_B, results_Speed, step=1):
    coords = meshy.coordinates()
    ds = np.linalg.norm(coords[1] - coords[0])
    
    Nt = len(results_Speed) - 1
    
    list_R = np.empty(0)
    for n in range(0, Nt+1, step):
        ub = np.array([np.dot(b,u0) for b in results_B[n].transpose()])
        
        speed_ub = np.array([sp*UB for sp, UB in zip(results_Speed[n], ub)])
        
        dR = np.abs(speed_ub)*speed_ub*speed_ub
    
        R = trapz(dR, dx=ds)
        
        list_R = np.append(list_R, R)
    
    return list_R
