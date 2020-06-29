"""
This scripts only loads the mesh and plot it (or export it if you like).
It can be used to quickly check is the mesh gives the right geometry.
"""

import matplotlib.pyplot as pp

from mpl_toolkits.mplot3d import Axes3D

from fenics import Mesh

from visualisation import prepare_frame
# =============================================================================
# Define mesh
name = 'straight_10'

mesh_file = '../xml_files/%s.xml' % name

mesh = Mesh(mesh_file) 

# It is possible to save the mesh in a XDMF file, or even a VTK file
#file_results = XDMFFile('output/%s.xdmf' % mesh_name)
#file_results.write(mesh)

coords = mesh.coordinates()

coords_x = coords[:,0]
coords_y = coords[:,1]
coords_z = coords[:,2]

cells = mesh.cells()
#==============================================================================
def plotIt():
    ax = Axes3D(pp.figure(figsize=(5,3.5)))
    prepare_frame(ax, lambd=15, phi=50)
    
    # Connectivity
    for c in cells:
        id1, id2 = c[0], c[1]
        ax.plot([coords_x[id1], coords_x[id2]],
                [coords_y[id1], coords_y[id2]],
                [coords_z[id1], coords_z[id2]],
                marker='o',
                linestyle='-',
                color='gray',
                alpha=0.5)

    width = 1.02
    xc = -width/4
    yc = -0.01
    zc = -0.5
    ax.set_xlim([xc, xc+width])
    ax.set_ylim([yc, yc+width])
    ax.set_zlim([zc, zc+width])

plotIt()
