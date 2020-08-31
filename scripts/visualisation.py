import numpy as np
import matplotlib.pyplot as pp

from mpl_toolkits.mplot3d import Axes3D
#==============================================================================
fontsize = 12

def prepare_frame(ax, lambd=15, phi=50):
    ax.xaxis.set_pane_color((1., 1., 1., 0.))
    ax.yaxis.set_pane_color((1., 1., 1., 0.))
    ax.zaxis.set_pane_color((1., 1., 1., 0.))

    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])

    ax.axis('off')

    ax.grid(False)

    ax.view_init(lambd, phi)

    # Reference frame
    ax.quiver(0., 0., 0., .1, 0., 0., color='green', linewidth=2)
    ax.quiver(0., 0., 0., 0., .1, 0., color='green', linewidth=2)
    ax.quiver(0., 0., 0., 0., 0., .1, color='green', linewidth=2)
    
    ax.text3D( 0.07, -0.05, 0, '$x$', (-1,1,0), color='green', fontsize=fontsize)
    ax.text3D(-0.05,  0.05, 0, '$y$', (-1,1,0), color='green', fontsize=fontsize)
    ax.text3D( 0.0, 0.0, 0.13, '$z$', (-1,1,0), color='green', fontsize=fontsize)

def plotIt(ax, meshy, results, to_scale=1, static=False, freq=2):
    if static:
        results_T, results_N, results_B, results_W, results_Fext = results
    else:
        results_T, results_N, results_B, results_W, results_Fext,\
        results_Speed = results
    
    # Unfolding the mesh input
    coordinates = meshy.coordinates()

    coords_x = coordinates[:,0]
    coords_y = coordinates[:,1]
    coords_z = coordinates[:,2]
    
    cells = meshy.cells()
    
    # Unfolding the results input
    x = results_W[:,0]
    y = results_W[:,1]
    z = to_scale*results_W[:,2]
    
    if not(static) and results_Speed.all() != None:
        scale = 0.1
        for posi in np.arange(-0.1,0.15,0.1):
            ax.quiver(coords_x[0], coords_y[0] - 0.1, coords_z[0] + posi,
                      scale*results_Speed[0], 0, 0, color='black')
        
    for c in cells:
        id1, id2 = c[0], c[1]
        ax.plot([coords_x[id1], coords_x[id2]],
                [coords_y[id1], coords_y[id2]],
                [coords_z[id1], coords_z[id2]],
                marker='.',
                linestyle='-',
                color='gray',
                alpha=0.5)

        ax.plot([x[id1], x[id2]],
                [y[id1], y[id2]],
                [z[id1], z[id2]],
                marker='.',
                linestyle='-',
                color='purple')

        # Material frame
        if results_T.all() != None and id1%freq == 0:
            scale = 0.1
            ax.quiver(x[id2], y[id2], z[id2],
                      scale*results_T[id2,0],
                      scale*results_T[id2,1],
                      scale*results_T[id2,2],
                      color='black', alpha=0.25)

        if results_N.all() != None and id1%freq == 0:
            scale = 0.1
            ax.quiver(x[id2], y[id2], z[id2],
                      scale*results_N[id2,0],
                      scale*results_N[id2,1],
                      scale*results_N[id2,2],
                      color='black', alpha=0.25)

        if results_B.all() != None and id1%freq == 0:
            scale = 0.1
            ax.quiver(x[id2], y[id2], z[id2],
                      scale*results_B[id2,0],
                      scale*results_B[id2,1],
                      scale*results_B[id2,2],
                      color='black', alpha=0.25)

        # External force
        if results_Fext.all() != None and id1%freq == 0:
            results_Fext_norms = np.linalg.norm(results_Fext, axis=1)
            max_Fext_norms     = np.max(results_Fext_norms)

            scale = 0.15/max_Fext_norms
                
            ax.quiver(x[id2], y[id2], z[id2],
                      scale*results_Fext[id2,0],
                      scale*results_Fext[id2,1],
                      scale*results_Fext[id2,2],
                      color='red', alpha=0.75)

    # This is a workaround to display a figure with equal axes (since there is
    # no ax.axis('equal') in 3D :/ )
    width = 1.02
#    width = 1.15

#    xc = -width/2
#    yc = -0.
#    zc = -0.5
    xc = -width/4
    yc = -0.01
    zc = -0.5
    ax.set_xlim([xc, xc+width])
    ax.set_ylim([yc, yc+width])
    ax.set_zlim([zc, zc+width])

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    
def save_sequence(freq, dt, meshy, results, name, to_scale=1):
    results_T, results_N, results_B, results_W, results_Fext, results_Speed \
    = results
    
    Nt = len(results_W)
    
    ax = Axes3D(pp.figure())

    print('\nSaving image sequence')
    
    for i, h in enumerate(range(freq,Nt,freq)):
        # A progression bar...
        print('%3.1f' % (100*h/Nt) + '% ' + '%s' % ('.'*(50*h//Nt)))
        
        prepare_frame(ax)
        
        results_n = results_T[h], results_N[h], results_B[h], results_W[h], \
                    results_Fext[h], results_Speed[h]
        
        plotIt(ax, meshy, results_n, to_scale=to_scale)
        pp.savefig('%s_%03d.png' % (name, i))
        
        ax.clear()
        
    print('Done')

def plot_vs_time(ax, time_results, dt, label, color='black'):
    Nt = len(time_results) - 1
    timeline = np.linspace(0, Nt*dt, Nt+1)
    
    ax.plot(timeline, time_results, linestyle='-', linewidth=1.5, color=color,
            label=label)
    
    ax.set_xlabel(r'$t/t_{\mathrm{s}}$', fontsize=fontsize)
    
    ax.legend(loc='best',
              fontsize=12,
              frameon=False,
              ncol=1,
              labelspacing=0.3,
              handlelength=1.5)

def plot_contour(ax, w, dt, Gamma):
    Nt = len(w) - 1
    timeline = np.linspace(0, Nt*dt, Nt+1)
    
    length = np.linspace(0, Gamma, len(w[0]))
    
    m = np.max(np.abs(w))
    
    cf = ax.contourf(timeline, length, w.transpose(), cmap='viridis',
                     levels=np.linspace(-m, m, 50))
    
    pp.colorbar(cf)

    ax.set_xlabel(r'$t/t_{\mathrm{s}}$', fontsize=12)
    ax.set_ylabel(r'$L/D$', fontsize=12)

def plot_profiles(ax, wx, wy):
    Nt = len(wx) - 1

    start = 2*Nt//5    
    step  = Nt//100 if Nt >= 100 else Nt
        
    for n in range(start, Nt+1, step):
        ax.plot(wx[n], wy[n], color='purple', linewidth=0.5)
        
    tip_curve_x = np.concatenate([wx[start:,-1], wx[start-1:start,-1]])
    tip_curve_y = np.concatenate([wy[start:,-1], wy[start-1:start,-1]])
    ax.plot(tip_curve_x, tip_curve_y,
            color='purple', linestyle='--', linewidth=0.5)
    
    ax.set_xlim([-1, 1])
    ax.axis('equal')
