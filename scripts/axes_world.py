import numpy as np
import matplotlib.pyplot as pp

from matplotlib import rc
#==============================================================================
rc('text', usetex=True)

markersize = 8
fontsize   = 12

def one_by_one(inch_x=3.5, inch_y=3.5):
    fig = pp.figure(figsize=(inch_x,inch_y))
    ax = fig.add_axes([0.15,0.13,0.8,0.81])
    
    pp.setp(ax.spines.values(), linewidth=0.7)

    ax.tick_params(axis='both', which='major', width=0.7)

    return ax

def one_by_two(inch_x=6.5,inch_y=3.5):
    fig = pp.figure(figsize=(inch_x,inch_y))
    
    ax_a = fig.add_axes([0.11,0.14,0.37,0.8])
    
    ax_b = fig.add_axes([0.61,0.14,0.37,0.8])
    
    ax_a.text(-0.26, 0.95, '$(a)$', fontsize=fontsize, transform=ax_a.transAxes)
    ax_b.text(-0.26, 0.95, '$(b)$', fontsize=fontsize, transform=ax_b.transAxes)
    
    pp.setp(ax_a.spines.values(), linewidth=0.7)
    pp.setp(ax_b.spines.values(), linewidth=0.7)
    
    ax_a.tick_params(axis='both', which='major', width=0.7)
    ax_b.tick_params(axis='both', which='major', width=0.7)
    
    return ax_a, ax_b

def two_by_one(inch_x=5, inch_y=3.5*2):
    fig = pp.figure(figsize=(inch_x,inch_y))
    
    ax_b = fig.add_axes([0.14,0.07,0.83,0.41])
    
    ax_a = fig.add_axes([0.14,0.56,0.83,0.41])
    
    ax_b.text(-0.13, 0.95, '$(b)$', fontsize=12, transform=ax_b.transAxes)
    ax_a.text(-0.13, 0.95, '$(a)$', fontsize=12, transform=ax_a.transAxes)
    
    pp.setp(ax_a.spines.values(), linewidth=0.7)
    pp.setp(ax_b.spines.values(), linewidth=0.7)
    
    ax_a.tick_params(axis='both', which='major', width=0.7)
    ax_b.tick_params(axis='both', which='major', width=0.7)
    
    return ax_a, ax_b

def slope_triangle(x0, y0, dx, dy, label_x, label_y, up, ax):
    x1 = x0*np.exp(dx)
    y1 = y0*np.exp(dy)

    if up == 'Up':
        triangle_x = [x0, x1, x0]
        triangle_y = [y0, y1, y1]

        text_x = x0*0.70, np.sqrt(y0*y1)
        text_y = np.sqrt(x0*x1), y1*1.05

        va_x = 'center'
        ha_x = 'left'

        va_y = 'bottom'
        ha_y = 'center'

    else:
        triangle_x = [x0, x1, x1]
        triangle_y = [y0, y0, y1]

        text_x = np.sqrt(x0*x1), y0*0.9
        text_y = x1*1.2, np.sqrt(y0*y1)

        va_x = 'top'
        ha_x = 'center'

        va_y = 'center'
        ha_y = 'left'

    ax.fill(triangle_x, triangle_y,
            edgecolor='dimgrey',
            facecolor='lightgrey',
            alpha=0.25)

    ax.text(text_x[0], text_x[1], r'$%s$' % label_x,
            verticalalignment=va_x,
            horizontalalignment=ha_x,
            fontsize=12)

    ax.text(text_y[0], text_y[1], r'$%s$' % label_y,
            verticalalignment=va_y,
            horizontalalignment=ha_y,
            fontsize=12)

