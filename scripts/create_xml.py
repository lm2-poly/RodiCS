import numpy as np
# =============================================================================
n_elements = 300
n_vertices = n_elements + 1

L = 1

list_x = np.linspace(0, 0, n_vertices)
list_y = np.linspace(0, L, n_vertices)
list_z = np.linspace(0, 0, n_vertices)

zipped = zip(list_x, list_y, list_z)

def save_xml():
    outfile = r'../xml_files/beam_mode1_%s.xml' % n_elements
    
    f = open(outfile, 'w')
    
    f.write(r'<dolfin xmlns:dolfin="http://www.fenicsproject.org">' + '\n')
    f.write(2*' ' + r'<mesh celltype="interval" dim="3">' + '\n')
    
    # Writing vertices
    f.write(4*' ' + r'<vertices size="%d">' % n_vertices  + '\n')
    for i, coords in enumerate(zipped):
        x, y, z = coords
        f.write(6*' ' + r'<vertex index="%d" x="%.6e" y="%.6e" z="%.6e"/>' \
                % (i, x, y, z) + '\n')
    f.write(4*' ' + r'</vertices>' + '\n')
    
    # Writing elements (called cells in FEniCS)
    f.write(4*' ' + r'<cells size="%d">' % n_elements + '\n')
    for k in range(n_elements):
        f.write(6*' ' + r'<line index="%d" v0="%d" v1="%d"/>' \
                % (k, k, k+1) + '\n')
    f.write(4*' ' + r'</cells>' + '\n')
    
    f.write(2*' ' + r'</mesh>' + '\n')
    f.write(r'</dolfin>' + '\n')
    
    f.close()
    
save_xml()
