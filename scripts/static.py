"""
This script simulates the static deformation of a rod using the FEniCS solver.

The constant bend_to_twist is the ratio of the flexural rigidity to the
torsional rigidity, ie. EI/GJ

"""

import numpy as np

from fenics import *
#==============================================================================
bend_to_twist = 1.5

def run_static(mesh, u0, which_force, force_mag, relaxation=0.8):
    cell = mesh.ufl_cell()

    coords   = mesh.coordinates()
    n_coords = mesh.num_vertices()

    facet_coords = np.array([[Facet(mesh,k).normal(i) for i in range(3)] \
                              for k in range(n_coords)])
    facet_coords[0] *= -1
    
    # Domains
    def root(x, on_boundary):
        return near(x[0], coords[0,0]) \
           and near(x[1], coords[0,1]) \
           and near(x[2], coords[0,2])
    
    def tip(x, on_boundary):
        return near(x[0], coords[-1,0]) \
           and near(x[1], coords[-1,1]) \
           and near(x[2], coords[-1,2])

    #####################################
    # Vectorial and mixed function spaces
    Ve = VectorElement('CG', cell, degree=1)
    
    element = MixedElement(Ve, Ve, Ve, Ve, Ve)
    V = FunctionSpace(mesh, element)

    Vt = FunctionSpace(mesh, Ve)

    v2d = vertex_to_dof_map(V)

    v2dt = vertex_to_dof_map(Vt)
    #####################################

    ############################################
    # Initial material frame and other functions
    tangent  = Function(Vt)
    normal   = Function(Vt)
    binormal = Function(Vt)
    Darboux  = Function(Vt)
    internal = Function(Vt)

    Tangent = tangent.vector()[:]
    for v in vertices(mesh):
        i = v.index()
        Tangent[v2dt[i*3:(i+1)*3]] = facet_coords[i]    
    tangent.vector()[:] = Tangent
    
    def tgrad(u):
        return dot(grad(u),tangent)
    
    u0_ufl = as_vector(u0)

    normal_ufl  = cross(tangent, u0_ufl)
    normal_ufl /= sqrt(dot(normal_ufl, normal_ufl))
    normal.assign(project(normal_ufl,Vt))

    binormal_ufl  = cross(tangent, normal)
    binormal_ufl /= sqrt(dot(binormal_ufl, binormal_ufl))
    binormal.assign(project(binormal_ufl,Vt))

    dtangentds = tgrad(tangent)
    dnormalds  = tgrad(normal)

    Darboux_ufl = cross(tangent, dtangentds) + dot(binormal,dnormalds)*tangent
    Darboux.assign(project(Darboux_ufl,Vt))
    
    Moment = Darboux + (1./bend_to_twist - 1)*dot(Darboux, tangent)*tangent
    dMomentds = tgrad(Moment)
    
    internal_ufl = cross(tangent, dMomentds)
    internal.assign(project(internal_ufl, Vt))

    print('--> Initial material frame defined.')
    ############################################

    ######################
    # Function declaration
    sol  = Function(V)
    sol_ = TestFunction(V)
    ######################

    #######################
    # Initial configuration
    S = sol.vector()[:]
    
    for v in vertices(mesh):
        i = v.index()
    
        t0_i     = tangent(coords[i])
        n0_i     = normal(coords[i])
        
        Omega0_i = Darboux(coords[i])
        Fint0_i  = internal(coords[i])
   
        w0_i = coords[i]
    
        S[v2d[i*15+ 0:i*15+ 3]] = t0_i
        S[v2d[i*15+ 3:i*15+ 6]] = n0_i
        S[v2d[i*15+ 6:i*15+ 9]] = Omega0_i
        S[v2d[i*15+ 9:i*15+12]] = Fint0_i
        S[v2d[i*15+12:i*15+15]] = w0_i
    
    # Done. Now we assign the initial values to sol
    sol.vector()[:] = S
    
    print('--> Initial configuration assigned.')
    #######################
    
    t , n , Omega , Fint , w  = split(sol)
    t_, n_, Omega_, Fint_, w_ = split(sol_)
    
    b = cross(t,n)
    
    #####################
    # Variational problem
    #
    # This is the part where you can implement the beam model of your choice.
    # Here we consider the Euler-Bernoulli beam theory, so we solve the
    # Kirchhoff equations (see Landau and Lifchitz, Theory of Elasticity)
    #
    # For information, Kirchhoff is German, and its name is pronouced as
    # 'ki' + 'r' + 'sh' (like 'sh' in should) + 'ho' (like 'ho' in hall) + 'ff'.
    
    # Equation dw/ds = t
    dwds = tgrad(w)
    L_dwds = dot(dwds - t, w_)
    
    # d?ds = Omega x ?, with ? = t, n
    dtds = tgrad(t)
    L_dtds = dot(dtds - cross(Omega, t), t_)
    
    dnds = tgrad(n)
    L_dnds = dot(dnds - cross(Omega, n), n_)
    
    # Constitutive equation M = Omega + (1./bend_to_twist - 1)*(Omega.t)t
    M = Omega + (1./bend_to_twist - 1)*dot(Omega, t)*t
    
    # dM = Fint x t
    dMds = tgrad(M)
    L_dMds = dot(dMds - cross(Fint, t), Omega_)
       
    # dFint/ds = -fext
    dFds = tgrad(Fint)
    
    class ForceMagnitude(UserExpression):
        def __init__(self, mag, **kwargs):
            self.mag = mag
    
            super().__init__(**kwargs)
    
        def eval(self, values, x):
            if not(np.array_equal(x,coords[0])): # ie. avoid the root
                values[0] = self.mag
            else:
                values[0] = 0
        
    if which_force == 'distributed':
        distr_mag = ForceMagnitude(force_mag, degree=3)
        fext      = distr_mag*u0_ufl
        
    if which_force == 'drag':
        un = np.dot(n, u0_ufl)
        ub = np.dot(b, u0_ufl)
        
        drag_mag = ForceMagnitude(force_mag, degree=3)
        fext     = drag_mag*(un*abs(un)*n + ub*abs(ub)*b)
        
    L_dFds = dot(dFds + fext, Fint_)
    
    # Final form
    F = (L_dtds + L_dnds + L_dMds + L_dFds + L_dwds)*dx(degree=1)
    
    print('--> Variational problem defined.')
    #####################

    #######################
    # Boundary conditions #
    def root(x, on_boundary):
        return near(x[0], coords[0,0]) \
           and near(x[1], coords[0,1]) \
           and near(x[2], coords[0,2])
    
    def tip(x, on_boundary):
        return near(x[0], coords[-1,0]) \
           and near(x[1], coords[-1,1]) \
           and near(x[2], coords[-1,2])

    t_root = tangent(coords[0])
    n_root = normal(coords[0])
    
    bc_t = DirichletBC(V.sub(0), Constant((t_root[0], t_root[1], t_root[2])),
                       root)
    bc_n = DirichletBC(V.sub(1), Constant((n_root[0], n_root[1], n_root[2])),
                       root)
    
    bc_Omega = DirichletBC(V.sub(2), Constant((0., 0., 0.)), tip)
    bc_Fint  = DirichletBC(V.sub(3), Constant((0., 0., 0.)), tip)    
    
    bc_w = DirichletBC(V.sub(4), Constant((0., 0., 0.)), root)
    
    bc = [bc_t, bc_n, bc_Omega, bc_Fint, bc_w]
    
    print('--> Boundary conditions defined.')
    #######################

    #########################
    # Solver initialisation #
    Jac     = derivative(F, sol)
    problem = NonlinearVariationalProblem(F, sol, bc, Jac)
    solver  = NonlinearVariationalSolver(problem)
    
    solver.parameters['newton_solver']['maximum_iterations']   = 50
    solver.parameters['newton_solver']['relative_tolerance']   = 1e-6
    solver.parameters['newton_solver']['absolute_tolerance']   = 1e-8
    solver.parameters['newton_solver']['relaxation_parameter'] = relaxation
    #########################

    ##################
    # Ready to solve #
    print('Solving...')
    
    solver.solve()
    
    print('Done')
    ##################

    ####################
    # Results extraction
    T = sol.sub(0, True)
    N = sol.sub(1, True)
    W = sol.sub(4, True)

    results_T = np.array([T(e) for e in coords])
    results_N = np.array([N(e) for e in coords])
    results_B = np.array([np.cross(T(e),N(e)) for e in coords])

    results_W = np.array([W(e) for e in coords])

    fext_proj = project(fext, Vt)
    results_Fext = np.array([fext_proj(e) for e in coords])
    ####################
    
    return np.array([results_T, results_N, results_B, results_W, results_Fext])
