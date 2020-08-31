"""
This script simulates the beam dynamics when put under hydrodynamic loads using
the FEniCS solver.

It contains three functions:

    initialise_results
    ------------------
        This is a way to save the results. Here I decided to save the vectors
        of the material frame (t, n, b, ), the displacement w, the external
        load fext, and the speed function U(t)/U0. Please feel free to modify
        the saving process as you like.
        
    run_dynamic
    -----------
        The core function of this script. Description is provided below.
        
The constants defined here are:
    - Cd:
        drag coefficient
    - bend_to_twist:
        the ratio of the flexural rigidity to the torsional rigidity EI/GJ
    - St:
        Strouhal number

"""

import time, datetime

import numpy as np

from fenics import *
# =============================================================================
Cd = 1.2
bend_to_twist = 1.5

St = 0.16

def initialise_results():
    return [np.empty(0), np.empty(0), np.empty(0), np.empty(0), np.empty(0),
            np.empty(0)]

def run_dynamic(results, mesh, u0, which_force, force_mag, dt, Nt, Ur, Gamma,
                speed_function, relaxation):
    """
    This function simulates the beam dynamics given a set of parameters.
    
    Parameters
    ----------
    results:
        the initialised array with empty subarrays of the results.
    
    mesh:
        a Mesh object.
        
    u0:
        unit vector of the flow stream direction.
        
    which_force:
        force considered (useful when dealing with many models).
        
    force_mag:
        the main force_mag involved in the hydrodynamic load.
        e.g. for the distributed load it was the coefficient g0 = f0L^3/EI.
        
    dt:
        time step.
        
    Nt:
        number of time steps.
        
    Ur:
        reduced velocity.
        
    Gamma:
        aspect ratio of the beam.
        
    relaxation:
        relaxation parameter used in the Newton solver. If the beam tends to be
        considerably curved, try to reduce this parameter.
        
    Returns
    -------
    Nothing! Because the initialised array of results is a Python list of numpy
    arrays, the output is automatically updated (results before simulations:
    empty ---> results after simulations: full of results). In this way, even
    if the simulation crashes in middle way, you can display your results and
    troubleshoot your problem.
    
    """
    
    cell = mesh.ufl_cell()
    
    coords   = mesh.coordinates()
    n_coords = mesh.num_vertices()

    facet_coords = np.array([[Facet(mesh,k).normal(i) for i in range(3)] \
                              for k in range(n_coords)])
    facet_coords[0] *= -1
        
    ##############################################
    # Scalar, vectorial, and mixed function spaces
    Ue = FiniteElement('CG', cell, degree=1)
    Ve = VectorElement('CG', cell, degree=1)
    
    element = MixedElement(Ve, Ve, Ve, Ve, Ve, Ve)
    V = FunctionSpace(mesh, element)

    Vs = FunctionSpace(mesh, Ue)
    Vt = FunctionSpace(mesh, Ve)

    v2d = vertex_to_dof_map(V)

    v2dt = vertex_to_dof_map(Vt)
    ##############################################

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
    
    sol_old = Function(V)

    print('--> Functions declared.')
    ######################

    #######################
    # Initial configuration
    S0 = sol_old.vector()[:]
    
    for v in vertices(mesh):
        i = v.index()
        
        t0_i     = tangent(coords[i])
        n0_i     = normal(coords[i])
        
        Omega0_i = Darboux(coords[i])
        Fint0_i  = internal(coords[i])
   
        w0_i = coords[i]
        v0_i = [0., 0., .0] # zero initial velocity
   
        S0[v2d[i*18+ 0:i*18+ 3]] = t0_i
        S0[v2d[i*18+ 3:i*18+ 6]] = n0_i
        S0[v2d[i*18+ 6:i*18+ 9]] = Omega0_i
        S0[v2d[i*18+ 9:i*18+12]] = Fint0_i
        S0[v2d[i*18+12:i*18+15]] = w0_i
        S0[v2d[i*18+15:i*18+18]] = v0_i
        
    # Done. Now we assign the initial values to sol_old
    sol_old.vector()[:] = S0
    
    print('--> Initial configuration assigned.')
    #######################
    
    t , n , Omega , Fint , w , v  = split(sol)
    t_, n_, Omega_, Fint_, w_, v_ = split(sol_)
    
    t_old, n_old, Omega_old, Fint_old, w_old, v_old = split(sol_old)
    
    b     = cross(t,n)
    b_old = cross(t_old, n_old)
    
    #####################
    # Variational problem
    #
    # This is the part where you can implement the beam model of your choice.
    # Here we consider the Euler-Bernoulli beam theory, so we solve the
    # Kirchhoff equations (see Landau and Lifchitz, Theory of Elasticity)
    #
    # For information, Kirchhoff is German, and its name is pronouced as
    # 'ki' + 'r' + 'sh' (like 'sh' in should) + 'ho' (like 'h'o in hall) + 'ff'.
    
    # Equation dw/ds = t
    dwds = tgrad(w)
    L_dwds = dot(dwds - t, w_)
    
    # d?ds = Omega x ?, with ? = t, n
    dtds = tgrad(t)
    L_dtds = dot(dtds - cross(Omega, t), t_)
    
    dnds = tgrad(n)
    L_dnds = dot(dnds - cross(Omega, n), n_)
    
    # Constitutive equation M = Omega + (1./gamma - 1)*(Omega.t)t
    M = Omega + (1./bend_to_twist - 1)*dot(Omega, t)*t
    
    # dM = Fint x t
    dMds = tgrad(M)
    L_dMds = dot(dMds - cross(Fint, t), Omega_)
    
    # Definition of UFL expressions by subclassing.
    class ForceMagnitude(UserExpression):
        def __init__(self, mag, **kwargs):
            self.mag = mag
            
            super().__init__(**kwargs)
    
        def eval(self, values, x):
            if not(np.array_equal(x,coords[0])): # ie. avoid the root.
                values[0] = self.mag
            else:
                values[0] = 0

    class Speed(UserExpression):
        def __init__(self, t_n, **kwargs):
            self.t_n  = t_n

            super().__init__(**kwargs)

        def eval(self, values, x):
            speed_up = speed_function(self.t_n)
            values[0] = speed_up

    class SpeedOld(UserExpression):
        def __init__(self, t_n, **kwargs):
            self.t_n  = t_n

            super().__init__(**kwargs)

        def eval(self, values, x):
            speed_up_old = speed_function(self.t_n - dt)
            values[0] = speed_up_old

    # Here comes the definition of your forces. We restricted ourselves to the
    # cases of a distributed force, and the drag with the added mass force.
    # You can think about gravity/buoyancy forces, electromagnetic forces, etc.
    speed     = Speed(t_n=0)
    speed_old = SpeedOld(t_n=0)
    if which_force == 'distributed':
        distr_mag = ForceMagnitude(force_mag, degree=3)
        fext      = speed*distr_mag*u0_ufl

    if which_force == 'drag':
        lmbda = Constant(St*Gamma/Ur)
        u_rel = speed*u0_ufl - lmbda*v
        
        ut = dot(t, u_rel)
        un = dot(n, u_rel)
        ub = dot(b, u_rel)
    
        u_rel_old = speed_old*u0_ufl - lmbda*v_old
        
        ub_old = dot(b_old, u_rel_old)
        
        drag_mag = ForceMagnitude(force_mag, degree=3)
        f_drag   = drag_mag*(un*abs(un)*n + ub*abs(ub)*b)

        term1_b = lmbda*(ub*b - ub_old*b_old)/dt
        term2_b = tgrad(project(ub*ut,Vs)*b)
        term3_b = Constant(1./2)*tgrad(project(ub*ub,Vs)*t)

        Cy = force_mag/Cd
        added_mag = ForceMagnitude(np.pi*Cy/Gamma/2., degree=3)
        f_added = added_mag*(term1_b + term2_b - term3_b)

        fext = f_drag + f_added

    # dFint/ds = acc - fext and vel = wdot
    dFds = tgrad(Fint)
    L_dFds = dot(dFds + fext - (v - v_old)/dt, Fint_)
    
    L_dwdt = dot((w - w_old)/dt - v, v_)
    
    # Final form
    F = (L_dwds + L_dtds + L_dnds + L_dMds + L_dFds + L_dwdt)*dx(degree=1)
    
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
    bc_v = DirichletBC(V.sub(5), Constant((0., 0., 0.)), root)
    
    bc = [bc_t, bc_n, bc_Omega, bc_Fint, bc_w, bc_v]
    
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

    # Not to display the information about Newton convergence.
    set_log_active(False)

    ###########################################################
    # Initialising the results with the initial configuration #
    results[0] = np.array([[sol_old(e)[0:3] for e in coords]])
    results[1] = np.array([[sol_old(e)[3:6] for e in coords]])
    results[2] = np.array([np.cross(results[0][0], results[1][0])])
    
    results[3] = np.array([[sol_old(e)[12:15] for e in coords]])
    
    fext_proj = project(fext, Vt)
    results[4] = np.array([[fext_proj(e) for e in coords]])

    speed_proj = project(speed, Vs)
    results[5] = np.array([[speed_proj(e) for e in coords]])
    ###########################################################
    
    ##################
    # Ready to solve #
    t_n = 0
    print('Solving...')

    start = time.time()
    for n in range(Nt):
        t_n += dt
    
        print('t_n = %e' % t_n)
    
        speed.t_n     = t_n
        speed_old.t_n = t_n
        
        solver.solve()
    
        # Extracting subfunctions
        T = sol.sub(0, True)
        N = sol.sub(1, True)
        W = sol.sub(4, True)
    
        current_t = np.array([[T(e) for e in coords]])
        current_n = np.array([[N(e) for e in coords]])
        current_b = np.array([[np.cross(T(e),N(e)) for e in coords]])
    
        current_w = np.array([[W(e) for e in coords]])
        
        fext_proj = project(fext, Vt)
        current_fext = np.array([[fext_proj(e) for e in coords]])

        speed_proj = project(speed, Vs)
        current_speed = np.array([[speed_proj(e) for e in coords]])
        
        results[0] = np.vstack([results[0], current_t])
        results[1] = np.vstack([results[1], current_n])
        results[2] = np.vstack([results[2], current_b])
    
        results[3] = np.vstack([results[3], current_w])
        
        results[4] = np.vstack([results[4], current_fext])

        results[5] = np.vstack([results[5], current_speed])

        sol_old.assign(sol) # DON'T FORGET ABOUT THIS LINE!
    ##################
        
    print('Done')
    
    cpu_time = time.time() - start
    print('cpu_time = %s seconds = %s (hh:mm:ss)' %\
          (cpu_time, datetime.timedelta(seconds=cpu_time)))

