"""
This scripts is similar to dynamic.py. Yet, since we simulate the vortex-
induced load, we need to solve for the new lift variable q, so we consider a
new mixed function space with additional variables.

The additional constants present in this script are:
    - Cl0:
        mean lift coefficient.
    - A, eps:
        coupling constants used in the model of Facchinetti et al. (2004)
"""

import time, datetime

import numpy as np

from fenics import *
# =============================================================================
St = 0.16
bend_to_twist = 1.5

# Aerodynamic coefficients
Cd  = 1.2
Cl0 = 0.3

# Coupling constants
A, eps = 12., 0.3

def initialise_results():
    """
    [t, n, b, w, fext, q, speed]
    """
    return [np.empty(0), np.empty(0), np.empty(0), np.empty(0), np.empty(0), \
            np.empty(0), np.empty(0)]

def run_wake_oscillator(results, mesh, u0, Cy, dt, Nt, Gamma, Ur,
                        speed_function, relaxation):
    
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
    
    element = MixedElement(Ve, Ve, Ve, Ve, Ve, Ve, Ue, Ue)
    V = FunctionSpace(mesh, element)

    Vt = FunctionSpace(mesh, Ve)
    Vs = FunctionSpace(mesh, Ue)

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
    
        t0_i = tangent(coords[i])
        n0_i = normal(coords[i])
           
        Omega0_i = Darboux(coords[i])
        Fint0_i  = internal(coords[i])
  
        w0_i = coords[i]
        v0_i = [0., 0., 0.]
        
        q0_i = 0.001*(2*np.random.rand() - 1)

        S0[v2d[i*20+ 0:i*20+ 3]] = t0_i
        S0[v2d[i*20+ 3:i*20+ 6]] = n0_i
        S0[v2d[i*20+ 6:i*20+ 9]] = Omega0_i 
        S0[v2d[i*20+ 9:i*20+12]] = Fint0_i 
        S0[v2d[i*20+12:i*20+15]] = w0_i
        S0[v2d[i*20+15:i*20+18]] = v0_i
        
        S0[v2d[i*20+18:i*20+19]] = q0_i
        
    # Done. Now we assign the initial values to sol_old
    sol_old.vector()[:] = S0
    
    print('--> Initial configuration assigned.')
    #######################
    
    t , n , Omega , Fint , w , v , q , p  = split(sol)
    t_, n_, Omega_, Fint_, w_, v_, q_, p_ = split(sol_)
    
    t_old, n_old, Omega_old, Fint_old, w_old, v_old, q_old, p_old = split(sol_old)
    
    b     = cross(t,n)
    b_old = cross(t_old,n_old)

    #####################
    # Variational problem
    #
    # This is the part where you can implement the beam model of your choice.
    # Here we consider the Euler-Bernoulli beam theory, so we solve the
    # Kirchhoff equations (see Landau and Lifchitz, Theory of Elasticity)
    #
    # For information, Kirchhoff is German, and its name is pronouced as
    # 'ki' + 'r' + 'sh' (like sh in should) + 'ho' (like ho in hole) + 'ff'.
    
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
       
    class ForceMagnitude(UserExpression):
        def __init__(self, mag, **kwargs):
            self.mag = mag
            
            super().__init__(**kwargs)
    
        def eval(self, values, x):
            if not(np.array_equal(x,coords[0])): # ie. avoid the root
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
            
    speed     = Speed(t_n=0)
    speed_old = SpeedOld(t_n=0)

    lmbda = Constant(St*Gamma/Ur)
    u_rel = speed*u0_ufl - lmbda*v
    
    ut = dot(t, u_rel)
    un = dot(n, u_rel)
    ub = dot(b, u_rel)

    uperp  = as_vector(normal(coords[0]))
    
    u_rel_old = speed_old*u0_ufl - lmbda*v_old
    
    ub_old = dot(b_old, u_rel_old)
    un_old = dot(n_old, u_rel_old)

    # Drag
    drag_mag = ForceMagnitude(Cd*Cy, degree=3)
    f_drag   = drag_mag*(un*abs(un)*n + ub*abs(ub)*b)

    term1_b = lmbda*(ub*b - ub_old*b_old)/dt
    term2_b = tgrad(project(ub*ut,Vs)*b)
    term3_b = Constant(1./2)*tgrad(project(ub*ub,Vs)*t)
    
    term1_n = lmbda*(un*n - un_old*n_old)/dt
    term2_n = tgrad(project(un*ut,Vs)*n)
    term3_n = Constant(1./2)*tgrad(project(un*un,Vs)*t)

    # Added mass force
    added_mag = ForceMagnitude(np.pi*Cy/Gamma/2., degree=3)
    f_added = added_mag*(  term1_b + term1_n \
                         + term2_b + term2_n \
                         - term3_b - term3_n) 

    # Vortex-induced lift
    viv_mag = ForceMagnitude(0.5*Cl0*Cy, degree=3)
    f_viv   = viv_mag*ub*ub*q*uperp
       
    fext = f_drag + f_viv + f_added
    
    # dFint/ds = m*acc - fext
    dFds = tgrad(Fint)
    L_dFds = dot(dFds + fext - (v - v_old)/dt, Fint_)
    
    L_dwdt = dot((w - w_old)/dt - v, v_)
    
    # van der Pol equations
    omegaf = Constant(2*np.pi*Ur)
    
    L_dqdt = (p - (q - q_old)/dt)*q_
    
    L_dpdt = ((p - p_old)/dt + \
              + Constant(eps)*omegaf*abs(ub)*(q*q - 1)*p \
              + omegaf*omegaf*ub*ub*q \
              - Constant(A*Gamma)*dot((v-v_old)/dt,uperp) )*p_
               
    # Final form
    L = L_dtds + L_dnds + L_dMds + L_dFds + L_dwds + L_dwdt + L_dqdt + L_dpdt
    F = L*dx(degree=1)
    
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

    bc_q = DirichletBC(V.sub(6), Constant(0.), root)
    bc_p = DirichletBC(V.sub(7), Constant(0.), root)
                       
    bc = [bc_t, bc_n, bc_Omega, bc_Fint, bc_w, bc_v, bc_q, bc_p]
    
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
    
    results[6] = np.array([[sol_old(e)[-2] for e in coords]])
    
    fext_proj = project(fext, Vt)
    results[4] = np.array([[fext_proj(e) for e in coords]])
    
    speed_proj = project(speed, Vs)
    results[5] = np.array([[speed_proj(e) for e in coords]])
    ###########################################################
    
    ##################
    # Ready to solve #
    start = time.time()

    t_n = 0
    print('Solving...')
    for n in range(Nt):
        t_n += dt
    
        print('t_n = %e' % t_n)
    
        speed.t_n     = t_n
        speed_old.t_n = t_n

        solver.solve()
    
        # Results
        T = sol.sub(0, True)
        N = sol.sub(1, True)
        W = sol.sub(4, True)
        
        Q = sol.sub(6, True)
    
        current_t = np.array([[T(e) for e in coords]])
        current_n = np.array([[N(e) for e in coords]])
        current_b = np.array([[np.cross(T(e),N(e)) for e in coords]])
    
        current_w = np.array([[W(e) for e in coords]])
        
        current_q = np.array([[Q(e) for e in coords]])
        
        fext_proj = project(fext, Vt)
        current_fext = np.array([[fext_proj(e) for e in coords]])
       
        speed_proj = project(speed, Vs)
        current_speed = np.array([[speed_proj(e) for e in coords]])

        results[0] = np.vstack([results[0], current_t])
        results[1] = np.vstack([results[1], current_n])
        results[2] = np.vstack([results[2], current_b])
    
        results[3] = np.vstack([results[3], current_w])
        
        results[6] = np.vstack([results[6], current_q])
        
        results[4] = np.vstack([results[4], current_fext])
        
        results[5] = np.vstack([results[5], current_speed])
        
        sol_old.assign(sol) # DON'T FORGET ABOUT THIS LINE!
    ##################
        
    print('Done')
    cpu_time = time.time() - start
    print('cpu_time = %s seconds = %s (hh:mm:ss)' %\
          (cpu_time, datetime.timedelta(seconds=cpu_time)))
