import ngsolve as ngs
from numpy.linalg import norm
from numpy import isnan
from ngsolve.webgui import Draw
from time import time

def solve(fes : ngs.FESpace,                                                        # finite element space
          residual : callable,                                                      # residual(state, test)
          residual_derivative : callable = None,                                    # residual_derivative(state, trial, test) (optional)
          initial_guess :  ngs.GridFunction |  ngs.CoefficientFunction = ngs.CF(0), # initial guess (CoefficientFunction or GridFunction)
          # Inspection parameters
          verbosity : int = 1,                                                      # verbosity level (0 - silent to 3 - detailed)
          draw : bool = False,                                                      # draw intermediate solutions
          # Newton parameters
          maxit_newton : int = 50,             # maximum number of Newton outer iterations
          tol : float = 1e-8,                  # (absolute) tolerance on Newton decrement : sqrt( < residual(uOld), du > )
          rtol_res : float = 1e-10,            # relative tolerance on the residual between 2 iterations (to save 1 useless iteration in case of linear problem)
          # Line search parameters
          linesearch : bool = True,            # flag to enable line search (recommended)
          maxit_linesearch : int = 20,         # maximum iteration number within the line search
          minstep_linesearch : float = 1e-12,  # minimum step size allowed in the line search 
          armijo_linesearch : float = 0.1,     # Armijo coefficient in [0, 1) such that |residual(u-step*du)|² < residual²(u) - armijo_linesearch*step*(|residual(u)|²)'(du)
          step_factor_linesearch : float = 0.3 # step size reduction factor in (0, 1) to reduce the step if too big 
          ) -> dict:
    
    """
    Solve a nonlinear PDE using Newton method.

    Parameters
    ----------
    fes : ngs.FESpace
        The finite element space.

    residual : callable
        Function taking (state, test function) and returning the residual form.

    residual_derivative : callable, optional
        Function taking (state, trial function, test function) and returning
        the bilinear form of the derivative. If None, symbolic differentiation is used.

    initial_guess : ngs.GridFunction or ngs.CoefficientFunction, optional
        Initial solution guess. Default is 0.

    verbosity : int, optional
        Verbosity level (0 = silent, 3 = very detailed). Default is 1.

    draw : bool, optional
        Whether to visualize intermediate results. Default is False.

    maxit_newton : int, optional
        Maximum number of Newton iterations. Default is 50.

    tol : float, optional
        Absolute convergence tolerance on the Newton decrement. Default is 1e-8.

    rtol_res : float, optional
        relative tolerance on the residual between 2 iterations (to save 1 useless 
        iteration in case of linear problem). Default is 1e-10.

    linesearch : bool, optional
        Enable or disable line search. Default is True.

    maxit_linesearch : int, optional
        Maximum number of line search iterations. Default is 20.

    minstep_linesearch : float, optional
        Minimum allowable step size during line search. Default is 1e-12.

    armijo_linesearch : float, optional
        Armijo condition coefficient for line search. Default is 0.1.

    step_factor_linesearch : float, optional
        Multiplicative factor to reduce step size in line search. Default is 0.3.

    Returns
    -------
    results : dict
        A dictionary containing:
        - "solution" : final solution (ngs.GridFunction)
        - "status" : integer code indicating termination reason (see below)
        - "linear_detected" : True if linear problem detected early
        - "iteration" : number of Newton iterations performed
        - "lastInverse" : last tangent matrix decomposition (for reuse or debugging)
        - "residual" : list of residual norms per iteration
        - "decrement" : list of Newton decrement values per iteration
        - "wall_time" : total computation time in seconds

    Status codes:
    -------------
    0 : ✅ SUCCESS — Newton converged successfully.
    1 : ❌ FAILURE — Maximum number of Newton iterations reached.
    2 : ❌ FAILURE — Line search failed: minimum step size reached.
    3 : ❌ FAILURE — Line search failed: max number of iterations reached.
    4 : ❌ FAILURE — NaN encountered in the residual.
    """

    # I) Initialization

    tStart = time()
    if verbosity >= 3 : print(f"-------------------- START NEWTON ---------------------")
    if verbosity >= 3 : print(f"Initializing ... ", end = "")
    du, v = fes.TnT()
    res2 = lambda sol : (norm(ngs.LinearForm(residual(sol, v)).Assemble().vec.FV().NumPy()[fes.FreeDofs()]))**2
    state, state_linesearch, descent = ngs.GridFunction(fes), ngs.GridFunction(fes), ngs.GridFunction(fes)
    state.Set(initial_guess)
    counter_newton = 0
    decrement_list = []
    res2_state = res2(state)
    residual_list = [ngs.sqrt(res2_state)]
    status = 0
    linear = False

    if draw : scene = Draw(state)
    if verbosity >= 3 : print(f"done ({(time()-tStart) * 1000 :.2f} ms).")
    if verbosity >= 2 : print(f"Initial residual : {residual_list[-1] :.5e}")
    if verbosity >= 3 : print(f"Start loop ... ")

    # II) Loop

    while 1:
        counter_newton += 1
        if verbosity >= 2 : print(f" It {counter_newton} -------------------------------------------------")

        # a) Assembly
        tStartAssembly = time()
        if verbosity >= 3 : print(f" - Assembly ... ", end = "")
        res = ngs.LinearForm(residual(state, v)).Assemble()
        if residual_derivative is None : # symbolic differentiation (recommended)
            dres = ngs.BilinearForm(residual(du, v))
            dres.AssembleLinearization(state.vec)
        else :
            dres = ngs.BilinearForm(residual_derivative(state, du, v)).Assemble()
        if verbosity >= 3 : print(f"done ({(time()-tStartAssembly) * 1000 :.2f} ms).")
        tStartSolve = time()
        if verbosity >= 3 : print(f" - Solve ... ", end = "")
        Kinv  = dres.mat.Inverse(freedofs=fes.FreeDofs(), inverse = "sparsecholesky") 
        descent.vec.data = Kinv * res.vec
        if verbosity >= 3 : print(f"done ({(time()-tStartSolve) * 1000 :.2f} ms).")

        decrement_list.append(ngs.sqrt(abs(ngs.Integrate(residual(state, descent), fes.mesh))))

        # b) Line search
        if linesearch :
            tStartLineSearch = time()
            if verbosity >= 2 : print(f" - Line search ... ")
            step = 1.
            counter_linesearch = 0
            state_linesearch.vec.data = state.vec - step * descent.vec
            res2_ls = res2(state_linesearch)
            if verbosity >= 2 : print(f"   it {counter_linesearch} | |residual|² = {res2_ls :.5e}| step = {step : .2e}")

            while not res2_ls < (1-2*armijo_linesearch*step) * res2_state : # enter the line search even if the residual is nan
                step *= step_factor_linesearch
                state_linesearch.vec.data = state.vec - step * descent.vec
                res2_ls = res2(state_linesearch)
                counter_linesearch += 1
                if verbosity >= 2 : print(f"   it {counter_linesearch} | |residual|² = {res2_ls :.5e}| step = {step : .2e}")

                if counter_linesearch >= maxit_linesearch:
                    if verbosity >= 1 : print(f"❌ FAILURE: maximal number of line search iterations reached !!")
                    status = 3
                    break 

                if step < minstep_linesearch:
                    if verbosity >= 1 : print(f"❌ FAILURE: minimal line search step reached !!")
                    status = 2
                    break 
            
            if verbosity >= 3 : print(f" - Line search done ({(time()-tStartLineSearch) * 1000 :.2f} ms).")

            if not status:
                state.vec.data = state_linesearch.vec
            
        else :
            state.vec.data = state.vec - descent.vec

        if isnan(res2_state):
            status = 4
            if verbosity >= 1 : 
                print(f"❌ FAILURE: NaN detected ", end = "")
                if linesearch : print("after line search ", end = "")
            print("!!")
            break

        if status:
            break

        # c) stop criterion
        res2_state = res2(state)
        residual_list.append(ngs.sqrt(res2_state))
        
        if verbosity >= 2 : print(f" - Crit.: |residual| = {residual_list[-1] : .5e}| decr = {decrement_list[-1] :.5e}")
        if draw : scene.Redraw(state)
        if verbosity >= 3 : print(f" - Newton iteration done ({(time()-tStartAssembly) * 1000 :.2f} ms).")


        if residual_list[-1] / residual_list[-2] < rtol_res:
            if verbosity >= 1 : print(f"Linear problem detected!")
            linear = True
            break
        
        if decrement_list[-1] < tol : 
            break

        if counter_newton >= maxit_newton: 
            if verbosity >= 1 : print(f"❌ FAILURE: maximum number of Newton iterations reached !!")
            status = 1
            break
    
    # III) Export results

    if verbosity >=2 and not status : print(f" ✅ SUCCESS: Newton has converged in {counter_newton} iterations.")  
    if verbosity >=2 :  print(f" Total wall time: {(time() - tStart) :.2f} s.")
    results = {"solution" : state, 
               "status" : status, 
               "linear_detected" : linear,
               "iteration": counter_newton, 
               "last_inverse" : Kinv, 
               "residual" : residual_list,
               "decrement": decrement_list,
               "wall_time" : time() - tStart}
    if verbosity >=2 : print(f" --------------------- END NEWTON --------------------- ")  
    return results

def solveAdjoint(state,
                 lastInverse : ngs.Vector = None,
                 rhs : callable = None,
                 expression : callable = None,
                 ) -> dict:
    fes = state.space
    v = fes.TestFunction()
    if lastInverse is not None:
        adjoint = ngs.GridFunction(fes)
        lf = ngs.LinearForm(rhs(state, v)).Assemble()
        if fes.is_complex: 
            adjoint.vec.data = lastInverse.H * lf.vec
        else :
            adjoint.vec.data = lastInverse.T * lf.vec
    elif expression is not None :
        adjoint = solve(fes, expression)
    return adjoint


###############################################################################################################################
# Tests

if __name__ == "__main__" : # simple tests
    from Geometry import transformer

    # geometry and mesh
    mesh = transformer(NCoils = 2)

    # matrial law
    nu0 = 1/(4e-7 * ngs.pi)
    exp5 = lambda x : 1 + x + 0.5*x**2 + (1/6)*x**3 + (1/24)*x**4  + (1/120)*x**5  # replace by ngs.exp(x) to get some nan that linesearch can handle
    dexp5 = lambda x : 1 + x + 0.5*x**2 + (1/6)*x**3 + (1/24)*x**4                 # replace by ngs.exp(x) to get some nan that linesearch can handle
    nu_NL = lambda b2 : 100 + 10 * exp5(1.8 * b2) # from OneLab electrical machine template (https://gitlab.onelab.info/doc/models/-/blob/master/ElectricMachines/BH.pro)
    dnu_NL = lambda b2 : 18 * dexp5(1.8 * b2)
    
    # space, residual and its derivative
    from ngsolve import grad, dx, OuterProduct, H1
    from ngsolve.comp import IntegrationRuleSpace

    FEM_order = 1 # try also with FEM_order > 1 (little bit longer computation then)
    fes = H1(mesh, order = FEM_order, dirichlet = "out")
    # Recall quadrature rule order in non-linear elements to ensure Newton convergence when order > 1
    intrules = IntegrationRuleSpace(fes.mesh, order = fes.globalorder - 1).GetIntegrationRules()

    def residual(state, v):
        jz = 10e6
        b0 = grad(state)
        res = grad(v) * nu0 * b0 * dx("Pp|Pm|Sp|Sm")
        res += grad(v) * nu_NL(b0**2) * b0 * dx("Omega_c", intrules = intrules)
        res +=  v * jz * dx("Pm") - v * jz * dx("Pp")
        return res.Compile()

    def residual_derivative(state, du, v):
        b0 = grad(state)
        dres = grad(v) * nu0 * grad(du) * dx("Pp|Pm|Sp|Sm")
        dres += grad(v) * nu_NL(b0**2) * grad(du) * dx("Omega_c", intrules = intrules)
        dres += grad(v) * (2 * dnu_NL(b0**2) * OuterProduct(b0, b0)* grad(du)) * dx("Omega_c", intrules = intrules)
        return dres.Compile()
    
    # Solve
    res_no_ls = solve(fes, residual, linesearch = False, maxit_newton=100, verbosity = 3) # without line search (sometimes faster but risky!)
    #########################################################################################################################################
    # Recommended:
    res_ls = solve(fes, residual, verbosity = 3) # simpler and sometimes faster (like here where no recompilation of the bilinear form is needed)
    #########################################################################################################################################
    res_ls_no_symd = solve(fes, residual, residual_derivative, verbosity = 3) # provide explicitely the derivative (sometimes faster)

    # Inspect
    print("Total computation time : ")
    print(f"Newton (no line-search)                    : {res_no_ls["wall_time"] :.3f} s")
    print(f"Newton + line search (symbolic derivative) : {res_ls["wall_time"]:.3f} s")
    print(f"Newton + line search (explicit derivative) : {res_ls_no_symd["wall_time"]:.3f} s")

    # Plot
    import matplotlib.pyplot as plt
    plt.figure()
    plt.semilogy(res_no_ls["decrement"], label = "Basic Newton")
    plt.semilogy(res_ls["decrement"], label = "Newton + line search (LS)")
    plt.semilogy(res_ls_no_symd["decrement"], label = "Newton + LS + explicit derivative")

    plt.grid(); plt.legend()
    plt.xlabel("Iterations")
    plt.ylabel("Newton decrement")
    plt.show()