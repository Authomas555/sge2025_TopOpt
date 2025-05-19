from ngsolve import *


R = CF(((0,1),(-1,0)), dims = (2,2))

def Curl(u):
    return R * grad(u)


def solveMag(nu, mesh, jz = 1e6):
    
    # 1) définition de l'espace fonctionnel
    fes = H1(mesh, order=1, dirichlet="out") 
    u = fes.TrialFunction()
    u_star = fes.TestFunction()
    
    # 2) Définition de la forme bilinéaire
    K = BilinearForm(fes, symmetric=True)
    K += Curl(u_star)* (nu*Curl(u)) *dx 

    # 3) Définition de la forme linéaire
    s = LinearForm(fes)
    s += u_star*jz*dx("Pp") - u_star*jz*dx("Pm")
    
    # 4) Assemblage -> on sort du monde continu pour discrétiser
    K.Assemble() # K est désormais une matrice
    s.Assemble() # s est désormais un vecteur
    
    # 5) Résolution -> la solution formelle s'écrit a = K^-1 * s 
    # (en réalité l'inverse n'est pas explicitement calculée, on utilise une décomposition de Cholesky)
    
    state = GridFunction(fes)  # On déclare que la solution est une fonction définie point par point.
    Kinv = K.mat.Inverse(inverse="sparsecholesky", freedofs=fes.FreeDofs())
    state.vec.data = Kinv * s.vec
    return state, Kinv


def solveAdjoint(Kinv, u, df):
    
    fes = u.space 
    u_star = fes.TestFunction()

    s = LinearForm(fes)
    s += df(u, u_star) 

    s.Assemble()
    
    adjointState = GridFunction(fes) 
    adjointState.vec.data = Kinv.T * s.vec
    return adjointState
