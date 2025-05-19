import matplotlib.pyplot as plt
import numpy as np
import sys, os
from ngsolve import *
from netgen.geom2d import SplineGeometry
from levelBaseLib import *
from solverlib import *
from advect_lib import *
def transformer(NCoils = 1, width = 0.6, height = 0.6, innerDiameterCoil = 4e-2, outerDiameterCoil = 6e-2, heightCoil = 4e-2, hmax = 1e-2):
    ## Geometry generation

    geo = SplineGeometry()
    geo.AddRectangle( p1 = (-width/2, -height/2), p2 = (width/2, height/2), leftdomain = 1, rightdomain = 0, bc = "out")

    if NCoils == 1:
        geo.AddRectangle( p1 = (-outerDiameterCoil/2, -heightCoil/2), p2 = (-innerDiameterCoil/2, heightCoil/2),
        leftdomain = 2, rightdomain = 1 )
        geo.AddRectangle( p1 = (innerDiameterCoil/2, -heightCoil/2), p2 = (outerDiameterCoil/2, heightCoil/2),
        leftdomain = 3, rightdomain = 1 )
        geo.SetMaterial( 2, "S1",bc='coil' )      # conducteur primaire positif
        geo.SetMaterial( 3, "S2",bc='coil' ) 
    elif NCoils == 2:
        th = (outerDiameterCoil - innerDiameterCoil) / 2
        ep = (outerDiameterCoil + innerDiameterCoil) / 2
        geo.AddRectangle( p1 = (-width/2 + (0.3-0.11), -heightCoil/2), p2 = (-width/2 + (0.3-0.11) + th, heightCoil/2),
        leftdomain = 2, rightdomain = 1,bc='coil' )
        geo.AddRectangle( p1 = (-width/2 + (0.3-0.11) + ep, -heightCoil/2), p2 = (-width/2 + (0.3-0.11) + ep + th, heightCoil/2),
        leftdomain = 3, rightdomain = 1,bc='coil' )
        geo.AddRectangle( p1 = (width/2 - (0.3-0.11) - th - ep, -heightCoil/2), p2 = (width/2 - (0.3-0.11) - ep, heightCoil/2),
        leftdomain = 4, rightdomain = 1,bc='coil' )
        geo.AddRectangle( p1 = (width/2 - (0.3-0.11) - th, -heightCoil/2), p2 = (width/2 - (0.3-0.11) , heightCoil/2),
        leftdomain = 5, rightdomain = 1,bc='coil' )
        geo.SetMaterial( 2, "P1" )      # conducteur primaire positif
        geo.SetMaterial( 3, "P2" )      # conducteur primaire négatif
        geo.SetMaterial( 4, "S1" )      # conducteur secondaire positif
        geo.SetMaterial( 5, "S2" )      # conducteur secondaire négatif

    # On nomme les différentes zones afin de les repérer plus tard :
    geo.SetMaterial( 1, "D" )       # zone de design
    ngmesh = geo.GenerateMesh( maxh = hmax ) # Maillage avec NetGen
    mesh = Mesh( ngmesh )    # Importation du maillage NetGen dans NgSolve
    return mesh
def eqFluxMoyen(mesh,state):
    return state*dx(mesh.Materials("S2")) - state*dx(mesh.Materials("S1")) 
def fluxMoyen(mesh,state):
    return Integrate (eqFluxMoyen(mesh,state),mesh) 
def duFluxMoyen(mesh,v):
    return v*dx(mesh.Materials("S2")) - v*dx(mesh.Materials("S1"))
def fObjectif(mesh,state):
    return  -1*fluxMoyen(mesh,state)
def dufObjectif(mesh,state,v,phi):
    return  -1*duFluxMoyen(mesh,v)
def eqObjectif(mesh,state):
    return -1*eqFluxMoyen(mesh,state)
def Lagrangian(solver,phi):
    return eqObjectif(solver.mesh,solver.gfu) + solver.getEquation(phi)
def ShapeDeriv(solver,phi):
    VEC = H1(mesh, order=1, dim=2,dirichlet='coil')
    PHI, X = VEC.TnT()
    dJOmega = LinearForm(VEC)
    dJOmega += Lagrangian(solver,phi).DiffShape(X) #Diff auto
    # Regularization
    b = BilinearForm(VEC)
    b += InnerProduct((0.6/100)**2*grad(X),grad(PHI))*dx + InnerProduct(X,PHI)*dx
    b.Assemble()
    dJOmega.Assemble()
    ########## rajout chelou
    #gfu = GridFunction(VEC)
    #gfu.vec.data = dJOmega.vec
    #gfu2 = GridFunction(VEC)
    #gfu2.Set(gfu / Norm(gfu) * atan(1e6*Norm(gfu) ))
    #dJOmega.vec.data = gfu2.vec
    ##########
    gfX = GridFunction(VEC)
    rhs = gfX.vec.CreateVector()
    rhs.data = dJOmega.vec - b.mat * gfX.vec
    update = gfX.vec.CreateVector()
    update.data = b.mat.Inverse(VEC.FreeDofs()) * rhs
    gfX.vec.data += update
    return gfX
def VolFrac(mesh,phi):
    return 1/TotalVolume * Integrate(InterpolateLevelSetToElems(phi,1,0,mesh)*dx(mesh.Materials("D")),mesh)
def ShapeDerivVolFrac(mesh,phi):
    VEC = H1(mesh, order=1, dim=2,dirichlet='coil')
    PHI, X = VEC.TnT()
    dJOmega = LinearForm(VEC)
    dJOmega += (1/TotalVolume * InterpolateLevelSetToElems(phi,1,0,mesh)*dx(mesh.Materials("D"))).DiffShape(X) #Diff auto
    # Regularization
    b = BilinearForm(VEC)
    b += InnerProduct((0.6/100)**2*grad(X),grad(PHI))*dx + InnerProduct(X,PHI)*dx
    b.Assemble()
    dJOmega.Assemble()
    gfX = GridFunction(VEC)
    rhs = gfX.vec.CreateVector()
    rhs.data = dJOmega.vec - b.mat * gfX.vec
    update = gfX.vec.CreateVector()
    update.data = b.mat.Inverse(VEC.FreeDofs()) * rhs
    gfX.vec.data += update
    return gfX

# Material properties
mu0 = 4e-7 * np.pi
mur = 1000
nuiron = 1/(mu0*mur)
J = 1e6
# Init state
Nx,Ny = 3,3
width,height = 0.6, 0.6
radiusHoles = 0.2*min(width,height)/max(Nx,Ny)
mesh = transformer(NCoils = 2,hmax=0.005)
TotalVolume = Integrate(CoefficientFunction(1)*dx(mesh.Materials("D")),mesh)
X,Y = np.meshgrid([-width/2 +(2*i+1)*(width)/(2*Nx) for i in range(Nx)],
                  [-height/2 + (2*i+1)*(height)/(2*Ny) for i in range(Ny)])
X,Y = X.flatten(),Y.flatten()
ls0 = GridFunction(H1(mesh))
with TaskManager():
    ls0.Set(sumRLS([hole([X[i],Y[i]],radiusHoles) for i in range(len(X))]))
Draw(mesh)
Draw(ls0,mesh,"Initial Design")
displayDict = {"D" : 1, "P1" : 0, "P2":0, "S1":0 ,"S2":0}
displayCF = CoefficientFunction([displayDict[mat] for mat in mesh.GetMaterials()])
Draw(IfPos(ls0,0,1)*displayCF,mesh,"Initial distribution")
FESolver = MagStat2D(mesh,dirichlet="out")
FESolver.set_ironBH(lambda b: nuiron)
FESolver.add_domain_tag([],['D'],['P1','P2','S1','S2'])
sourceDict = {"D" : 0, "P1" : J, "P2":-J, "S1":0 ,"S2":0}
jd = CoefficientFunction([sourceDict[mat] for mat in mesh.GetMaterials()])
FESolver.set_source(lambda v: jd*v)
FESolver.SolveState(ls0)
FESolver.SolveAdjoint(ls0,lambda *args : -1*dufObjectif(*args))
FESolver.show()
Draw(ShapeDeriv(FESolver,ls0),mesh,"ShapeDerivative")
def innerMinimization(phi, lam, beta, alpha, n_max = 50, alpha_min = 1e-4):
    from copy import copy
    n = 0
    FESolver.SolveState(phi)
    obj = fObjectif(FESolver.mesh,FESolver.gfu)
    objectiveHistory = [obj + lam * (VolFrac(FESolver.mesh,phi) - 0.1) + beta * (VolFrac(FESolver.mesh,phi) - 0.1) **2]
    singleObj = [obj]
    rhoHistory = [copy(phi)]
    massHistory = [VolFrac(FESolver.mesh,phi)]
    while ( n < n_max and alpha > alpha_min):
    
        # 1) Calcul de l'état physique :
        FESolver.SolveState(phi)
        
        # 2) Calcul de l'état adjoint :
        FESolver.SolveAdjoint(phi,lambda *args : -1*dufObjectif(*args))
        
        # 3) Calcul du gradient :
        gradient = -ShapeDeriv(FESolver,phi) + (lam  + 2* beta  * (VolFrac(FESolver.mesh,phi) - 0.1) )*ShapeDerivVolFrac(mesh,phi)

        # 4) Mise à jour :
        phi_temp = advectNGPy(phi,gradient ,step=alpha,dirichlet="",out=True)
        n += 1
        
        # 5) Contrôle du pas :
        FESolver.SolveState(phi_temp)
        singleObj.append(fObjectif(FESolver.mesh,FESolver.gfu))
        objectiveHistory.append(singleObj[-1]+ lam * (VolFrac(FESolver.mesh,phi_temp) - 0.1) + beta * (VolFrac(FESolver.mesh,phi_temp) - 0.1) **2)
        massHistory.append(VolFrac(FESolver.mesh,phi_temp))
        print(f'it n°{n} | L = {objectiveHistory[-1]} |f = {singleObj[-1]} | m = {massHistory[-1]} ', end = '\r')
        
        if objectiveHistory[-1] >= objectiveHistory[-2]:
            alpha = alpha/2
            objectiveHistory.pop(); massHistory.pop(); singleObj.pop()
        elif objectiveHistory[-1] < objectiveHistory[-2]:
            alpha = alpha*1.2
            if n%4 == 0:
                obj = redist(phi_temp)
                phi.vec.data = obj.Do()
            else:
                phi = phi_temp
            rhoHistory.append(copy(phi))
            Draw(IfPos(phi,0,1)*displayCF,mesh,"Distribution")
    Draw(gradient,mesh,"ShapeDeriv")
    return phi, alpha, massHistory, objectiveHistory,singleObj

# Optimization 
 ## Initialisation
lam = 0.2*1e-7
mass_target = 0.1 # Fraction massique cible (à faire varier)
beta = 1e-7          # Coefficient de pénalisation (à faire varier)
alpha = 1       # Pas initial
alpha_min = 1e-4  # Pas minimal
N_max = 15      # Nombre d'itérations de l'algorithme d'optimisation
volFracTarget = 0.1
phi = GridFunction(ls0.space)
phi.Set(ls0)
input('Start [press enter to continue]')
Masslist = []
LagrList = []
objList = []
lamList = [lam]
N = 0
while ( N < N_max):
    phi, alpha, lagrList, massList,singleObj = innerMinimization(phi, lam, beta, alpha = 1, n_max = 50, alpha_min = 1e-4)
    Draw(IfPos(phi,0,1)*displayCF,mesh,"Distribution")
    lam += 1000* 2* beta*(VolFrac(FESolver.mesh,phi) - mass_target)
    if (VolFrac(FESolver.mesh,phi) - mass_target)<0 and lam > 0 :
        lam += 9000* 2* beta*(VolFrac(FESolver.mesh,phi) - mass_target)
    lamList.append(lam)
    LagrList += lagrList
    Masslist += massList
    objList += singleObj
    print('Outer step',lam,N)
    N += 1
#vtk = VTKOutput(mesh,[phi,IfPos(phi,0,1)*displayCF,
#                          FESolver.gfu,FESolver.gfp,CoefficientFunction((grad(FESolver.gfu)[1],-grad(FESolver.gfu)[0])),
#                          -ShapeDerivative + (1 +volFracTarget -constraintList[-1])**2 *(constraintList[-1]>volFracTarget)*CShapeDerivative]
#                          ,
#                         ['phi','Distribution','State','Adjoint','Bfield',"ShapeDerivative"],f"C:/Users/Thomas/Desktop/LevelSet_SGE/data/optimIt_{i+1:03d}",subdivision=2)
#vtk.Do()
print(LagrList)
print(objList)
print(Masslist)
print(lamList)
input('stop')