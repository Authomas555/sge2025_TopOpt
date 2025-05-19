import warnings
import numpy as np
from ngsolve.solvers import Newton
from ngsolve.solve import Draw
from ngsolve.comp import H1,VectorH1,GridFunction,BilinearForm,SymbolicBFI,LinearForm,SymbolicLFI
from ngsolve.fem import sqrt,CoefficientFunction
from ngsolve.utils import grad,dx,Sym,Id,Trace,atan2,x,y,cos,sin,Grad,ds
from ngsolve.bla import InnerProduct
from pyngcore import TaskManager
from levelBaseLib import InterpolateLevelSetToElems


mu0 = 4*np.pi/(10**7) ; nu0 = 1/mu0

def jacobi2D(u_x,u_y):
    '''Implementation of Grad of 2D vector ~easier to use for separate components'''
    d_xx,d_xy = grad(u_x)[0],grad(u_x)[1]
    d_yx,d_yy = grad(u_y)[0],grad(u_y)[1]
    return CoefficientFunction((d_xx,d_xy,d_yx,d_yy),dims=(2,2))

def strainTostress(strain,E,nu):
    mu = E/(2*(1+nu)) ; lb = E*nu/((1+nu)*(1-2*nu))
    return 2*mu*strain + lb*Trace(strain)*Id(2)

def concatenate(list_tag):
    s = ''
    for e in list_tag : s += e + '|'
    return s[:-1]


class MagStat2D:
    
    def __init__(self,mesh,dirichlet=None):
        self.mesh = mesh
        self.fes = H1(mesh,dirichlet=dirichlet,order=1)
        self.gfu,self.gfp = GridFunction(self.fes),GridFunction(self.fes)

    def add_domain_tag(self,iron_tag=[], working_area_tag = [], air_tag=[]):
        self.tags = {'iron':concatenate(iron_tag),'air':concatenate(air_tag),'working_area':concatenate(working_area_tag)}
        missing_dom = list(set(self.mesh.GetMaterials())-set(iron_tag+air_tag+working_area_tag))
        if missing_dom != [] : 
            print('Warning: domains ' + str(missing_dom) + ' initialized as air')
            self.tags['air'] += '|'+concatenate(missing_dom)
    
    def set_ironBH(self,nuB):
        self.nu_iron = nuB

    def set_source(self,source_state):
        self.source_state = source_state 

    def SolveState(self,phi=None):
        # Get test and trial function from FE-space
        u,v = self.fes.TnT()
        # Compute field and material properties
        Bnorm = sqrt(1e-12+grad(u)*grad(u))
        nuB = self.nu_iron(Bnorm)
        # Initialize forms
        ## Bilinear 
        K_matrix = BilinearForm(self.fes)
        K_matrix += SymbolicBFI(nu0*grad(u)*grad(v),definedon=~self.mesh.Materials(self.tags['iron']+'|'+self.tags['working_area']))
        K_matrix += SymbolicBFI(nuB*grad(u)*grad(v),definedon=self.mesh.Materials(self.tags['iron']))
        if phi is not None : K_matrix += SymbolicBFI((nu0  + (nuB-nu0)*InterpolateLevelSetToElems(phi,1,0,self.mesh))*grad(u)*grad(v),definedon=self.mesh.Materials(self.tags['working_area']))
        ## Linear
        K_matrix += -self.source_state(v)*dx
        self.gfu.Set(0)
        # Solve system
        with TaskManager():  Newton(K_matrix, self.gfu, freedofs = self.fes.FreeDofs() ,inverse="sparsecholesky",printing=False)

    def SolveAdjoint(self,phi,source_function):
        # Get test and trial function from FE-space
        u,v = self.fes.TnT()
        # Compute field and material properties
        Bnorm = sqrt(1e-12+grad(u)*grad(u))
        nuB = self.nu_iron(Bnorm)
        # Initialize forms
        ## Bilinear 
        K_matrix = BilinearForm(self.fes)
        K_matrix += SymbolicBFI(nu0*grad(u)*grad(v),definedon=~self.mesh.Materials(self.tags['iron']+'|'+self.tags['working_area']))
        K_matrix += SymbolicBFI(nuB*grad(u)*grad(v),definedon=self.mesh.Materials(self.tags['iron']))
        K_matrix += SymbolicBFI((nu0  + (nuB-nu0)*InterpolateLevelSetToElems(phi,1,0,self.mesh))*grad(u)*grad(v),definedon=self.mesh.Materials(self.tags['working_area']))
        # Linear form
        duJ = LinearForm(self.fes) ; duJ += source_function(self.mesh,self.gfu,v,phi) 
        # Assemble and solve
        with TaskManager():
            K_matrix.AssembleLinearization(self.gfu.vec) ; duJ.Assemble()
            self.gfp.vec.data = -(K_matrix.mat.Inverse(self.fes.FreeDofs(), inverse="sparsecholesky").T * duJ.vec)

    def getEquation(self,phi):
        Bnorm = sqrt(1e-12+grad(self.gfu)*grad(self.gfu))
        nuB = self.nu_iron(Bnorm)
        eq = nu0*grad(self.gfu)*grad(self.gfp)*dx(definedon=~self.mesh.Materials(self.tags['iron']+'|'+self.tags['working_area']))
        eq += nuB*grad(self.gfu)*grad(self.gfp)*dx(definedon=self.mesh.Materials(self.tags['iron']))
        eq += (nu0  + (nuB-nu0)*InterpolateLevelSetToElems(phi,1,0,self.mesh))*grad(self.gfu)*grad(self.gfp)*dx(definedon=self.mesh.Materials(self.tags['working_area']))
        eq += -1*self.source_state(self.gfp)*dx
        return eq
                

    def show(self):
        B = CoefficientFunction((grad(self.gfu)[1],-grad(self.gfu)[0])) 
        Draw(self.gfu,self.mesh,'State') ; Draw(B,self.mesh,'Bfield')
        Draw(self.gfp,self.mesh,'Adjoint')