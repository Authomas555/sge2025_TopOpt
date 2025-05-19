from ngsolve import *

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