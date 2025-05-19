# A basic lib that adds elementary tools for LS
from ngsolve.fem import IfPos,sqrt,specialcf,CoefficientFunction
from ngsolve.bla import Norm
from ngsolve.utils import x,y,grad,PyDet
from ngsolve.comp import GridFunction,L2,VectorL2,H1,VectorH1,BilinearForm,SymbolicBFI,VOL,BND
from ngsolve import TaskManager
import numpy as np


# Create initial levelset
def hole(c,radius):
    r = (x-c[0])**2 + (y-c[1])**2
    return IfPos(r-radius**2,sqrt(r-radius**2),-sqrt(radius**2-r))
def ellipse(c,radius,eps):
    r = (eps*x-c[0])**2 + (y-c[1])**2
    return IfPos(r-radius**2,sqrt(r-radius**2),-sqrt(radius**2-r))
def sumLS(ls1,ls2):
    return IfPos(ls1-ls2,ls2,ls1)
def sumRLS(list_ls):
    if len(list_ls) == 1 : return list_ls[0]
    curr_ls = list_ls[0]
    for i in range(1,len(list_ls)) : curr_ls = sumLS(curr_ls,list_ls[i])
    return curr_ls
#Returns volume fraction of negative part of triangle : from NGSOLVE tutorials
def GetVolumeFraction(psi0, psi1, psi2, EPS):
    if psi0 < -EPS and psi1 < -EPS and psi2 < -EPS:
        return 1
    elif psi0 > EPS and psi1 > EPS and psi2 > EPS:
        return 0
    
    ### if p0 is cut...
    if abs(psi0)<EPS and abs(psi1)>=EPS and abs(psi2)>=EPS: #psi0 is on vertex, psi1 and psi2 not
        if np.sign(psi1)!=np.sign(psi2):        # edge p1p2 is cut
            s = 1 * (-psi1)/(psi2 - psi1)
            if psi1 > 0:
                s = 1-s
        else:   #element is only cut in p0 (tangent)
            s = -0.5 * np.sign(psi1) + 0.5      # if psi1 <0 (-->psi2<0) then s=1, else (i.e. psi1, psi2>0) then s=0)
        return s
    elif abs(psi0)<EPS and abs(psi1)<EPS and abs(psi2)>=EPS: #edge p0p1 is on interface
        if psi2>0:
            s=0
        else:
            s=1
        return s
    elif abs(psi0)<EPS and abs(psi1)>=EPS and abs(psi2)<EPS: #edge p0p2 is on interface
        if psi1>0:
            s=0
        else:
            s=1
        return s
            
    ### if p1 is cut...
    if abs(psi1)<EPS and abs(psi2)>=EPS and abs(psi0)>=EPS: #psi1 is on vertex, psi2 and psi0 not
        if np.sign(psi2)!=np.sign(psi0):        # edge p2p0 is cut
            s = 1 * (-psi2)/(psi0 - psi2)
            if psi2 > 0:
                s = 1-s
        else:   #element is only cut in p1 (tangent)
            s = -0.5 * np.sign(psi2) + 0.5      # if psi2 <0 (-->psi0<0) then s=1, else (i.e. psi2, psi0>0) then s=0)
        return s
    elif abs(psi1)<EPS and abs(psi2)<EPS and abs(psi0)>=EPS: #edge p1p2 is on interface
        if psi0>0:
            s=0
        else:
            s=1
        return s
    elif abs(psi1)<EPS and abs(psi2)>=EPS and abs(psi0)<EPS: #edge p1p0 is on interface
        if psi2>0:
            s=0
        else:
            s=1
        return s
            
    ### if p2 is cut...
    if abs(psi2)<EPS and abs(psi0)>=EPS and abs(psi1)>=EPS: #psi2 is on vertex, psi0 and psi1 not
        if np.sign(psi0)!=np.sign(psi1):        # edge p0p1 is cut
            s = 1 * (-psi0)/(psi1 - psi0)
            if psi0 > 0:
                s = 1-s
        else:   #element is only cut in p1 (tangent)
            s = -0.5 * np.sign(psi0) + 0.5      # if psi0 <0 (-->psi1<0) then s=1, else (i.e. psi0, psi1>0) then s=0)
        return s
    elif abs(psi2)<EPS and abs(psi0)<EPS and abs(psi1)>=EPS: #edge p2p0 is on interface
        if psi1>0:
            s=0
        else:
            s=1
        return s
    elif abs(psi2)<EPS and abs(psi0)>=EPS and abs(psi1)<EPS: #edge p1p0 is on interface
        if psi0>0:
            s=0
        else:
            s=1
        return s
         
    if abs(psi0)>=EPS and abs(psi1)>=EPS and abs(psi2)>=EPS:    
        if np.sign(psi1) == np.sign(psi2) and np.sign(psi0) != np.sign(psi1):
            s = psi0*psi0 / ((psi1 - psi0)*(psi2 - psi0))
            if psi0 > 0:
                s = 1-s
        elif np.sign(psi2) == np.sign(psi0) and np.sign(psi1) != np.sign(psi2):
            s = psi1*psi1 / ( (psi2-psi1)*(psi0-psi1) )
            if psi1 > 0:
                s = 1-s
        elif np.sign(psi0) == np.sign(psi1) and np.sign(psi2) != np.sign(psi0):
            s = psi2*psi2 / ( (psi0-psi2)*(psi1-psi2) )
            if psi2 > 0:
                s = 1-s
        else:
            print(psi0, psi1, psi2)
            print("Error. This should not happen.")
    else:       #all three values are below EPS
        print(psi0, psi1, psi2)
        print("Error. This shoul not happen. Think about re-initializing your level set function as a signed distance function.")
    return s
def InterpolateLevelSetToElems(phi, val1, val2, mesh, EPS=1e-9):
    s  = 0
    nCutElems = 0
    phi_proj = GridFunction(L2(mesh))
    for i,el in enumerate(mesh.ngmesh.Elements2D()):

        psi0 = phi.vec[el.vertices[0].nr-1]
        psi1 = phi.vec[el.vertices[1].nr-1]
        psi2 = phi.vec[el.vertices[2].nr-1]
        
        if psi0 < 0 and psi1 < 0 and psi2 < 0:
            phi_proj.vec[i] = val1
        elif psi0>0 and psi1>0 and psi2 > 0:
            phi_proj.vec[i] = val2
        else:
            s = GetVolumeFraction(psi0, psi1, psi2, EPS)  #s... volume fraction of negative part of triangle
            nCutElems += 1
            phi_proj.vec[i] = val2 + s*(val1-val2) #== s*val1 + (1-s)*val2
    return phi_proj

def __ref_to_mesh(tri_coords,ref_coords):
    x_coord = tri_coords[0][0] + ref_coords[0] * (tri_coords[1][0]-tri_coords[0][0]) + ref_coords[1] * (tri_coords[2][0]-tri_coords[0][0])
    y_coord = tri_coords[0][1] + ref_coords[0] * (tri_coords[1][1]-tri_coords[0][1]) + ref_coords[1] * (tri_coords[2][1]-tri_coords[0][1])
    return x_coord,y_coord

def GetIso(phi,data_to_int,mask):
    mesh = phi.space.mesh
    pnts = mesh.ngmesh.Points()
    grad_phi = GridFunction(VectorL2(mesh)) ; grad_phi.Set(grad(phi))
    jacob = GridFunction(L2(mesh)) ; jacob.Set(PyDet(specialcf.JacobianMatrix(mesh.dim)))
    source_vec = GridFunction(VectorH1(mesh))
    for i,el in enumerate(mesh.ngmesh.Elements2D()):
        if mask.vec[i] == 0: continue
        psi0 = phi.vec[el.vertices[0].nr-1]
        psi1 = phi.vec[el.vertices[1].nr-1]
        psi2 = phi.vec[el.vertices[2].nr-1]
        if psi0>0 and psi1>0 and psi2>0 : continue#phi_proj.vec[i] = 0 #Check if iso-0 intersects the triangle
        elif psi0<0 and psi1<0 and psi2<0 : continue#phi_proj.vec[i] = 0
        else: 
            if np.sign(psi0) == np.sign(psi1): #cut seg 2 and seg 3
                pts = [(psi2/(psi2-psi1),1-psi2/(psi2-psi1)),(0,psi0/(psi0-psi2))]
            if np.sign(psi1) == np.sign(psi2): #cut seg 3 and seg 1
                pts = [(0,psi0/(psi0-psi2)),(psi0/(psi0-psi1),0)]
            if np.sign(psi0) == np.sign(psi2): #cut seg 1 and seg 2
                pts = [(psi0/(psi0-psi1),0),(psi2/(psi2-psi1),1-psi2/(psi2-psi1))]
            tri_coords = [pnts[el.vertices[i].nr] for i in range(3)]
            x1,y1,x2,y2 = *__ref_to_mesh(tri_coords,pts[0]), *__ref_to_mesh(tri_coords,pts[1]) #get intersec in absolute mesh coordinate from ref element cut
            if grad_phi.components[0].vec[i]*(y2-y1) + grad_phi.components[1].vec[i]*(x1-x2) > 0 : #check if the normal from the segment is outward facing
                x1,x2,y1,y2 = x2,x1,y2,y1
            Na,Nb,Nc = line_to_integralP0y(np.sign((x2-x1))*(y2-y1)/(x2-x1),y1-(y2-y1)*x1/(x2-x1), y1,y2,x1,x2) #mass term for each node P0
            for m in range(2):
                source_vec.components[m].vec[el.vertices[0].nr-1] += data_to_int.components[m].vec[i] * Na#*abs(jacob.vec[i])
                source_vec.components[m].vec[el.vertices[1].nr-1] += data_to_int.components[m].vec[i] * Nb#*abs(jacob.vec[i])
                source_vec.components[m].vec[el.vertices[2].nr-1] += data_to_int.components[m].vec[i] * Nc#*abs(jacob.vec[i])
            norm = sqrt((y2-y1)**2 + (x1-x2)**2) #this is the length of the cut
    return source_vec


def GetIsoTest(phi):
    mesh = phi.space.mesh
    pnts = mesh.ngmesh.Points()
    phi_proj = GridFunction(L2(mesh))
    grad_phi = GridFunction(VectorL2(mesh)) ; grad_phi.Set(grad(phi))
    jacob = GridFunction(L2(mesh)) ; jacob.Set(PyDet(specialcf.JacobianMatrix(mesh.dim)))
    source = GridFunction(H1(mesh))
    source_vec = GridFunction(VectorH1(mesh))
    list_coords,list_normal,list_vec,list_tri_norm = [],[],[],[]
    for i,el in enumerate(mesh.ngmesh.Elements2D()):
        #if mask.vec[i] == 0: continue
        psi0 = phi.vec[el.vertices[0].nr-1]
        psi1 = phi.vec[el.vertices[1].nr-1]
        psi2 = phi.vec[el.vertices[2].nr-1]
        if psi0>0 and psi1>0 and psi2>0 : continue#phi_proj.vec[i] = 0 #Check if iso-0 intersects the triangle
        elif psi0<0 and psi1<0 and psi2<0 : continue#phi_proj.vec[i] = 0
        else: 
            if np.sign(psi0) == np.sign(psi1): #cut seg 2 and seg 3
                pts = [(psi2/(psi2-psi1),1-psi2/(psi2-psi1)),(0,psi0/(psi0-psi2))]
            if np.sign(psi1) == np.sign(psi2): #cut seg 3 and seg 1
                pts = [(0,psi0/(psi0-psi2)),(psi0/(psi0-psi1),0)]
            if np.sign(psi0) == np.sign(psi2): #cut seg 1 and seg 2
                pts = [(psi0/(psi0-psi1),0),(psi2/(psi2-psi1),1-psi2/(psi2-psi1))]
            tri_coords = [pnts[el.vertices[i].nr] for i in range(3)]
            list_tri_norm.append((((tri_coords[0][0]+tri_coords[1][0]+tri_coords[2][0])/3),((tri_coords[0][1]+tri_coords[1][1]+tri_coords[2][1])/3),tri_coords[2][1]-tri_coords[1][1],tri_coords[2][0]-tri_coords[1][0]))
            x1,y1,x2,y2 = *__ref_to_mesh(tri_coords,pts[0]), *__ref_to_mesh(tri_coords,pts[1]) #get intersec in absolute mesh coordinate from ref element cut
            if grad_phi.components[0].vec[i]*(y2-y1) + grad_phi.components[1].vec[i]*(x1-x2) > 0 : #check if the normal from the segment is outward facing
                x1,x2,y1,y2 = x2,x1,y2,y1
            phi_proj.vec[i] = 1
            Na,Nb,Nc = line_to_integralP0y(np.sign((x2-x1))*(y2-y1)/(x2-x1),y1-(y2-y1)*x1/(x2-x1), y1,y2,x1,x2) #mass term for each node P0
            #Na,Nb,Nc = line_to_integralP0x(np.sign(x1-x2)*(y2-y1)/(x2-x1),y1-(y2-y1)*x1/(x2-x1), x1,x2) #mass term for each node P0
            source.vec[el.vertices[0].nr-1] += Na*abs(jacob.vec[i]) #multiply by jacobian to scale properly
            source.vec[el.vertices[1].nr-1] += Nb*abs(jacob.vec[i])
            source.vec[el.vertices[2].nr-1] += Nc*abs(jacob.vec[i])
            for m in range(2):
                source_vec.components[m].vec[el.vertices[0].nr-1] += grad_phi.components[m].vec[i] * Na*abs(jacob.vec[i])
                source_vec.components[m].vec[el.vertices[1].nr-1] += grad_phi.components[m].vec[i] * Nb*abs(jacob.vec[i])
                source_vec.components[m].vec[el.vertices[2].nr-1] += grad_phi.components[m].vec[i] * Nc*abs(jacob.vec[i])
            norm = sqrt((y2-y1)**2 + (x1-x2)**2) #this is the length of the cut
            normal = ((y2-y1)/norm,(x1-x2)/norm)
            list_coords.append(((x1+x2)/2,(y1+y2)/2))
            list_vec.append((x1,y1,x2,y2))
            list_normal.append(normal)
    return phi_proj,list_coords,list_normal,source,list_vec,source_vec,list_tri_norm

def line_to_integralP0x(a,b,x_i,x_j):
    det = np.sqrt(1+a**2)
    d_x1,d_x2 = (x_j-x_i),(x_j**2-x_i**2)
    Na = (1-b)*d_x1 - (1+a)*d_x2/2
    Nb = d_x2/2
    Nc = a*d_x2/2 + b*d_x1
    return Na*det,Nb*det,Nc*det


def line_to_integralP0y(a,b,y_i,y_j,x_i,x_j):
    if a == 0 :
        d_x1,d_x2 = x_j-x_i,x_j**2-x_i**2
        Na = (1-b)*d_x1 - d_x2/2
        Nb = d_x2/2
        Nc = b*d_x1
        det = 1
    else :  
        if a<0 : y_i,y_j = y_j,y_i
        det = np.sqrt(1 + 1/a**2)
        d_y1,d_y2 = y_j-y_i,y_j**2-y_i**2
        Na = (1+b/a)*d_y1 - (1+1/a)*d_y2/2
        Nb = 1/a * d_y2/2 - b/a * d_y1
        Nc = d_y2/2
    return Na*det,Nb*det,Nc*det 

def line_to_integralP1(a,b,x_i,x_j):
    dx_1,dx_2,dx_3 = (x_j-x_i), (x_j**2-x_i**2), (x_j**3-x_i**3)
    NaNa = (1-b)*dx_1+ (1-b)(1+a)*d_x2 + (1+a)*dx_3/3    
    NaNb = (1-b)*dx_2/2 - (1+a)*d_x3/3
    NaNc = (1-b)*b*dx_1 - (a*b+1)*dx_2/2 - (1+a)*a*d_x3/2
    NbNb = dx_3/3
    NbNc = a*dx_3/3 + b*dx_2/2
    NcNc = a**2*dx_3/3 + a*b*dx_2 + b**2*dx_1
    return NaNa,NbNb,NcNc,NaNb,NaNc,NbNc

def advectNGPy(phi,vel,step=1,dirichlet='dirichlet',out=False):
    mesh = phi.space.mesh
    fes = L2(mesh,order=1,dirichlet=dirichlet)
    phi_temp = GridFunction(fes) ; phi_temp.Set(phi)
    u,v = fes.TrialFunction(), fes.TestFunction()
    n = specialcf.normal(mesh.dim)
    gf_ms = GridFunction(L2(mesh))
    vel = vel
    gf_ms.Set(specialcf.mesh_size/(1e-12+Norm(vel)))
    dt = min(np.array(gf_ms.vec))
    c = BilinearForm(fes)
    c += SymbolicBFI( vel * grad(u) * v)
    bn = vel*n
    vin = IfPos(bn,v.Other(v),v)
    c += SymbolicBFI( bn*(u.Other(u) - u) * vin, VOL, skeleton=True)
    c += SymbolicBFI( IfPos(bn, 0, bn) * (u.Other(u) - u) * v, BND, skeleton=True)
    res = phi_temp.vec.CreateVector()
    with TaskManager():
        c.Apply(phi_temp.vec,res)
        fes.SolveM(res,rho=CoefficientFunction(1.0))
    phi_temp.vec.data += dt * step * res
    if out:
        phi_out = GridFunction(phi.space)
        phi_out.Set(phi_temp)
        return phi_out
    else:
        phi.Set(phi_temp)

def redistNGPy(phi):
    fes = L2(mesh,order=1)
    phi_temp = GridFunction(fes) ; phi_temp.Set(phi)
    u,v = fes.TrialFunction(), fes.TestFunction()
    n = specialcf.normal(mesh.dim)
    sign_phi = IfPos(phi_temp,1,-1)
    vel = sign_phi*grad(phi_temp)/Norm(grad(phi_temp))
    gf_ms = GridFunction(L2(mesh))
    gf_ms.Set(specialcf.mesh_size/Norm(vel))
    dt = min(np.array(gf_ms.vec))
    c = BilinearForm(fes)
    c += SymbolicBFI( vel * grad(u) * v)
    bn = vel*n
    vin = IfPos(bn,v.Other(v),v)
    c += SymbolicBFI( bn*(u.Other(u) - u) * vin, VOL, skeleton=True)
    c += SymbolicBFI( IfPos(bn, 0, bn) * (u.Other(u) - u) * v, BND, skeleton=True)
    res = phi_temp.vec.CreateVector()
    c.Apply(phi_temp.vec,res)
    fes.SolveM(res,rho=CoefficientFunction(1.0))
    phi_temp.vec.data += dt * res