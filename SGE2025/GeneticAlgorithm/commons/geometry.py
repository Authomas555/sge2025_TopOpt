from ngsolve import *
from netgen.geom2d import SplineGeometry
from netgen.geom2d import CSG2d, Circle, Rectangle

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



def transformerCG(NCoils = 1, width = 0.6, height = 0.6, innerDiameterCoil = 4e-2, outerDiameterCoil = 6e-2, heightCoil = 4e-2, hmax = 1e-2):
    geo = CSG2d()
    recDomain = Rectangle( pmin=(-width/2, -height/2), pmax=(width/2, height/2), bc="out")
    if NCoils == 1:
        recCoils1 = Rectangle( pmin=(-outerDiameterCoil/2, -heightCoil/2), pmax=(-innerDiameterCoil/2, heightCoil/2), mat="S1", bc="coil")
        recCoils2 = Rectangle( pmin=(innerDiameterCoil/2, -heightCoil/2), pmax=(outerDiameterCoil/2, heightCoil/2), mat="S2" , bc="coil")
        geo.Add(recCoils1)
        geo.Add(recCoils2)
        remainDom = (recDomain-recCoils1 - recCoils2).Mat("D")
        geo.Add(remainDom)
    elif NCoils == 2:
        th = (outerDiameterCoil - innerDiameterCoil) / 2
        ep = (outerDiameterCoil + innerDiameterCoil) / 2
        recCoils1 = Rectangle( pmin=(-width/2 + (0.3-0.11), -heightCoil/2), pmax=(-width/2 + (0.3-0.11) + th, heightCoil/2), mat="P1", bc="coil")
        recCoils2 = Rectangle( pmin=(-width/2 + (0.3-0.11) + ep, -heightCoil/2), pmax=(-width/2 + (0.3-0.11) + ep + th, heightCoil/2), mat="P2" , bc="coil")
        recCoilp1 = Rectangle( pmin=(width/2 - (0.3-0.11) - th - ep, -heightCoil/2), pmax=(width/2 - (0.3-0.11) - ep, heightCoil/2), mat="S1", bc="coil")
        recCoilp2 = Rectangle( pmin=(width/2 - (0.3-0.11) - th, -heightCoil/2), pmax=(width/2 - (0.3-0.11) , heightCoil/2), mat="S2" , bc="coil")
        geo.Add(recCoils1)
        geo.Add(recCoils2)
        geo.Add(recCoilp1)
        geo.Add(recCoilp2)
        remainDom = (recDomain-Rectangle(pmin=(-7*width/30,-5*height/30),pmax=(7*width/30,5*height/30))).Mat("D")
        for i in range(-7,7):
            for j in range(-5,5):
                x0,y0 =  i*width/30, j*height/30
                x1,y1 = x0 + width/30,y0 + height/30
                geo.Add((Rectangle(pmin=(x0,y0), pmax=(x1,y1)) - recCoils1  - recCoils2 - recCoilp1 - recCoilp2).Mat(f"D_{i+7}_{j+5}").Maxh(1e-2))
                
        geo.Add(remainDom)
    ngmesh = geo.GenerateMesh( maxh = hmax ) # Maillage avec NetGen
    mesh = Mesh( ngmesh )    # Importation du maillage NetGen dans NgSolve
    return mesh