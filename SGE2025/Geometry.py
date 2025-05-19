from ngsolve import *
from netgen.geom2d import SplineGeometry

def transformer(NCoils = 1, width = 0.6, height = 0.6, innerDiameterCoil = 4e-2, outerDiameterCoil = 6e-2, heightCoil = 4e-2, hmax = 1e-2):
    ## Geometry generation

    geo = SplineGeometry()
    geo.AddRectangle( p1 = (-width/2, -height/2), p2 = (width/2, height/2), leftdomain = 1, rightdomain = 0, bc = "out")

    if NCoils == 1:
        geo.AddRectangle( p1 = (-outerDiameterCoil/2, -heightCoil/2), p2 = (-innerDiameterCoil/2, heightCoil/2),
        leftdomain = 2, rightdomain = 1 )
        geo.AddRectangle( p1 = (innerDiameterCoil/2, -heightCoil/2), p2 = (outerDiameterCoil/2, heightCoil/2),
        leftdomain = 3, rightdomain = 1 )
        geo.SetMaterial( 2, "Sp",bc='coil' )      # conducteur primaire positif
        geo.SetMaterial( 3, "Sm",bc='coil' ) 
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
        geo.SetMaterial( 2, "Pp" )      # conducteur primaire positif
        geo.SetMaterial( 3, "Pm" )      # conducteur primaire négatif
        geo.SetMaterial( 4, "Sm" )      # conducteur secondaire négatif
        geo.SetMaterial( 5, "Sp" )      # conducteur secondaire positif

    # On nomme les différentes zones afin de les repérer plus tard :
    geo.SetMaterial( 1, "Omega_c" )       # zone de design
    ngmesh = geo.GenerateMesh( maxh = hmax ) # Maillage avec NetGen
    mesh = Mesh( ngmesh )    # Importation du maillage NetGen dans NgSolve
    return mesh