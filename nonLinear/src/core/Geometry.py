import ngsolve
from ngsolve import Mesh
from netgen.geom2d import SplineGeometry, CSG2d, Rectangle

def transformer(NCoils :int = 1, 
                width : float = 0.6, 
                height : float= 0.6, 
                innerDiameterCoil : float = 4e-2, 
                outerDiameterCoil : float = 6e-2, 
                heightCoil : float = 4e-2, 
                hmax : float = 1e-2) -> ngsolve.Mesh :
    """
    Generate a finite element mesh for a 2D transformer cross-section using Netgen/NGSolve.

    Parameters
    ----------
    NCoils : int, optional
        Number of coils to include in the geometry.
        - If 1, the secondary coil is included in the optimization domain (multi-material topology optimization).
        - If 2, the secondary coil is fixed (air/iron topology optimization).
    width : float, optional
        Total width of the design domain (default is 0.6).
    height : float, optional
        Total height of the design domain (default is 0.6).
    innerDiameterCoil : float, optional
        Inner diameter of each coil (default is 0.04).
    outerDiameterCoil : float, optional
        Outer diameter of each coil (default is 0.06).
    heightCoil : float, optional
        Height of the coils (default is 0.04).
    hmax : float, optional
        Maximum mesh element size (default is 0.01).

    Returns
    -------
    mesh : ngsolve.Mesh
        A generated 2D finite element mesh including design and coil subdomains,
        suitable for electromagnetic topology optimization.

    Notes
    -----
    - The function uses `SplineGeometry` to define subdomains and materials.
    - For `NCoils == 1`, two coil regions (positive and negative primary windings) are included in the design domain.
    - For `NCoils == 2`, four fixed coils are added (positive/negative primary and secondary windings).
    - Materials are labeled as: "Omega_c" (design region), "Sp"/"Sm" (coils), and "Pp"/"Pm" for additional coil phases.

    """
    geo = SplineGeometry()
    geo.AddRectangle( p1 = (-width/2, -height/2), p2 = (width/2, height/2), leftdomain = 1, rightdomain = 0, bc = "out")

    if NCoils == 1:    # in case of multi-material topology optimization (including the secondary coil in the optimization)
        geo.AddRectangle( p1 = (-outerDiameterCoil/2, -heightCoil/2), p2 = (-innerDiameterCoil/2, heightCoil/2),
        leftdomain = 2, rightdomain = 1, bc='coil'  )
        geo.AddRectangle( p1 = (innerDiameterCoil/2, -heightCoil/2), p2 = (outerDiameterCoil/2, heightCoil/2),
        leftdomain = 3, rightdomain = 1, bc='coil' )
        geo.SetMaterial( 2, "Sp")      # Positive primary winding
        geo.SetMaterial( 3, "Sm")      # Negative primary winding

    elif NCoils == 2:  # in case of iron/air topology optimization (fixed secondary coil)
        th = (outerDiameterCoil - innerDiameterCoil) / 2
        ep = (outerDiameterCoil + innerDiameterCoil) / 2
        geo.AddRectangle( p1 = (-width/2 + (0.3-0.11), -heightCoil/2), 
                         p2 = (-width/2 + (0.3-0.11) + th, heightCoil/2),
        leftdomain = 2, rightdomain = 1, bc='coil' )
        geo.AddRectangle( p1 = (-width/2 + (0.3-0.11) + ep, -heightCoil/2), 
                         p2 = (-width/2 + (0.3-0.11) + ep + th, heightCoil/2),
        leftdomain = 3, rightdomain = 1, bc='coil' )
        geo.AddRectangle( p1 = (width/2 - (0.3-0.11) - th - ep, -heightCoil/2), 
                         p2 = (width/2 - (0.3-0.11) - ep, heightCoil/2),
        leftdomain = 4, rightdomain = 1, bc='coil' )
        geo.AddRectangle( p1 = (width/2 - (0.3-0.11) - th, -heightCoil/2), 
                         p2 = (width/2 - (0.3-0.11) , heightCoil/2),
        leftdomain = 5, rightdomain = 1, bc='coil' )
        geo.SetMaterial( 2, "Pp" )       # Positive primary winding
        geo.SetMaterial( 3, "Pm" )       # Negative primary winding
        geo.SetMaterial( 4, "Sm" )       # Positive secondary winding
        geo.SetMaterial( 5, "Sp" )       # Negative secondary winding

    geo.SetMaterial( 1, "Omega_c" )       # Design zone
    mesh = Mesh(geo.GenerateMesh( maxh = hmax )) # Generate the mesh with NetGen
    return mesh


def transformer_bit_array(width : float = 0.6, 
                          height : float = 0.6, 
                          innerDiameterCoil : float = 4e-2, 
                          outerDiameterCoil : float = 6e-2, 
                          heightCoil : float = 4e-2, 
                          hmax : float = 1e-2) -> ngsolve.Mesh :
    """
    Create a 2D parametric mesh for a transformer design using a bit array layout.

    This function defines a rectangular design region containing multiple small 
    subregions ("pixels" or "bits") for binary topology optimization, along with 
    fixed coil regions. The geometry is defined using Constructive Solid Geometry (CSG).

    Parameters
    ----------
    width : float, optional
        Total width of the design domain (default is 0.6).
    height : float, optional
        Total height of the design domain (default is 0.6).
    innerDiameterCoil : float, optional
        Inner diameter of the coils (default is 0.04).
    outerDiameterCoil : float, optional
        Outer diameter of the coils (default is 0.06).
    heightCoil : float, optional
        Height of the coils (default is 0.04).
    hmax : float, optional
        Maximum mesh element size (default is 0.01).

    Returns
    -------
    mesh : ngsolve.Mesh
        A generated 2D finite element mesh with coils and design subdomains labeled 
        for use in binary topology optimization problems.

    Notes
    -----
    - The mesh includes four coil rectangles labeled "P1", "P2", "S1", and "S2".
    - The design region is split into a grid of small rectangles (bit array), 
      excluding the coil areas.
    - Each subrectangle is assigned a unique material name ("D_ij") for individual 
      control in optimization.

    """
    
    geo = CSG2d()

    # Outer domain rectangle (full design area)
    recDomain = Rectangle(pmin=(-width/2, -height/2), pmax=(width/2, height/2), bc="out")

    # Coil parameters: thickness and offset
    th = (outerDiameterCoil - innerDiameterCoil) / 2
    ep = (outerDiameterCoil + innerDiameterCoil) / 2

    # Define 4 rectangular coils and assign materials
    recCoils1 = Rectangle(pmin=(-width/2 + (0.3-0.11), -heightCoil/2),
                          pmax=(-width/2 + (0.3-0.11) + th, heightCoil/2),
                          mat="P1", bc="coil")
    recCoils2 = Rectangle(pmin=(-width/2 + (0.3-0.11) + ep, -heightCoil/2),
                          pmax=(-width/2 + (0.3-0.11) + ep + th, heightCoil/2),
                          mat="P2", bc="coil")
    recCoilp1 = Rectangle(pmin=(width/2 - (0.3-0.11) - th - ep, -heightCoil/2),
                          pmax=(width/2 - (0.3-0.11) - ep, heightCoil/2),
                          mat="S1", bc="coil")
    recCoilp2 = Rectangle(pmin=(width/2 - (0.3-0.11) - th, -heightCoil/2),
                          pmax=(width/2 - (0.3-0.11), heightCoil/2),
                          mat="S2", bc="coil")

    # Add the coil shapes to the geometry
    geo.Add(recCoils1)
    geo.Add(recCoils2)
    geo.Add(recCoilp1)
    geo.Add(recCoilp2)

    # Create a central "remainDom" area (non-pixelized), labeled as material "D"
    Nsize = 20 # 30
    remainDom = (recDomain - recCoils1 - recCoils2 - recCoilp1 - recCoilp2 - Rectangle(pmin=(-7*width/Nsize, -5*height/Nsize),
                                       pmax=(7*width/Nsize, 5*height/Nsize))).Mat("D")

    # Subdivide the central area into a grid of bit rectangles
    for i in range(-7, 7): # TODO : generalize to arbitrary bit array size
        for j in range(-5, 5):
            x0, y0 = i * width / Nsize, j * height / Nsize
            x1, y1 = x0 + width / Nsize, y0 + height / Nsize

            # Create small subrectangle and subtract coil areas from it
            bit = Rectangle(pmin=(x0, y0), pmax=(x1, y1))
            bit = bit - recCoils1 - recCoils2 - recCoilp1 - recCoilp2
            bit = bit.Mat(f"D_{i+7}_{j+5}").Maxh(hmax)
            geo.Add(bit)

    geo.Add(remainDom) # Add remaining background domain
    mesh = Mesh(geo.GenerateMesh(maxh=hmax))  # Generate the mesh with NetGen
    return mesh


if __name__ == "__main__" : 
    print('testing core/Geometry.py ...')
    import numpy as np
    import netgen.gui 
    from netgen.webgui import Draw
    width = np.random.rand()+0.1
    height = np.random.rand()+0.1
    innerDiameterCoil = width*0.1
    outerDiameterCoil = width*0.13
    heightCoil = height * 0.1
    mesh1 = transformer(NCoils = 1, width = width, height = height, innerDiameterCoil = innerDiameterCoil, outerDiameterCoil = outerDiameterCoil, heightCoil = heightCoil, hmax = 1e-2)
    Draw(mesh1,"Transformer with 1 coil")
    input("Press Enter to inspect mesh2...")
    mesh2 = transformer(NCoils = 2, width = width, height = height, innerDiameterCoil = innerDiameterCoil, outerDiameterCoil = outerDiameterCoil, heightCoil = heightCoil, hmax = 1e-2)
    Draw(mesh2,"Transformer with 2 coil")
    input("Press Enter to inspect mesh3...")
    mesh3 = transformer_bit_array(width = width, height = height, innerDiameterCoil = innerDiameterCoil, outerDiameterCoil = outerDiameterCoil, heightCoil = heightCoil, hmax = 1e-2)
    Draw(mesh3,"Transformer bit-array")
    input("Press Enter to close the GUI and terminate the script...")
    print('done !')