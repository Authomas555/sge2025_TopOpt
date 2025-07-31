from ngsolve import atan
import matplotlib.pyplot as plt
import numpy as np
import inspect

interpolationLibrary = {}
# Feel free to add your own interpolation methods here, following the same structure as below. Keys should be lower case and unique.
# Conditions : 
# - f differentiable (C1 is better)
# - f(0) = 0
# - f(1) = 1
# - f'(x) >= 0 on [0,1]
# - give a default values for parameters (other than x)

# other interpolations that depend on the values to interpolate (not compatible) : 
# - geometric interpolation (Dyck & Lowther 1996, Labbé 2010)
# - sequence MIS (Sanogo & Messine, 2018)

##############################################################################################################################
interpolationLibrary["power"] = {"function" : lambda x, p = 1 : x ** p,
                                 "derivative" : lambda x, p = 1 : p * x ** (p-1),
                                 "condition" : (lambda p : p>0, "p should be strictly positive"),
                                 "name" : "power law",
                                 "latex" : "x^p",
                                 "reference" : "Bendsøe (1989), Mlejnek (1992), Bendsøe and Sigmund (1999, 2003)",
                                 "alias" : ["power law", "pow", "simp"]}

##############################################################################################################################

interpolationLibrary["ramp"] =  {"function" : lambda  x, p = 0 : x / (1 + p*(1-x)),
                                 "derivative" : lambda x, p = 0 : (p + 1)/(1 + p*(1-x))**2,
                                 "name" : "rational function",
                                 "condition" : (lambda p : p>=0, "p should be positive"),
                                 "latex" : "\\frac{x}{1 + p(1-x)}" ,                       # Python doesn't like the backslash
                                 "reference" : "Stolpe and Svanberg (2001), Hansen (2005)",
                                 "alias" : ["rational", "rat"]}
# aliases
interpolationLibrary["rational"] = interpolationLibrary["ramp"]

##############################################################################################################################

interpolationLibrary["polynomial"] =  {"function" : lambda  x, p = 1, a = 3 : x/a + (a-1)/a*x**p,
                                 "derivative" : lambda x, p = 1, a = 3 : 1/a + (p*x**(p - 1)*(a - 1))/a,
                                 "name" : "rational function",
                                 "condition" : (lambda p : p>0, "p should be strictly positive"),
                                 "latex" : "\\frac{x}{a} + \\frac{a-1}{a} x^p" ,
                                 "reference" : "Jihong Zhu (2008, PhD Thesis)"}
# aliases
interpolationLibrary["zhu"] = interpolationLibrary["polynomial"]
interpolationLibrary["poly"] = interpolationLibrary["polynomial"]

##############################################################################################################################

interpolationLibrary["atan"] =  {"function" : lambda  x, p = 1e-6 : (1+atan(p*(2*x-1))/atan(p))/2,
                                 "derivative" : lambda x, p = 1e-6 : p / (atan(p) * (p**2*(2*x - 1)**2 + 1)),
                                 "name" : "rational function",
                                 "condition" : (lambda p : p>0, "p should be strictly positive"),
                                 "latex" : "\\frac{1+\\atan(p(2x-1))}{2\\atan(p)}",
                                 "reference" : "Lukáš (2006, An Integration of Optimal Topology and Shape Design for Magnetostatics)"}
# aliases
interpolationLibrary["lukas"] = interpolationLibrary["atan"]
interpolationLibrary["arctan"] = interpolationLibrary["atan"]

##############################################################################################################################


# Add aliases to the library
keys = list(interpolationLibrary.keys())
for key in keys:
    if "alias" in interpolationLibrary[key]:
        for alias in interpolationLibrary[key]["alias"]:
            interpolationLibrary[alias] = interpolationLibrary[key]

###############################################################################################################################

def get_default_args(func):
    """ Get default argument of a function """
    signature = inspect.signature(func)
    return {
        k: v.default
        for k, v in signature.parameters.items()
        if v.default is not inspect.Parameter.empty
    }

class Interpolation:
    """ Density interpolation class """
    name : str
    function : callable
    derivative : callable
    parameters : dict
    latex : str
    reference : str

    def __init__(self, type, **kwargs):
        type = type.lower()
        if type not in interpolationLibrary:
            raise ValueError(f"Interpolation type '{type}' not recognized. Available : \n{'\n'.join(interpolationLibrary.keys())}.")
        self.type = type
        # check condition on p
        if kwargs is not None:
            if "p" in kwargs : 
                condition, message = interpolationLibrary[type]["condition"]
                if not condition(kwargs["p"]):
                    raise ValueError(message)

        self.parameters = get_default_args(interpolationLibrary[type]["function"])
        for keys in kwargs.keys():
            self.parameters[keys] = kwargs[keys]

        self.function = lambda x : interpolationLibrary[type]["function"](x, **self.parameters).real
        self.derivative = lambda x :  interpolationLibrary[type]["derivative"](x, **self.parameters).real
        self.name = interpolationLibrary[type]["name"]
        self.latex = interpolationLibrary[type]["latex"]
        self.reference = interpolationLibrary[type]["reference"]
    
    def __call__(self, x):
        return self.function(x)

    def __str__(self):
        return self.latex + " | " + self.__str_parameters()
    
    def __str_parameters(self):
        return  ", ".join([f"{k}={v}" for k, v in self.parameters.items()])

    def plot(self, num_points = 100, label = "default", param = "default"):
        x = np.linspace(0, 1, num_points)
        if label == "default": label = self.name
        if param == "default": param = " | " + self.__str_parameters()
        plt.plot(x, self(x), label = label + param)

        

if __name__ == "__main__" : # simple tests
    I = Interpolation("zhu", p = 3)
    I.plot()
    print(I.latex)

