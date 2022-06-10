import sympy as sy
import numpy as ny
import matplotlib.pyplot as plt

def trapecio(f,a,b):
    x = sy.symbols('x')
    fx = sy.sympify(f)
    I =((sy.N(fx.subs(x,a))+sy.N(fx.subs(x,b)))*(b-a))/2
    d2f = sy.diff(fx,x,2)
    error = error_trapecio(d2f,a,b)
    return I

def simpson(f,a,b):
    x = sy.symbols('x')
    fx = sy.sympify(f)
    x1 = (a+b)/2
    I = ((b-a)/6)*(sy.N(fx.subs(x,a))+4*sy.N(fx.subs(x,x1))+sy.N(fx.subs(x,b)))
    d4f = sy.diff(fx,x,4)
    error = error_simpson(d4f,a,b)
    return I

def error_trapecio(f,a,b):
    x = sy.symbols('x')
    fx = sy.sympify(f)
    d1f = sy.diff(fx,x)
    p_c = sy.solve(d1f,x)
    posibles_maximos = []
    for i in p_c:
        if a < i < b:
            posibles_maximos.append(abs(sy.N(fx.subs(x,i))))
    posibles_maximos.append(abs(sy.N(fx.subs(x,a))))
    posibles_maximos.append(abs(sy.N(fx.subs(x, b))))
    vmax = (max(posibles_maximos))
    error = ((b-a)**3)/12 * vmax
    return error

def error_simpson(f,a,b):
    x = sy.symbols('x')
    fx = sy.sympify(f)
    d1f = sy.diff(fx,x)
    p_c = sy.solve(d1f,x)
    posibles_maximos = []
    for i in p_c:
        if a < i < b:
            posibles_maximos.append(abs(sy.N(fx.subs(x,i))))
    posibles_maximos.append(abs(sy.N(fx.subs(x,a))))
    posibles_maximos.append(abs(sy.N(fx.subs(x, b))))
    vmax = (max(posibles_maximos))
    error = ((b-a)**5)/2880 * vmax
    return error

def simpson_compuesto(f,a,b,m):
    I = 0
    h = (b-a)/(m-1)
    x0 = a
    for i in range(m-1):
        xi = x0 + i*h
        xii = x0 + (i+1)*h
        I = I + simpson(f,xi,xii)
    return I

def trapecio_compuesto(f,a,b,m):
    I = 0
    h = (b-a)/(m-1)
    x0 = a
    for i in range(m-1):
        xi = x0 + i*h
        xii = x0 + (i+1)*h
        I = I + trapecio(f,xi,xii)
    return I


def cuadraturas_gausseanas(f,n):
    x = sy.symbols('x')
    fx = sy.sympify(f)
    I = 0
    xi = [[0],[0.57735,-0.57735],[0.77459,-0.77459,0],[0.86113,-0.86113,0.33998,-0.33998],
         [0.90617,-0.90617,0.53846,-0.53846,0],[0.53846,-0.53846,0.6612,-0.6612,0.23861,-0.23861],
         [0.9491,-0.9491,0.74153,-0.74153,0.40584,-0.40584,0],
         [0.96028,-0.96028,0.79666,-0.79666,0.52553,-0.52553,0.18343,-0.18343],
         [0.96816,-0.96816,0.83603,-0.83603,0.61337,-0.61337,0.32425,-0.32425,0],
         [0.9739,-0.9739,0.86506,-0.86506,0.6794,-0.6794,0.43339,-0.43339,0.14887,-0.14887]]
    wi = [[2],[1, 1], [0.55556, 0.55556, 0.88889], [0.34785, 0.34785, 0.65214, 0.65214],
         [0.23692, 0.23692, 0.47862, 0.47862, 0.56889], [0.17132, 0.17132, 0.36076, 0.36076, 0.46791, 0.46791],
         [0.12948, 0.12948, 0.2797, 0.2797, 0.38183, 0.38183, 0.41795],
         [0.10122, 0.10122, 0.22238, 0.22238, 0.3137, 0.3137, 0.36268, 0.36268],
         [0.08127, 0.08127, 0.18064, 0.18064, 0.26061, 0.26061, 0.31234, 0.31234, 0.33023],
         [0.06667, 0.06667, 0.14945, 0.14945, 0.21908, 0.21908, 0.26926, 0.26926, 0.29552, 0.29552]]
    for i in range(len(xi[n-1])):
        I = I + (wi[n-1][i]*sy.N(fx.subs(x,xi[n-1][i])))
    return I

def c_g_general(f,n,a,b):
    x = sy.symbols('x')
    fx = sy.sympify(f)
    gx = (b-a)/2 * (fx.subs(x,((b-a)*x + b + a)/2))
    I = cuadraturas_gausseanas(gx,n)
    return I