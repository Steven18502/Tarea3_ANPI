import sympy as sy
import numpy as ny
import matplotlib.pyplot as plt
from sympy.calculus.util import continuous_domain


def trapecio(f, a, b):
    """
    Ejecuta la aproximacion de la integral de la funcion f en el intervalo [a, b] por el
    metodo del trapecio

    :param str f: funcion en formato de string
    :param int/float a: limite inferior integral
    :param int/float b: limite superior integral
    """
    x = sy.symbols('x')
    fx = sy.sympify(f)

    signo = 1
    if (a > b):  # inversion de signo en caso de a > b
        a_temp = a
        a = b
        b = a_temp
        signo = -1

    if es_continua(fx, x, a, b):  # evaluar continuidad de funcion
        if (a != b):
            I = ((sy.N(fx.subs(x, a))+sy.N(fx.subs(x, b)))*(b-a))/2
            d2f = sy.diff(fx, x, 2)
            error = error_trapecio(d2f, a, b)
            return (I*signo, error)
        else:
            return (0, 0)  # caso a == b
    else:
        return "La funcion no es continua en el intervalo [" + str(a) + "," + str(b) + "]"


def error_trapecio(f, a, b):
    """
    Calcula cota de error para el metodo del trapecio

    :param str f: funcion en formato de string
    :param int/float a: limite inferior integral
    :param int/float b: limite superior integral
    """
    valor_max = vmax(f, a, b)
    error = ((b-a)**3)/12 * valor_max
    return error


def vmax(f, a, b):
    """
    Calcula el valor maximo de f(x) en un rango [a,b]
    """
    x = sy.symbols('x')
    fx = sy.sympify(f)
    # primera derivada para hallar puntos criticos de la funcion
    d1f = sy.diff(fx, x)
    p_c = sy.solve(d1f, x)
    posibles_maximos = []
    for i in p_c:
        if a < i < b:
            posibles_maximos.append(abs(sy.N(fx.subs(x, i))))
    posibles_maximos.append(abs(sy.N(fx.subs(x, a))))
    posibles_maximos.append(abs(sy.N(fx.subs(x, b))))
    return max(posibles_maximos)


def simpson(f, a, b):
    """
    Ejecuta la aproximacion de la integral de la funcion f en el intervalo [a, b] por el
    metodo de Simpson

    :param str f: funcion en formato de string
    :param int/float a: limite inferior integral
    :param int/float b: limite superior integral
    """
    x = sy.symbols('x')
    fx = sy.sympify(f)

    signo = 1
    if (a > b):  # inversion de signo en caso de a > b
        a_temp = a
        a = b
        b = a_temp
        signo = -1

    if es_continua(fx, x, a, b):  # evaluar continuidad de funcion
        if (a != b):
            x1 = (a+b)/2
            I = ((b-a)/6)*(sy.N(fx.subs(x, a))+4 *
                           sy.N(fx.subs(x, x1))+sy.N(fx.subs(x, b)))
            d4f = sy.diff(fx, x, 4)
            error = error_simpson(d4f, a, b)
            return (I*signo, error)
        else:
            return (0, 0)  # caso a == b
    else:
        return "La funcion no es continua en el intervalo [" + str(a) + "," + str(b) + "]"


def error_simpson(f, a, b):
    """
    Calcula cota de error para el metodo de Simpson

    :param str f: funcion en formato de string
    :param int/float a: limite inferior integral
    :param int/float b: limite superior integral
    """
    valor_max = vmax(f, a, b)
    error = ((b-a)**5)/2880 * valor_max
    return error


def regla_boole(f, a, b):
    """
    Ejecuta la aproximacion de la integral de la funcion f en el intervalo [a, b] por el
    metodo de la Regla de Boole

    :param str f: funcion en formato de string
    :param int/float a: limite inferior integral
    :param int/float b: limite superior integral
    """
    x = sy.symbols('x')
    fx = sy.sympify(f)

    signo = 1
    if (a > b):  # inversion de signo en caso de a > b
        a_temp = a
        a = b
        b = a_temp
        signo = -1

    if es_continua(fx, x, a, b):  # evaluar continuidad de funcion
        if (a != b):
            h = (b-a)/4
            f = []
            for i in range(5):
                f.append(sy.N(fx.subs(x, a + i*h)))
            I = 2*h/45*(7*f[0] + 32*f[1] + 12*f[2] + 32*f[3] + 7*f[4])
            d6f = sy.diff(fx, x, 6)
            error = error_boole(d6f, a, b, h)
            return (I*signo, error)
        else:
            return (0, 0)  # caso a == b
    else:
        return "La funcion no es continua en el intervalo [" + str(a) + "," + str(b) + "]"


def error_boole(f, a, b, h):
    """
    Calcula cota de error para el metodo de Boole

    :param str f: funcion en formato de string
    :param int/float a: limite inferior integral
    :param int/float b: limite superior integral
    :param float h
    """
    valor_max = vmax(f, a, b)
    error = (8*(h)**7)/945 * valor_max
    return error


def trapecio_compuesto(f, a, b, m):
    """
    Ejecuta la aproximacion de la integral de la funcion f en el intervalo [a, b] por el
    metodo de Simpson compuesto

    :param str f: funcion en formato de string
    :param int/float a: limite inferior integral
    :param int/float b: limite superior integral
    :param int m: numero de puntos
    """
    x = sy.symbols('x')
    fx = sy.sympify(f)

    signo = 1
    if (a > b):  # inversion de signo en caso de a > b
        a_temp = a
        a = b
        b = a_temp
        signo = -1

    if es_continua(fx, x, a, b):  # evaluar continuidad de funcion
        if (a != b):
            I = 0
            error = 0
            h = (b-a)/(m-1)
            x0 = a
            for i in range(m-1):
                xi = x0 + i*h
                xii = x0 + (i+1)*h
                result_trapecio = trapecio(f, xi, xii)
                I += result_trapecio[0]
                error += result_trapecio[1]
            return (I*signo, error)
        else:
            return (0, 0)  # caso a == b
    else:
        return "La funcion no es continua en el intervalo [" + str(a) + "," + str(b) + "]"


def simpson_compuesto(f, a, b, m):
    """
    Ejecuta la aproximacion de la integral de la funcion f en el intervalo [a, b] por el
    metodo de Simpson compuesto

    :param str f: funcion en formato de string
    :param int/float a: limite inferior integral
    :param int/float b: limite superior integral
    :param int m: numero de puntos
    """
    x = sy.symbols('x')
    fx = sy.sympify(f)

    signo = 1
    if (a > b):  # inversion de signo en caso de a > b
        a_temp = a
        a = b
        b = a_temp
        signo = -1

    if es_continua(fx, x, a, b):  # evaluar continuidad de funcion
        if (a != b):
            I = 0
            error = 0
            h = (b-a)/(m-1)
            x0 = a
            for i in range(m-1):
                xi = x0 + i*h
                xii = x0 + (i+1)*h
                result_simpson = simpson(f, xi, xii)
                I += result_simpson[0]
                error += result_simpson[1]
            return (I*signo, error)
        else:
            return (0, 0)  # caso a == b
    else:
        return "La funcion no es continua en el intervalo [" + str(a) + "," + str(b) + "]"


def cuadraturas_gausseanas(f, n):
    x = sy.symbols('x')
    fx = sy.sympify(f)
    I = 0
    xi = [[0], [0.57735, -0.57735], [0.77459, -0.77459, 0], [0.86113, -0.86113, 0.33998, -0.33998],
          [0.90617, -0.90617, 0.53846, -0.53846, 0], [0.53846, -
                                                      0.53846, 0.6612, -0.6612, 0.23861, -0.23861],
          [0.9491, -0.9491, 0.74153, -0.74153, 0.40584, -0.40584, 0],
          [0.96028, -0.96028, 0.79666, -0.79666,
              0.52553, -0.52553, 0.18343, -0.18343],
          [0.96816, -0.96816, 0.83603, -0.83603,
              0.61337, -0.61337, 0.32425, -0.32425, 0],
          [0.9739, -0.9739, 0.86506, -0.86506, 0.6794, -0.6794, 0.43339, -0.43339, 0.14887, -0.14887]]
    wi = [[2], [1, 1], [0.55556, 0.55556, 0.88889], [0.34785, 0.34785, 0.65214, 0.65214],
          [0.23692, 0.23692, 0.47862, 0.47862, 0.56889], [
              0.17132, 0.17132, 0.36076, 0.36076, 0.46791, 0.46791],
          [0.12948, 0.12948, 0.2797, 0.2797, 0.38183, 0.38183, 0.41795],
          [0.10122, 0.10122, 0.22238, 0.22238, 0.3137, 0.3137, 0.36268, 0.36268],
          [0.08127, 0.08127, 0.18064, 0.18064, 0.26061,
              0.26061, 0.31234, 0.31234, 0.33023],
          [0.06667, 0.06667, 0.14945, 0.14945, 0.21908, 0.21908, 0.26926, 0.26926, 0.29552, 0.29552]]
    for i in range(len(xi[n-1])):
        I = I + (wi[n-1][i]*sy.N(fx.subs(x, xi[n-1][i])))
    return I


def cuadraturas_gausseanas_general(f, a, b, n):
    """
    Ejecuta la aproximacion de la integral de la funcion f en el intervalo [a, b] por el
    metodo de Cuadraturas Gauss-Legendre

    :param str f: funcion en formato de string
    :param int/float a: limite inferior integral
    :param int/float b: limite superior integral
    :param int n: grado polinomio Legendre
    """
    x = sy.symbols('x')
    fx = sy.sympify(f)

    signo = 1
    if (a > b):  # inversion de signo en caso de a > b
        a_temp = a
        a = b
        b = a_temp
        signo = -1

    if valida_orden(n):
        if es_continua(fx, x, a, b):  # evaluar continuidad de funcion
            if (a != b):
                gx = (b-a)/2 * (fx.subs(x, ((b-a)*x + b + a)/2))
                I = cuadraturas_gausseanas(gx, n)
                return I*signo
            else:
                return (0, 0)  # caso a == b
        else:
            return "La funcion no es continua en el intervalo [" + str(a) + "," + str(b) + "]"
    else:
        return "Orden 'n' debe ser un numero entero: 1 <= n <= 10"


def valida_orden(n):
    if isinstance(n, int):
        if (1 <= n <= 10):
            return True
    return False


def es_continua(f, x, a, b):
    """
    Determina si la funcion es continua en el intervalo [a,b]

    :param Sympy.sympify f
    :param Sympy.symbols x
    :param float a
    :param float b
    """
    dominio_funcion = continuous_domain(f, x, sy.S.Reals)
    intervalo = sy.Interval(a, b)
    interseccion = sy.Intersection(dominio_funcion, intervalo)
    return interseccion == intervalo


# Main
f_x = "1/x"
a, b = 1, 2

# Simples
result_trapecio = trapecio(f_x, a, b)
print("Metodo trapecio: I = " +
      str(result_trapecio[0]) + "\tError = " + str(result_trapecio[1]) + "\n")

result_simpson = simpson(f_x, a, b)
print("Metodo de Simpson: I = " +
      str(result_simpson[0]) + "\tError = " + str(result_simpson[1]) + "\n")

result_boole = regla_boole(f_x, a, b)
print("Metodo de Boole: I = " +
      str(result_boole[0]) + "\tError = " + str(result_boole[1]) + "\n")

# Compuestos
m = 10
result_trapecio_comp = trapecio_compuesto(f_x, a, b, m)
print("Metodo de trapecio compuesto: I = " +
      str(result_trapecio_comp[0]) + "\tError = " + str(result_trapecio_comp[1]) + "\n")

result_simpson_comp = simpson_compuesto(f_x, a, b, m)
print("Metodo de Simpson compuesto: I = " +
      str(result_simpson_comp[0]) + "\tError = " + str(result_simpson_comp[1]) + "\n")

orden = 8
print("Metodo de cuadraturas guasseanas: I = " +
      str(cuadraturas_gausseanas_general(f_x, a, b, orden)))
