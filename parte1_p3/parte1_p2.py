"""
Fabián Crawford Barquero
Irene Muñoz Castro
Luis Morales Rodríguez
Steven Badilla Soto
Adrián Trejos Salazar
"""

import numpy as np
from sympy import *
import sympy as sp

def simpson_compuesto(function, a, b, n):
     
    """
    Ejecuta la aproximacion de la integral de la funcion f en el intervalo [a, b] por el
    metodo de Simpson compuesto

    :param str f: funcion en formato de string
    :param int/float a: limite inferior integral
    :param int/float b: limite superior integral
    :param int m: numero de puntos
    """
    f = sp.lambdify([sp.Symbol('x')], function)
    xk=a
    sumImpar=0
    sumPar=0
    k=1
    h=(b-a)/(n-1)
    par=0
    impar=0
    while k<n:
        xk=a+k*h
        if k!=n-1 and k%2 == 0:
            sumPar= sumPar + f(xk)
            par=par+1
        else:
            if k!=n-1:
                sumImpar= sumImpar + f(xk)
                impar=impar+1

        k=k+1

    #Aproximación de la integral
    I=(h/3)*(f(a) + 4*sumImpar + 2*sumPar + f(b))

    #Calculando la cota de Error
    x = Symbol('x')
    f1=diff(function,x)
    f2=diff(f1,x)
    f3=diff(f2,x)
    f4=diff(f3,x)
    f4=sp.lambdify([sp.Symbol('x')], f4)
    val1=f4(a)
    val2=f4(b)
    if abs(val1)>abs(val2):
        M=abs(val1)
    else: M= abs(val2)
    er=((b-a)*M*h**4)/180

    return [I , er]


def simpson_simple(f, a, b):
    """
    Ejecuta la aproximacion de la integral de la funcion f en el intervalo [a, b] por el
    metodo de Simpson

    :param str f: funcion en formato de string
    :param int/float a: limite inferior integral
    :param int/float b: limite superior integral
    """
    x = Symbol('x')
    function = sympify(f)
    h = (b - a) / 2
    x1 = (a + b) / 2

    I = ((h / 3) * (function.subs(x, a) + (4 * function.subs(x, x1)) + function.subs(x, b)))
    I = I.evalf()
    fourthDeriv = diff(f, x, x, x, x)
    if(fourthDeriv.subs(x, a) > fourthDeriv.subs(x, b)):
        er = ((h**5) / 90) * abs(fourthDeriv.subs(x, b))
        er = er.evalf()

    else:
        er = ((h**5) / 90) * abs(fourthDeriv.subs(x, a))
        er = er.evalf()

    return [I, er]


def boole(f, a, b):
    """
    Ejecuta la aproximacion de la integral de la funcion f en el intervalo [a, b] por el
    metodo de la Regla de Boole

    :param str f: funcion en formato de string
    :param int/float a: limite inferior integral
    :param int/float b: limite superior integral
    """
    x = Symbol('x')
    function = sympify(f)
    evalValue = []
    h = (b - a) / 4
    for i in range(0, 5):
        evalValue.append(a + (i * h))
    I = ((2 * h) / (45)) * ((7 * function.subs(x, evalValue[0])) +
                            (32 * function.subs(x, evalValue[1])) +
                            (12 * function.subs(x, evalValue[2])) +
                            (32 * function.subs(x, evalValue[3])) +
                            (7 * function.subs(x, evalValue[4])))
    sixth_deriv = diff(function, x, x, x, x, x, x)
    if(sixth_deriv.subs(x, a) > sixth_deriv.subs(x, b)):
        er = -(8*(h**7) / 945) * abs(sixth_deriv.subs(x, b))
        er = er.evalf()
    else:
        er = -(8*(h**7) / 945) * abs(sixth_deriv.subs(x, a))
        er = er.evalf()

    return [I, er]



def trapecio_compuesto(function, a, b, n):
    """
    Ejecuta la aproximacion de la integral de la funcion f en el intervalo [a, b] por el
    metodo de Simpson compuesto

    :param str f: funcion en formato de string
    :param int/float a: limite inferior integral
    :param int/float b: limite superior integral
    :param int m: numero de puntos
    """
    f = sp.lambdify([sp.Symbol('x')], function)
    xk=a
    sum=0
    k=1
    h=(b-a)/(n-1)

    while k<n:
        xk=a+k*h
        if k!=n-1:
            sum= sum + f(xk)
        k=k+1

    #Aproximación de la integral
    I=(h/2)*(f(a) + 2*sum + f(b))

    #Calculando la cota de Error
    x = Symbol('x')
    f1=diff(function,x)
    f2=diff(f1,x)
    f2=sp.lambdify([sp.Symbol('x')], f2)
    val1=f2(a)
    val2=f2(b)
    if abs(val1)>abs(val2):
        M=abs(val1)
    else: M= abs(val2)
    er=((b-a)*M*h**2)/12

    return [I , er]


def trapecio_simple(f, a, b):
    """
    Ejecuta la aproximacion de la integral de la funcion f en el intervalo [a, b] por el
    metodo del trapecio

    :param str f: funcion en formato de string
    :param int/float a: limite inferior integral
    :param int/float b: limite superior integral
    """
    x = Symbol('x')
    function = sympify(f)
    h = b - a
    I = (( h / 2) * (function.subs(x, a) + function.subs(x, b)))
    I = I.evalf()
    functionSecondDeriv = diff(function, x, x)
    if(functionSecondDeriv.subs(x, a) > functionSecondDeriv.subs(x, b)):
        er = ((h**3) / 12) * abs(functionSecondDeriv.subs(x, b))
        er = er.evalf()
    else:
        er = ((h**3) / 12) * abs(functionSecondDeriv.subs(x, a))
        er = er.evalf()

    return [I, er]


def cuadraturas_gaussianas(f, a, b, n):
    """
    Ejecuta la aproximacion de la integral de la funcion f en el intervalo [a, b] por el
    metodo de Cuadraturas Gauss-Legendre

    :param str f: funcion en formato de string
    :param int/float a: limite inferior integral
    :param int/float b: limite superior integral
    :param int n: grado polinomio Legendre
    """
    x = Symbol('x')
    ceros=ceros_pol(n)[0]
    p=ceros_pol(n)[1]
    w=pesos(ceros, p, n)
    fun=sp.sympify(f)
    sum=0
    k=0
    for item in ceros:
        argument=((b-a)*item+b+a)/2
        sum = sum + w[k]*fun.subs(x,argument)
        k=k+1
    I=((b-a)/2)*sum

    #Cota de error
    x = Symbol('x')
    fdiff=diff(f,x, 2*n)
    fun=sp.lambdify([sp.Symbol('x')], fdiff)
    val1=fun(a)
    val2=fun(b)
    if abs(val1)>abs(val2):
        M=abs(val1)
    else: M= abs(val2)
    error=( ( ((b-a)**(2*n+1)) *np.math.factorial(n)**4)/( (2*n+1)*(np.math.factorial(2*n))**3 ) )*M

    return [I ,error]

def ceros_pol(n):
    x = Symbol('x')
    q=diff((x**2-1)**n, x, n)
    p=(1/(np.math.factorial(n)*2**n))*q
    ceros=solve(p,x)
    return [ceros, p]

def pesos(ceros, p, n):
    x = Symbol('x')
    pd=diff(p,x)
    pd=lambdify(sp.Symbol('x'), pd)

    w=[]
    for item in ceros:
        wk=2/((1-item**2)*(pd(item))**2)
        w.append(wk)
    return w


def validar_simbolos(simbolos):
    i = 0
    while i < len(simbolos):
        if(str(simbolos[i]) != "x"):
            return False
        else:
            i = i + 1
    return True

def validar_funciones(argumentos):
    i = 0
    funciones = ["sin","pi", "cos", "tan","cot", "sec","csc","sinc",
                 "asin","acos","atan","acot","asec","acsc","atan2",
                 "sinh", "cosh", "tanh","coth", "sech","csch", "sqrt",
                 "asinh","acosh","atanh","acoth","asech","acsch","log", "exp"]
    lista = []
    i = 0
    while i < len(argumentos):
        j = 0
        funcion_actual = str(argumentos[i])
        funcion = ""
        while j < len(str(funcion_actual)):
            if(funcion_actual[j] != '('):
                funcion = funcion + str(funcion_actual[j])
                j = j + 1
            else:
                break
        i = i + 1
        lista = lista + [str(funcion)]

    i = 0
    while i < len(lista):
        if lista[i] == '':
            del lista[i]
        i = i + 1

    for x in lista:
        if x in funciones:
            continue
        else:
            return False
    return True
