# rk4.py · TFG: el sistema de Lorenz con Euler y RK4 de paso fijo,
# y la comprobación de su orden de error.
# Versión del 2026-10-28, sin terminar: lo que falta, al final del fichero.
# Requisitos: Python 3.11 y NumPy (usado: 2.4.6).

import numpy as np


# ===== 1. El sistema =====

# Sistema de Lorenz: la derivada del estado (x, y, z).
# Entra: x, el estado, array de forma (3,); r, el parámetro de control
# (28 si no se da). Sale: array nuevo de forma (3,).
# sigma = 10 y b = 8/3, fijos.
def lorenz(x,r=None):
    
    sigma=10.0
    b=8/3
    
    if r is None:
        r=28.0
    
    xpunto=sigma*(x[1]-x[0])
    ypunto=r*x[0] -x[1] -x[0]*x[2]
    zpunto=x[0]*x[1] -b*x[2]
    
    return np.array([xpunto,ypunto,zpunto])


# ===== 2. Los métodos de un paso =====

# Un paso de Euler: estado nuevo = x + h*f(x).
# Entra: f, función que da la derivada de un estado; x, el estado,
# array de forma (d,); h, el paso.
# Sale: el estado un tiempo h después, array nuevo de forma (d,);
# x no cambia.
def paso_euler(f,x,h):
    return x + f(x)*h


# Un paso de RK4: cuatro pendientes (al principio, dos en el punto
# medio y al final) con pesos 1, 2, 2, 1.
# Entra y sale lo mismo que en paso_euler.
def paso_rk4(f,x,h):
    k1= f(x)
    k2= f(x+k1*(h/2))
    k3= f(x+k2*(h/2))
    k4= f(x+k3*h)
    return x + (k1+2*k2+2*k3+k4)*(h/6)


# ===== 3. El integrador de paso fijo =====

# n pasos de tamaño h = T/n con el método g, desde x0.
# Entra: g, paso_euler o paso_rk4; f, como en los pasos; x0, el estado
# inicial, array de forma (d,); T, el tiempo final; n, el número de
# pasos (entero).
# Sale: el estado en T, array de forma (d,); x0 no cambia.
# La trayectoria, de forma (n+1, d), donde cada fila es un estado en cada instante, se guarda pero no se devuelve.
def integrador(g,f,x0,T,n):
    
    h=T/n
    x=x0
    trayectoria=np.zeros((n+1,len(x0)))
    trayectoria[0]= x0
    
    for i in range(n):
        x=g(f,x,h)
        trayectoria[i+1]=x
    
    return x


# ===== 4. La medida del orden =====

# Error en T del método g para cada número de pasos de la lista m.
# Entra: g, f, x0 y T, como en integrador; m, lista de números de pasos;
# xT, el estado de referencia en T, array de forma (d,).
# Sale: error (la norma euclídea de xT menos el estado calculado) y haches (el paso T/n), dos arrays de forma (len(m),).

def errores(g,f,x0,T,m,xT):
    
    x=[]
    haches=[]
    error=[]
    
    for i in m:
        haches.append(T/i)
        y=integrador(g,f,x0,T,i)
        x.append(y)
    
    x=np.array(x)
    
    for i in range(len(x)):
        error_i=np.linalg.norm(xT-x[i])
        error.append(error_i)
    
    error=np.array(error)
    haches=np.array(haches)
    
    return error, haches


# Pendiente de la recta log-log entre (h1, E1) y (h2, E2):
# log(E2/E1) / log(h2/h1).
# ValueError si algún valor es menor o igual que 0, o si h1 = h2.
def pendiente(h1,E1,h2,E2):
    
    valores=[h1, E1, h2, E2]
    verdad=[]
    
    for i in range(4):
        verdad.append(valores[i]<=0)
        
    if any(verdad) == True :
        raise ValueError("valores negativos o iguales a 0")
    
    elif abs(h1-h2)<1e-12:
        raise ValueError("valores h iguales")
        
    else:
        mlog=np.log((E2/E1))/np.log((h2/h1))
        
        return mlog


# Pendientes locales entre puntos consecutivos de (h, err)
# Entra: err y h, arrays de forma (k,). Sale: array de forma (k-1,).
def pendienteslocales(err,h):
    diferenciasr=np.diff(np.log(err))
    diferenciash=np.diff(np.log(h))
    return diferenciasr/diferenciash


#------------------------------ no terminado ------------------------------
# Falta la comprobación del orden en Lorenz. 
#   1. elegir el tiempo final T, el punto de partida y los n de cada método;
#   2. calcular la solución de referencia y estimar su error;
#   3. errores y pendientes locales de Euler y RK4, en una tabla
#   4. la figura log-log del error frente a h
# Por decidir: si integrador devuelve también la trayectoria o deja de
# guardarla.