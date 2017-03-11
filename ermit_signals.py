#!/usr/bin/python3
#-*- coding: utf-8 -*-

import import_lib as il
#вычисление производных
def dfdx_n(x1,n, flag=True):
    x = il.sym.Symbol('x')
    if(flag):
        f = il.sym.exp(-x*x/2.0)
    else:
        f = pow(x,n)*il.sym.exp(-x/1.0)
    if(n==False):
        return f.subs(x,x1)
    else:
        diff_f = il.sym.diff(f, x,n)
        return diff_f.subs(x,x1)
        
#вычисление полинома Эрмита  
def hermitov_poly(x1,n):    
    result = pow(-1,n)*il.m.exp(pow(x1,2)/2.0)*dfdx_n(x1,n)
    return result

#вычисление полинома и ряда Лаггера
def lager_poly(x1, n):
    k = (1./il.m.factorial(n))*il.m.exp(x1)
    m = dfdx_n(x1,n, flag=False)
    return k*m
    
def lager_row(x1,n):
    sum,k = 0,0
    if(n==k):
        sum = 1
        return sum
    else:
        while(k<=n):
            sum+=il.m.factorial(n)*pow(x1,k)*pow(-1,k)/float(il.m.factorial(n-k)*il.m.factorial(k)*il.m.factorial(k))
            k+=1
    return sum
#вычисление классов сигнала инвариантных к повороту плоскости неопределенности
def signals_u(t,n):
    K = pow(2,0.25)/pow(il.m.factorial(n),0.5)
    M = il.m.exp(-il.m.pi*pow(t,2))
    return K*M*hermitov_poly(2*pow(il.m.pi,0.5)*t,n)

#преобразование Фурье от исходного сигнала
def spectr_u(f,n):
    return pow((-1j),n)*signals_u(f,n)

def Ermit(t,f,n):
    return il.m.exp(-il.m.pi*(t**2 + f**2)/2.)*lager_row(il.m.pi*(t**2 + f**2),n)
if __name__ == "__main__":
    X = [x for x in il.np.arange(0,5,0.2)]
    Y = [f for f in il.np.arange(0,5,0.2)]
    xgrid, ygrid = il.np.meshgrid(X, Y)
    Z = [[Ermit(t,f,10) for t in X ] for f in Y]
    import pylab
    from mpl_toolkits.mplot3d import Axes3D
    fig = pylab.figure()
    axes = Axes3D(fig)
    axes.plot_surface(xgrid, ygrid,Z, rstride=1, cstride=1)
    pylab.show()    
"""
    
    Y = [signals_u(t,4) for t in T]
    MAX_Y = max(Y)
    Y = [round(y/MAX_Y,2) for y in Y]
    Z = il.np.fft.fft(abs(il.np.array(Y)))
    il.plt.plot(T,Z)
    il.plt.grid()
    il.plt.show()
"""
