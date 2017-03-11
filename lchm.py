#!/usr/bin/python3
#-*- coding: utf-8 -*-
import import_lib as il
import cmath as cm
#ЛЧМ с V-образной ЧМ
T = 10e-3
delta_f = 1e+4
b = 2*T*delta_f
print "b = ", b
k = 1/pow(2*T,0.5)
def signal_u1(t):
    return k*cm.exp(1j*b*t*t)
def signal_u2(t):
    return k*cm.exp(-1j*b*t*t)
def signal_u(t):
    if -T<t<T:
        return signal_u1(t)+signal_u2(t)
    else:
        return 0
def X11(t,f):
    k1 = 0.5*(1-abs(t)/T)*il.m.sin(il.m.pi*(f-b*t/il.m.pi)*(T - abs(t)))
    k2 = cm.exp(1j*(il.m.pi*f*(T+t)-b*t*T))
    if (il.m.pi*(f-b*t/il.m.pi)*(T - abs(t)))==0:
        return 1
    else:
        return k1*k2/(il.m.pi*(f-b*t/il.m.pi)*(T - abs(t)))
def X22(t,f):
    return X11(t,-f)*cm.exp(1j*2*il.m.pi*f*t)
def Lchm_V(t,f):
    return X11(t,f) + X22(t,f)
    
if __name__ == "__main__":
    time = [t for t in il.np.arange(-6,6,100*T)]
    freq = [f for f in il.np.arange(-2e3,2e3,1/T)]
    xgrid, ygrid = il.np.meshgrid(time, freq)
    Z = [[abs(Lchm_V(t,f)) for t in time ] for f in freq]
    import pylab
    from mpl_toolkits.mplot3d import Axes3D
    fig = pylab.figure()
    axes = Axes3D(fig)
    axes.plot_surface(xgrid, ygrid,Z, rstride=1, cstride=1)
    pylab.show()