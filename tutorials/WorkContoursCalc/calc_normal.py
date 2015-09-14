## Calculate normals of Yield locus given in terms of their cartesian coordinates (x,y)

import numpy as np
import matplotlib.pyplot as plt


def main(fn='./IFsteel/dat/Hill48.txt',locus_xy=None,ths=[90,45,0]):
    thetas = np.array(ths)*np.pi/180.
    X,Y,normals = read(fn,locus_xy)
    X_at_th = np.interp(thetas,normals,X)
    Y_at_th = np.interp(thetas,normals,Y)
    return X_at_th, Y_at_th

def read(fn='./IFsteel/dat/Hill48.txt',locus_xy=None):
    if type(locus_xy)==type(None):
        locus_xy = np.loadtxt(fn,skiprows=1,delimiter=',').T
    x,y = locus_xy
    flt = []
    for i in xrange(len(x)):
        if x[i]>0 and y[i]>0:
            flt.append(True)
        else: flt.append(False)
    flt=np.array(flt)
    x=x[flt]; y=y[flt]
    r,th=xy2rt(x,y)
    diffx = np.diff(x)
    diffy = np.diff(y)

    slopes = diffx/diffy
    ths=[]
    for i in xrange(len(x)-1):
        th_ = np.arctan2(diffy[i],diffx[i])+np.pi/2.
        th_ = th_180(th_)
        ths.append(th_)
    th_av=[]
    for i in xrange(len(ths)-1):
        th_av.append((ths[i]+ths[i+1])/2.)
    th_av=np.array(th_av)
    X = x[1:-1]
    Y = y[1:-1]
    TH = th[1:-1]
    for i in xrange(len(X)):
        x,y = X[i],Y[i]
        th = th_av[i]
        r=0.05
        dx=np.cos(th)*r
        dy=np.sin(th)*r
        r_old = np.sqrt(x**2+y**2)
        r_new = np.sqrt((x+dx)**2+(y+dy)**2)
        if r_new<r_old:
            th = th_av[i] + np.pi
            if th>np.pi:
                th = -(2*np.pi-th)
            th_av[i]=th
            dx=np.cos(th)*r
            dy=np.sin(th)*r
    return X,Y,th_av

def ths_180(ths):
    for i in xrange(len(ths)):
        ths[i] = th_180(ths[i])
    return ths[i]

def th_180(th):
    if th>np.pi:
        th = -(2*np.pi-th)
    if th<-np.pi:
        th = th + np.pi*2
    return th
    
def xy2rt(x,y):
    th=np.arctan2(y,x)
    r=np.sqrt(x**2+y**2)
    return r, th


def find_stress(x,y,*ths):
    TH   = np.arctan2(y,x)
    _xs_ = np.interp(ths,TH,x)
    _ys_ = np.interp(ths,TH,y)
    return _xs_, _ys_

def ex(thetas=[0,45,90]):
    fig=plt.figure(1)
    fig.clf()
    ax=fig.add_subplot(111)
    ax.set_aspect('equal')

    ths = np.linspace(-np.pi,np.pi,100)
    x,y = np.cos(ths), np.sin(ths)
    ax.plot(x,y,'.')
    thetas=np.array(thetas)*np.pi/180.
    xs,ys=find_stress(x,y,*thetas)

    print 'xs:','-'*20
    print xs
    print 'ys:','-'*20
    print ys
    ax.plot(xs,ys,'o')
    ax.set_xlim(-1.1,1.1)
    ax.set_ylim(-1.1,1.1)


    ## Actual
    x,y = np.cos(thetas),np.sin(thetas)
    ax.plot(x,y,'+')
    
