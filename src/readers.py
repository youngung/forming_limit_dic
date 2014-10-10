### collection of some readers
import numpy as np
from os import sep
def read_disc(fn='BB/2012JUL11/node_data.csv',ndat_col=13,
              ixy_flip=False):
    """
    Read dat extracted from VIC3D
    """
    # import pandas as pd
    # pd.read_csv()
    dat=np.genfromtxt(fn,delimiter=',',skiprows=1)
    shp = dat.shape
    # print shp
    # print ndat_col
    # print float(shp[1])/float(ndat_col)
    ## raise IOError
    dat=np.reshape(dat, (shp[0],shp[1]/ndat_col,ndat_col))

    if ixy_flip:
        ## xy coordinate
        dum = dat[:,:,0][::].copy() ## x coordinate
        dat[:,:,0] = dat[:,:,1][::].copy() ## y cooridnate
        dat[:,:,1] = dum[::].copy()

        ## strain coordinate
        dum= dat[:,:,3][::].copy()
        dat[:,:,3] = dat[:,:,4][::].copy()
        dat[:,:,4] = dum[::].copy()

        ## strain rate coordinate
        dum = dat[:,:,9][::].copy()
        dat[:,:,9] = dat[:,:,10][::].copy()
        dat[:,:,10]= dum[::].copy()
    return dat

def simple_3col(
        fn='../dat/IFsteel/wXRD/BB/BB_'+\
        'extracted_node_dat_revisited_201'+\
        '4OCT09_sepsteelB1-2_BB_2012Jul11-0008_0.csv'):
    """
    Read 3col dat (exx,eyy,exy)
    """
    dat=np.genfromtxt(fn,delimiter=',').T
    dat=dat[:3]
    exx,eyy,exy = dat
    ezz = -exx-eyy
    ezz = abs(ezz)
    ar = []
    for i in range(len(ezz)):
        if np.isnan(ezz[i]): ar.append(-1000)
        else: ar.append(ezz[i])

    ind = np.argmax(ar)
    return exx[ind],eyy[ind],-ezz[ind]

def write_dic_csvs(path='../dat/IFsteel/wXRD/BB'):
    """
    Find max exx+eyy and write it to 'DIC_ezz.txt'

    Arguments
    =========
    path='../dat/IFsteel/wXRD/BB'
    """
    from glob import glob
    fns = glob('%s%s%s'%(path,sep,'*sep*2012Jul*.csv'))
    fns.sort()
    dic_inds=[]

    Ezz = []
    Exx = []
    Eyy = []
    for i in range(len(fns)):
        fn  = fns[i]
        exx,eyy,ezz = simple_3col(fn)
        dic_ind = int(fn.split(sep)[-1].split('-')[-1].\
                      split('.csv')[0].split('_')[0])
        dic_inds.append(dic_ind)
        print dic_ind, fn

        # if i==0:
        #     Exx.append(exx)
        #     Eyy.append(eyy)
        #     Ezz.append(ezz)
        try:
            if dic_ind==dic_inds[i-1]+1:
                Exx.append((exx_old + exx )/2.)
                Eyy.append((eyy_old + eyy )/2.)
                Ezz.append((ezz_old + ezz )/2.)
        except: pass
        exx_old = exx
        eyy_old = eyy
        ezz_old = ezz

    Exx.append(exx_old)
    Eyy.append(eyy_old)
    Ezz.append(ezz_old)

    print len(Exx)

    f=open('DIC_eps.txt','w')
    f.write('%7s %7s %7s\n'%('Exx','Eyy','Ezz'))
    for i in range(len(Exx)):
        f.write('%2i  %+.4f %+.4f %+.4f\n'%(i+1,Exx[i],Eyy[i],Ezz[i]))
    f.close()
