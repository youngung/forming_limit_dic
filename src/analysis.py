### read DIC profile csv file and carry out analysis
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os, time
from MP.lib.axes_label import __deco_fld__ as deco_fld
from MP.lib.mpl_lib import wide_fig as wf
from MP.lib.mpl_lib import rm_all_lab as ral
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_pdf import PdfPages
from readers import read_disc
from MP.lib import mpl_lib
from MP import ssort, progress_bar
from numpy import random
uet = progress_bar.update_elapsed_time
sort=ssort.shellSort

### data analyzed by and extracted from a newer Vic3D
def main_disc(fn='BB/2011NOV10/node_disc_36.csv',
              dist_bb=[10],
              dist_rb=[5],
              fang=90.,ixy_flip=False,ioff=True,
              ntot_col=14,
              icol_ezz=12,
              fig=None,
              iopt=0,
              fout=None):
    """
    fn       = file name: " /dum/[strainpath]/[date].csv "
    dist_bb  = distances from the critical point
    dist_rb  = distances from the critical point
    fang     = fracture angle
    ixy_flip = flag if x-y cooridnate should be switched
    ioff     = interactive mode switch for MPL
    ntot_col = total number of column belonging to a dat block
    icol_ezz = 12
    fig      = None
    iopt     = 0
    fout     = None

    ## below is the column information for a block
    icol:     0 1 2 3   4   5   6  7  8  9      10     11     12
    variable: x,y,z,exx,eyy,exy,e1,e2,vm,exxdot,eyydot,exydot,e33

    ## some other forms
    icol:     0 1 2 3   4   5   6  7    8       9     10    11   12
    variable: x,y,z,exx,eyy,exy,e1,e2,exxdot,eyydot,exydot,ezz,ezzdot
    """

    ## output file name
    strain_path = fn.split(os.sep)[1]
    date = fn.split(os.sep)[2].split('.')[0]
    fn_out = '%s_%s'%(strain_path,date)
    pdf_all = PdfPages('%s.pdf'%fn_out)

    ## MPL interactive mode switch
    if ioff: plt.ioff()
    else: plt.ion()

    ## Read data from VIC-3d software
    dat=read_disc(fn=fn,ndat_col=ntot_col,
                  ixy_flip=ixy_flip)

    ## when using debugging: reduce the data sets.
    nst,nxy,ncol=dat.shape
    # dat = dat[:,random.choice(nxy, nxy/10),:]
    # nst,nxy,ncol=dat.shape

    ## find the index for a reduced set of
    ## xy coordinates.
    # dat_r = dat[:,random.choice(nxy, nxy/20),:]
    dat_r = dat[::].copy()
    n_image = find_indx_at_fracture(dat_r,iopt=iopt)

    ## n_image-10 to n_image
    fig_p   = plt.figure(figsize=(36,24))
    fig_ezz = plt.figure(figsize=(36,24))
    gs=gridspec.GridSpec(
        40,40,wspace=4.0,hspace=2.0)
    X=[];Y=[];Z=[];X_=[];Y_=[];Z_=[]
    inds=[];m_=[];M_=[];Norm=[];Norm_=[]
    for i in xrange(10):
        mn =+np.inf; mx =-np.inf
        mn_=+np.inf; mx_=-np.inf
        ind = n_image-10+i+1
        inds.append(ind)
        xs,ys,zs = analysis_xy_disc_e33dot(
            dat=dat,istp=ind,icol=12)
        xs_,ys_,ezz = analysis_xy_disc_ezz(
            dat=dat,istp=ind)
        zs=np.log10(zs)
        for j in xrange(len(zs)):
            if zs[j]==-np.inf or zs[j]==np.inf: zs[j]=np.nan
            if mx<zs[j]:   mx=zs[j]
            elif mn>zs[j]: mn=zs[j]
        for j in xrange(len(ezz)):
            if ezz[j]==-np.inf or ezz[j]==np.inf: ezz[j]=np.nan
            if mx_<ezz[j]:   mx_=ezz[j]
            elif mn_>ezz[j]: mn_=ezz[j]

        norm=mpl.colors.Normalize(vmin=mn,vmax=mx)
        cmap,m=mpl_lib.norm_cmap(mn=mn,mx=mx,cm_name='brg')
        m_.append(m)
        norm_=mpl.colors.Normalize(vmin=mn_,vmax=mx_)
        cmap,m=mpl_lib.norm_cmap(mn=mn_,mx=mx_,cm_name='brg')
        M_.append(m)

        X.append(xs);Y.append(ys);Z.append(zs)
        X_.append(xs_);Y_.append(ys_);Z_.append(ezz)
        Norm.append(norm);Norm_.append(norm_)

    # mx=max(Z);mn=min(Z)
    print

    ## Construct figures
    axps = []; axcbs = [] ; n_ax = 0
    ax_ezzs=[]; axcbs_ezz=[]
    for ix in xrange(4):
        for iy in xrange(4):
            if n_ax<10:
                _ax_=fig_p.add_subplot(
                    gs[10*ix:10*(ix+1),10*iy  :10*iy+8])
                _axcb_=fig_p.add_subplot(
                    gs[10*ix:10*(ix+1),10*iy+8:10*iy+9])
                axps.append(_ax_); axcbs.append(_axcb_)
                _ax_=fig_ezz.add_subplot(
                    gs[10*ix:10*(ix+1),10*iy  :10*iy+8])
                _axcb_=fig_ezz.add_subplot(
                    gs[10*ix:10*(ix+1),10*iy+8:10*iy+9])
                ax_ezzs.append(_ax_); axcbs_ezz.append(_axcb_)
            n_ax=n_ax+1

    ## Plot
    t0 = time.time()
    for i in xrange(10):
        axp = axps[i]
        xs,ys,zs = X[i],Y[i],Z[i]
        cs = []
        for j in xrange(len(xs)):
            cs.append(m_[i].to_rgba(zs[j]))
        uet(time.time()-t0,head='time spent for plotting')
        axp.scatter(xs,ys,marker='.',color=cs)
        axp.set_title('image #: %i'%inds[i])

        ax_ezz=ax_ezzs[i]
        xs,ys,zs = X_[i],Y_[i],Z_[i]
        cs = []
        ## color for each dot
        for j in xrange(len(xs)):
            cs.append(M_[i].to_rgba(zs[j]))
        uet(time.time()-t0,head='time spent for plotting')
        ax_ezz.scatter(xs,ys,marker='.',color=cs)
        ax_ezz.set_title('image #: %i'%inds[i])
        mpl_lib.add_cb(ax=axcbs[i],filled=True,
                       format='%7.2e',norm=Norm[i],
                       ylab=r'$\dot{\bar{E}}_{33}$')
        mpl_lib.add_cb(ax=axcbs_ezz[i],filled=True,
                       format='%7.4f',norm=Norm_[i],
                       ylab=r'$\bar{E}_{33}$')

    print
    ral(axps)#fig_p.axes)
    ral(ax_ezzs)#fig_ezz.axes)
    fig_p.savefig('ezz_dot.pdf')
    fig_p.savefig('ezz_dot.png')
    fig_ezz.savefig('ezz.pdf')
    fig_ezz.savefig('ezz.png')
    plt.close(fig_p);plt.close(fig_ezz)
    uet(time.time()-t0,head='time spent for plotting')
    print

    # ## through thickness strain rate minimum
    fig1 = plt.figure(figsize=(6,4))
    gs   = gridspec.GridSpec(10,20,wspace=0.2,top=0.9,left=0.2,right=0.7)
    ax   = fig1.add_subplot(gs[:,0:12])
    axb  = fig1.add_subplot(gs[:,14:16]) ## color bar

    fig2,fig3,mdat,c=plot_xy_disc_e33dot(
        dat,istp=n_image,icol=icol_ezz,
        _fig_ezz_=fig1,fang=fang,dist_bb=dist_bb,
        dist_rb=dist_rb)
    ax.set_xlim(-40,40);ax.set_ylim(-40,40)
    ax.set_aspect('equal')
    ax.set_xlabel('X [mm'); ax.set_ylabel('Y [mm')

    col=['path','date','dist [mm]','Exx','Eyy','Std1','Std2','color']
    row=['%9s','%9s','%9i','%9.3f','%9.3f','%9.3f','%9.3f','%9s']
    header =''; fmt    =''
    for i in xrange(len(col)):
        header='%s %s'%(header,'%9s'%col[i])
        fmt   ='%s %s'%(fmt,row[i])
    print header
    if fout!=None: fout.write(header+'\n')
    db=np.array(dist_bb)
    dum = (db[::-1]).tolist()+[0]
    ndat = mdat.shape[-1]
    for i in xrange(ndat):
        dat = mdat[:,:,i]
        ex = dat[0,0]; ey = dat[0,1]
        exe= dat[1,0]; eye= dat[1,1]
        fig.axes[0].errorbar(ey,ex,xerr=eye,yerr=exe,fmt='.',
                             color=c[i])
        line= fmt%(strain_path,date,dum[i],ex,ey,exe,eye,c[i])
        print line
        if fout!=None: fout.write(line + '\n')

    pdf_all.savefig(fig1); pdf_all.savefig(fig2)
    pdf_all.savefig(fig3)
    deco_fld(ax=fig.axes[0],iopt=2,ft=15,iasp=False)
    fig.axes[0].set_xlim(-0.1,)
    fig.axes[0].set_ylim(-0.1,)
    plt.draw()
    pdf_all.close()
    ## if ioff: plt.close('all')
    return dat_r, fig1, fig2, fig3

def main_disc_progr(dat,iopt,n=5):
    nstp, nxy, ncol = dat.shape
    n_image = find_indx_at_fracture(dat,iopt=iopt)

def find_coordinate_index(
    dat_disc,ind_image=-1,xy=[0,0],is_max=True,
    icol=11,dist_filter=1.):
    """
    Find index that cooresponds
    to the given 'xy' coordinate

    Arguments
    =========
    dat_disc
    ind_image   = -1
    xy          = [0,0]
    is_max      = True
    icol        = 11
    dist_filter = 1.
    """
    x0,y0 = xy
    d     = dat[istp,:,:]
    x,y   = dat[istp,:,0],dat[istp,:,1]
    dist  = []; ind_ref = []
    for i in xrange(len(x)):
        if np.isnan(x0) or np.isnan(y0): pass
        else:
            d=np.sqrt((x0-x)**2 + (y0-y)**2)
            if d<dist_filter:
                dist.append(d)
                ind_ref[i]
    if len(ind_ref)==0: return -1
    dist_sorted, indc = sort(dist)
    return in_ref[indc[0]]

def find_indx_at_fracture(dat,iopt=0):
    """
    Find the DIC index giving the maximum
    thinining rate

    Arguments
    =========
    dat
    iopt = 0
    """
    nstp, nxy, ncol = dat.shape
    if iopt==0: idx=9; idy=10
    elif iopt==1: idx=8; idy=9
    exx_dot = dat[:,:,idx]
    eyy_dot = dat[:,:,idy]
    ezz_dot = exx_dot+eyy_dot  ## absolute value
    EZ=[]
    t0=time.time()
    for istp in xrange(nstp):
        inx = np.nanargmax(ezz_dot[istp,:])
        EZ.append(ezz_dot[istp,inx])
        uet(second=time.time()-t0,
            head='Elapsed time in finding indx')

    return np.argmax(EZ)


def analysis_xy_disc_ezz(
        dat,istp=-1):
    d=dat[istp,:,:]
    n=d.shape[0]
    exx=d[:,3]
    eyy=d[:,4]
    ezz=-exx-eyy
    for i in xrange(len(ezz)):
        if ezz[i]==np.inf or ezz[i]==-np.inf:
            ezz[i] = np.nan
    xs,ys,zs=[],[],[]

    for i in xrange(n):
        if not(np.isnan(ezz[i])):
            x,y=d[i,0],d[i,1]
            xs.append(x)
            ys.append(y)
    return xs,ys,ezz

def analysis_xy_disc_e33dot(
        dat,istp=-1,icol=11):
    """
    Conduct the analysis on the disc.

    1. Find the spot that gives the maximum 'thinning rate'
    2. provide 'min-max' range of thickness strain rate

    Arguments
    =========
    dat
    istp      = -1
    _fig_ezz_ = None
    icol      = 11
    """
    d=dat[istp,:,:]
    #dz=np.log10(abs(d[:,icol]))
    dz=abs(d[:,icol])
    for i in xrange(len(dz)):
        if dz[i]==np.inf or dz[i]==-np.inf:
            dz[i] = np.nan
    n=d.shape[0]
    xs,ys,zs=[],[],[]
    for i in xrange(n):
        x,y,z=d[i,0],d[i,1],-d[i,icol]
        if not(np.isnan(z)):
            xs.append(x)
            ys.append(y)
            zs.append(z)
    return xs,ys,zs


def plot_xy_disc_e33dot(
        dat,istp=-1,_fig_ezz_=None,icol=11,fang=90,
        dist_bb=[5,10], ## backbone
        dist_rb=[1,5]):
    """
    Extract DIC data from 'disc' area

    Data are analyzed in the grid
    made of backbone and ribs
    * The backbone grid is an analogy to the
     standard method using circular grid

    Arguments
    =========
    dat
    istp    = -1 ## DIC image index
    _fig_ezz_ = None
    icol    = column for ezz
    fang    = fracture angle


    dist_bb = backbone bb
    dist_rb = rib      rb
    """
    d = dat[istp,:,:]
    dz = np.log10(abs(d[:,icol]))
    mx = -np.inf; mn = +np.inf
    for i in xrange(len(dz)):
        if dz[i]==np.inf or dz[i]==-np.inf:
            dz[i] = np.nan
        if dz[i]>mx: mx = dz[i]
        if dz[i]<mn: mn = dz[i]

    norm=mpl.colors.Normalize(vmin=mn, vmax=mx)
    cmap, m = mpl_lib.norm_cmap(
        mn=mn,mx=mx,cm_name='brg')
    n=d.shape[0]

    ## find values in the back-bone branch shape
    ids = []; xs = []; ys = []; zs = []
    ax =_fig_ezz_.axes[0]
    axb=_fig_ezz_.axes[1]
    for i in xrange(n):
        x,y,z=d[i,0],d[i,1],-d[i,icol]
        if not(np.isnan(z)):
            zs.append(d[i,icol])
            ids.append(i)
            xs.append(d[i,0])
            ys.append(d[i,1])
        c = m.to_rgba(np.log10(z))
        ax.plot(x,y,'.',ms=1.5,
                color=c,mfc=c,mec='None')

    mpl_lib.add_cb(ax=axb,filled=True, format='%7.2e',
                   norm = norm, ylab=r'$\dot{\bar{E}}_{33}$')

    val, ind = sort(zs);
    imn = ids[ind[0]]
    ax.plot(d[imn,0],d[imn,1],'gx') ## minimum
    ax.set_xlim(-0.5,0.8)
    ax.set_ylim(-0.5,0.8)

    ## draw a backbone grid
    x0,y0=d[imn,0],d[imn,1] # center of backbone

    db=np.array(dist_bb); dr=np.array(dist_rb)
    DB=[]; DR=[]
    for i in xrange(len(db)): DB.append(-db[-i-1])
    DB.append(0)
    for i in xrange(len(db)): DB.append(db[i])
    for i in xrange(len(dr)): DR.append(-dr[-i-1])
    DR.append(0)
    for i in xrange(len(dr)): DR.append(dr[i])

    db=DB[::];dr=DR[::]
    xy_db_dr = np.zeros((len(db),len(dr),2))
    for i in xrange(len(db)):
        for j in xrange(len(dr)):
            xy_db_dr[i,j,:] = db[i], dr[j]

    ## rotate xy by fang
    for i in xrange(len(db)):
        for j in xrange(len(dr)):
            y,x = xy_db_dr[i,j,:]
            xy_db_dr[i,j,:] = rot_xy([x,y],fang)

    ## translate
    xy_db_dr[:,:,0] = xy_db_dr[:,:,0] + x0
    xy_db_dr[:,:,1] = xy_db_dr[:,:,1] + y0

    ## print xy_db_dr.shape
    ## return xy_db_dr

    if len(db)>3:
        nw = len(db)/3; nh = 4

    fig = plt.figure(figsize=(nh*3,nw*3))
    fig.suptitle(r'$\bar{E}_{33}$')
    fig1 = plt.figure(figsize=(nh*3,nw*3))
    fig1.suptitle(r'$\bar{E}_{33}$')
    gs  = gridspec.GridSpec(
        nw,nh,wspace=0.0,
        hspace=0.2,left=0.2,top=0.80)
    for i in xrange(nw*nh): fig.add_subplot(gs[i])
    gs  = gridspec.GridSpec(
        nw,nh,wspace=0.0,
        hspace=0.2,left=0.2,top=0.80)
    for i in xrange(4): fig1.add_subplot(gs[i])

    for _ax_ in fig.axes: _ax_.locator_params(nbins=4)
    for _ax_ in fig1.axes: _ax_.locator_params(nbins=4)

    iax=0; sym=['k^','k+','kd','kx','k<','k>','k*',
                'k1','k2','k3','k4','k8','ks','kp',
                'kh','kH','kD']

    ## thickness strain rate
    epst_dot_av = []; epst_dot_std = []
    exx_av      = []; exx_std      = []
    eyy_av      = []; eyy_std      = []
    ezz_av      = []; ezz_std      = []

    fout=open('xy_indx_ribs.txt','w')

    col = ('ib','ir','x[mm]','y[mm]','xyind')
    fout.write('%7s %7s %7s %7s %7s'%col); fout.write('\n')
    form = '%7i %7i %7.3f %7.3f %7i\n'
    form_nan = '%7i %7i %7s %7s %7s\n'
    t0 = time.time()
    for i in xrange(len(db)): ## bone
        indices_along_ribs = []; coord = []
        for j in xrange(len(dr)):   ## ribs
            x,y = xy_db_dr[i,j,:]   ## Coordinates of the necked spot
            ind = find_nearest([x,y],xs,ys,err=1.5)

            if not(np.isnan(ind)):
                indices_along_ribs.append(ids[ind])
                _x_ = d[ind,0]; _y_ = d[ind,1] ## actual coordinates
                coord.append([_x_,_y_]) ## actual coordinates
                fout.write(form%(i,j,_x_,_y_,ind))
            else: fout.write(form_nan%(i,j,'nan','nan','nan'))

        uet(time.time()-t0,head='Time spent for finding indices for BR')

        epsx = []; epsy = []; epsz = []; epst_dot = [];
        for k in xrange(len(indices_along_ribs)):
            indx = indices_along_ribs[k]
            ## ax.plot(x,y,sym[i],ms=3)
            ex=dat[:,indx,3]; ey=dat[:,indx,4]
            ## et=dat[:,indx,12] ## et=ez or ez_dot
            ez = -ex-ey

            ex_d=dat[:,indx,9]; ey_d=dat[:,indx,10]
            ez_d=ex_d+ey_d ## positive value

            epsx.append(ex); epsy.append(ey); epsz.append(ez)

            epst_dot.append(ez_d)
            fig.axes[iax].plot(ez_d) ##

        iax = iax + 1
        x, y = np.array(coord).T[:]
        ## Locations of ribs
        ax.plot(x,y,'-',color='gray',alpha=0.5)

        epst_dot=np.array(epst_dot).T
        epsx    =np.array(epsx).T
        epsy    =np.array(epsy).T
        epsz    =np.array(epsz).T

        D    = []; D_e   = [];E_ex = []; E_exe = [];
        E_ey = []; E_eye = [];E_ez = []; E_eze = [];
        for j in xrange(len(epst_dot)): ## number of ribs
            D.append(np.mean(epst_dot[j]))
            D_e.append(np.std(epst_dot[j]))
            E_ex.append(np.mean(epsx[j]))
            E_exe.append(np.std(epsx[j]))
            E_ey.append(np.mean(epsy[j]))
            E_eye.append(np.std(epsy[j]))
            E_ez.append(np.mean(epsz[j]))
            E_eze.append(np.std(epsz[j]))

        epst_dot_av.append(D)
        epst_dot_std.append(D_e)
        exx_av.append(E_ex)
        exx_std.append(E_exe)
        eyy_av.append(E_ey)
        eyy_std.append(E_eye)
        ezz_av.append(E_ez)
        ezz_std.append(E_eze)

    print
    fout.close()

    ## find maximum D value's index to trim the data
    ndat = (len(epst_dot_av)-1)/2+1
    ref=epst_dot_av[ndat-1] ## data at necking center spot
    _ref_=[]; _ind_ =[]
    for i in xrange(len(ref)):
        if not(np.isnan(ref[i])):
            _ref_.append(ref[i]); _ind_.append(i)
    _val_, _i_ = sort(_ref_)
    ix = _ind_[_i_[-1]]

    mdat = np.zeros((2,2,ndat)); c    = []
    for i in xrange(ndat):
        i0 =  i; i1 = -i -1
        ## thickness strain rate
        d1=epst_dot_av[i0];  d2=epst_dot_av[i1]
        e1=epst_dot_std[i0]; e2=epst_dot_std[i1]
        d =[]; e=[]

        ## epsilon xx
        ex1  = exx_av[i0];  ex2  = exx_av[i1]
        exe1 = exx_std[i0]; exe2 = exx_std[i1]
        ex=[]; exe=[]

        ## epsilon yy
        ey1  = eyy_av[i0];  ey2  = eyy_av[i1]
        eye1 = eyy_std[i0]; eye2 = eyy_std[i1]
        ey=[]; eye=[]

        ## epsilon zz
        ez1  = ezz_av[i0];  ez2  = ezz_av[i1]
        eze1 = ezz_std[i0]; eze2 = ezz_std[i1]
        ez=[]; eze=[]

        for j in xrange(len(d1)):  # along the time axis
            d.append((d1[j]+d2[j])/2.)
            e.append((e1[j]+e2[j])/2.)
            ex.append((ex1[j]+ex2[j])/2.)
            exe.append((exe1[j]+exe2[j])/2.)
            ey.append((ey1[j]+ey2[j])/2.)
            eye.append((eye1[j]+eye2[j])/2.)
            ez.append((ez1[j]+ez2[j])/2.)
            eze.append((eze1[j]+eze2[j])/2.)

        xind=np.arange(len(d))+1
        p = 0.7
        i0 = int(float(len(xind))*p)
        ## find maximum value and index?

        l, =fig1.axes[0].plot(xind,d,'-',label=abs(db[-i-1]))

        fig1.axes[1].plot(xind[i0:ix+1],d[i0:ix+1],'+')
        fig1.axes[2].errorbar(xind[i0:ix+1],d[i0:ix+1],
                              yerr=e1[i0:ix+1],fmt='+')

        l, =fig1.axes[3].plot(ex[i0:ix+1],ey[i0:ix+1])
        fig1.axes[3].errorbar(
            ex[ix],ey[ix],color=l.get_color(),
            xerr=exe[ix],yerr=eye[ix])

        mdat[0,0,i] = ex[ix]; mdat[0,1,i] = ey[ix]
        mdat[1,0,i] = exe[ix]; mdat[1,1,i] = eye[ix]
        c.append(l.get_color())

    ## deco
    from MP.lib.mpl_lib import tune_xy_lim
    fig1.axes[2].set_yscale('log')
    # fig1.axes[2].set_ylim(1e-3,)
    # tune_xy_lim(fig.axes);
    tune_xy_lim(fig1.axes)
    ral(fig.axes[1:]); ral(fig1.axes[1:])
    fig1.axes[0].legend(loc='best',ncol=2,fontsize=5,
                        bbox_to_anchor=(0.5,1.2)).\
        get_frame().set_alpha(0.5)
    return fig,fig1,mdat,c

def rot_xy(xy,th):
    """
    Arguments
    =========
    xy
    th
    """
    c=np.cos(th*np.pi/180.)
    s=np.sin(th*np.pi/180.)
    rot = np.array([[c,-s],[s,c]])
    return np.dot(rot,xy)

def find_nearest(xy0,xs,ys,err=None):
    """
    Find index of the nearest (xy0[0],xy0[1])
    coordinates of the given sets of coordinates
    (xs[..], ys[..])

    Arguments
    =========
    xy0
    xs
    ys
    err = None
    """
    ds=calc_dist(xy0,xs,ys)
    val,ind=sort(ds)
    if err!=None:
        if val[0]>err:
            return np.nan
    return ind[0]

def calc_dist(xy0,xs,ys):
    ds=[]
    for i in xrange(len(xs)):
        xy1=[xs[i],ys[i]]
        ds.append(__calc_dist__(xy0,xy1))
    return ds

def __calc_dist__(xy0,xy1):
    x0,y0=xy0; x1,y1=xy1
    return np.sqrt((x0-x1)**2 + (y0-y1)**2)

def ex02(
        ## non-smooth
        fns=[#'Bsteel_non_smooth/BB/2011NOV10.csv',
             'Bsteel_non_smooth/BB/2012JUL11.csv',
            #'Bsteel_non_smooth/PSRD/2011NOV17.csv',
             'Bsteel_non_smooth/PSRD/2012JUL12.csv',
             'Bsteel_non_smooth/PSRD/2012JUL13.csv',
            #'Bsteel_non_smooth/PSTD/2011NOV17.csv',
             'Bsteel_non_smooth/PSTD/2012JUL17.csv'],
        fangs=[#90,
               90,
            #0,
               90,
               90,
            #90,
               0],
        dist_bb=[2,4,6],
        dist_rb=[2,3,4,5,6,7,8],

        ## dist=[5,10,15,20],
        ixy_flip=[#False,
                  False,
            #False,
            True,
            True,
            #False,
            False],
        ntot_col=[#14,
                  14,
            #14,
                  14,
                  14,
            #14,
                  14],
        iopt=0,
        icol_ezz = 12, ## total ezz strain
        ioff=True):
    """
    fn       = file name: " /dum/[strainpath]/[date].csv "

    fang     = fracture angle

    ioff     = interactive mode switch for MPL
    ntot_col = total number of column belonging to a dat block

    dist_bb     = distances from the critical point
    dist_rb     = distances from the critical point

    ixy_flip = flag if x-y cooridnate should be switched
    ntot_col = [14, ..] # number of variables aligned for columns


    ## below is the column information for a block
    icol:     0 1 2 3   4   5   6  7  8  9      10     11     12
    variable: x,y,z,exx,eyy,exy,e1,e2,vm,exxdot,eyydot,exydot,e33

    ## Or below is the column information for a block
    icol:     0 1 2 3   4   5   6  7    8      9     10   11  12
    variable: x,y,z,exx,eyy,exy,e1,e2,exxdot,eyydot,exydot,e33,e33_dot
    """

    from MP.lib.mpl_lib import wide_fig as wf
    import matplotlib.pyplot as plt
    from MP.lib.axes_label import __deco_fld__ as deco_fld

    fn = open('EXP_FLD.txt','w')
    fig4=wf(nw=1,nh=1)
    data=[]
    for i in xrange(len(fns)):
        #try:
            ## print fns[i]
        d, fig1, fig2, fig3 = main_disc(
            fn=fns[i],
            dist_bb=dist_bb,
            dist_rb=dist_rb,
            fang=fangs[i],
            ixy_flip=ixy_flip[i],
            ntot_col=ntot_col[i],
            icol_ezz=icol_ezz,iopt=iopt,
            ioff=ioff,fig=fig4,fout=fn)
        data.append(d)
        # except ValueError:
        #     print 'error occured in %s'%fns[i]
    print 'Analysis results have been saved to %s'%fn.name
    fn.close()
    fig1.savefig('fig1.pdf')
    fig2.savefig('fig2.pdf')
    fig3.savefig('fig3.pdf')
    fig4.axes[0].set_xlim(-0.5,0.5)
    fig4.axes[0].set_ylim(-0.5,0.5)
    fig4.savefig('EXP_FLD.pdf')
    return data
