### read DIC profile csv file and carry out analysis
import numpy as np
import matplotlib.pyplot as plt
def read(fn='x_vs_exx_profile.csv'):
    ## dat=np.loadtxt(fn,skiprows=1,dtype='float').T
    f = open(fn,'r')
    lines = f.readlines()
    dats = []
    for i in range(len(lines)-1):
        d=map(float,lines[i+1].split('\n')[0].split(',')[:-1])
        dats.append(d)

    dats=np.array(dats).T
    x=dats[::2]
    ys=dats[1::2]
    x=x[0]

    return x, ys

def find_lost_correlation(ys,nlast=15,offset=0.1,ndat=4):
    ## find the one that starts loosing correlation...
    i0=None
    ysi = ys[::-1][:nlast]
    ysi = ysi[::-1]
    for i in range(nlast):
        d=ysi[i]
        d_mean=np.mean(d)

        ## print 'd_mean value:', d_mean
        for j in range(len(d)):
            if d[j]<offset*d_mean:
                ## print d[j]
                if i0==None:
                    i0=nlast-i
    ind0=-i0+1-ndat
    ind1=-i0+1
    return ind0,ind1

def plot(fn='x_vs_exx_profile.csv',x=None,ys=None,ifig=1):
    xy_dist=fn.split('_')[0].upper()
    eps_xy=fn.split('_')[2][1:]
    from MP.lib.mpl_lib import wide_fig as wf
    fig=wf(nw=2,nh=3,iarange=True,w0=0.2,ws=0.6,w1=0.2,
           uw=3.0,uh=3.0)
    ax1=fig.axes[0];ax2=fig.axes[1];ax3=fig.axes[2];
    ax4=fig.axes[3];ax5=fig.axes[4];ax6=fig.axes[5];
    ax2.set_yscale('log');ax4.set_yscale('log');ax6.set_yscale('log')

    xylab(ax1,iopt=0,xy=eps_xy,xyp=xy_dist)
    xylab(ax3,iopt=0,xy=eps_xy,xyp=xy_dist)
    xylab(ax5,iopt=0,xy=eps_xy,xyp=xy_dist)
    xylab(ax2,iopt=1,xy=eps_xy,xyp=xy_dist)
    xylab(ax4,iopt=1,xy=eps_xy,xyp=xy_dist)
    xylab(ax6,iopt=1,xy=eps_xy,xyp=xy_dist)

    ## x,ys=read(fn=fn)
    for i in range(len(ys)):
        ax1.plot(x,ys[i],'x')

    nlast=10
    syms=['o','x','+','d','^','*','>','<','s','1']

    for i in range(nlast):
        ax3.plot(x,ys[::-1][i],syms[i],label=i+1)

    yst=ys.T; ystd=np.diff(yst); ysd = ystd.T
    for i in range(len(ysd)):
        ax2.plot(x,ysd[i],'x-')

    for i in range(nlast):
        try:
            ax4.plot(x,ysd[::-1][i],'%s-'%syms[i],label=i+1)
        except ValueError: pass

    ax3.legend(loc='best',ncol=2,fontsize=7)

    # ## find the one that starts loosing correlation...
    ndat=4
    ind0,ind1 = find_lost_correlation(ys,nlast=15,offset=0.1,ndat=ndat)
    ##

    ysi=ys[ind0:ind1]; ysdi=ysd[ind0:ind1]
    ##
    ind_lost0=None;ind_lost1=None
    for i in range(len(ysi[-1])):
        if ysi[-1][i]==0: ## replace lost correlation datum to np.nan
            ysi[-1][i]=np.nan
            if ind_lost0==None:
                ind_lost0=i-1

    for i in range(len(ysi[-1][::-1])):
        if np.isnan(ysi[-1][::-1][i]):
            ind_lost1=len(ysi[-1])-i
            break

    for i in range(5):
        ax5.plot(x[ind_lost0-i],ysi[-1][ind_lost0-i],'kx')
    #ax5.plot(x[ind_lost1],ysi[-1][ind_lost1],'kx')

    print ind_lost0, ind_lost1

    for i in range(ndat):
        ax5.plot(x,ysi[i])
        ax6.plot(x,ysdi[i])

    return fig,ind_lost0, ind_lost1

def xylab(ax,iopt=0,ft=12,xy='xx',xyp='X'):
    ax.set_xlabel('%s Position [mm]'%xyp,dict(fontsize=ft))
    if iopt==0:
        ax.set_ylabel(r'$\bar{E}_{\mathrm{%s}}$'%xy,
                      dict(fontsize=ft))
    if iopt==1:
        ax.set_ylabel(r'$\dot{\bar{E}}_{\mathrm{%s}}$'%xy,
                      dict(fontsize=ft))
    if iopt==2:
        ax.set_ylabel(r'$|\dot{\bar{E}}_{\mathrm{%s}}|$'%xy,
                      dict(fontsize=ft))

def main_old(fns=['x_vs_exx_profile.csv','x_vs_eyy_profile.csv']):
    from MP.lib.mpl_lib import wide_fig as wf
    from MP.lib.axes_label import __deco_fld__ as deco_fld
    dat=[]
    for i in range(len(fns)):
        x,ys=read(fns[i])
        dat.append(ys)
        fig, i0,i1=plot(fns[i],x=x,ys=ys,ifig=i)
    exx,eyy = dat;exx=exx.T;eyy=eyy.T
    fig=wf(nw=1,nh=1,iarange=True,w0=0.2,ws=0.6,w1=0.2,
           uw=3.0,uh=3.0)
    #return dat
    for i in range(20):
        fig.axes[0].plot(eyy[i0-i],exx[i0-i],'k-',alpha=0.5)
        fig.axes[0].plot(eyy[i0-i][-1],exx[i0-i][-1],'rx')

    deco_fld(fig.axes[0],iopt=2)
    fig.axes[0].set_xlim(-0.3,0.5)
    fig.axes[0].set_ylim(0.,0.8)

### data analyzed by and extracted from a newer Vic3D
def read_lines(fn='BB/2012JUL11/epsilon_xlines.csv',
               xy_trans=False):
    from StringIO import StringIO
    lines=open(fn,'r').read()
    lines=lines.split('\n')
    nline_per_block=None
    i=0
    while nline_per_block==None:
        if len(lines[i])==0:
            nline_per_block=i
        else:i=i+1

    ## print lines[1]
    header = lines[1].split(',')
    n_probe_line=0
    for i in range(len(header)):
        if len(header[i])>0:
            n_probe_line=n_probe_line+1

    nblock = (len(lines)-1)/nline_per_block - 1
    npositions = nline_per_block - 3

    # print '# blocks:', nblock
    # print '# nline_per_block:', nline_per_block
    # print '# positions :', npositions
    # print '# lines :', n_probe_line

    i0=0
    i1=nline_per_block

    master=np.zeros((3,n_probe_line,npositions,nblock))
    for i in range(nblock):             ## per a frame
        LINE = lines[i0:i1]
        LINE = LINE[3:]
        for j in range(len(LINE)):      ## per position point
            dat=LINE[j].split(',')
            ## print len(dat)
            for k in range(len(dat)/3): ## per line profile
                try:
                    x  =float(dat[k*3])
                    exx=float(dat[k*3+1])
                    eyy=float(dat[k*3+2])
                except:
                    x = np.nan
                    exx=np.nan
                    eyy=np.nan

                master[0,k,j,i] = x
                master[1,k,j,i] = exx
                master[2,k,j,i] = eyy

        i0 = i0 + nline_per_block+1
        i1 = i1 + nline_per_block+1

    if xy_trans:
        dum=master[1,:,:,:].copy()
        master[1,:,:,:] = master[2,:,:,:][::]
        master[2,:,:,:] = dum[::]

    return master

def find_lost_of_correlation(ydat):
    x=ydat[1,:,:]
    xt=x.T
    i0 = None; j0 = None
    for i in range(len(xt)):
        for j in range(len(xt[i])):
            if np.isnan(xt[i][j]):
                i0 = i
                j0 = j
                return i0-1, j0-1
    return i0,j0

def read_disc(fn='BB/2012JUL11/node_data.csv',ndat_col=13,
              ixy_flip=False):
    dat=np.genfromtxt(fn,delimiter=',',skiprows=1)
    shp = dat.shape
    print shp
    print ndat_col
    print float(shp[1])/float(ndat_col)
    ## raise IOError
    dat=np.reshape(dat, (shp[0],shp[1]/ndat_col,ndat_col))


    if ixy_flip:
        dum = dat[:,:,0][::].copy() ## x coordinate
        dat[:,:,0] = dat[:,:,1][::] .copy() ## y cooridnate
        dat[:,:,1] = dum[::].copy()
    return dat

def main(fn='BB/2012JUL11/epsilon_xlines.csv',
         n_adjacent = 1,n_last_point = 1,xy_trans=False):
    dat_string = ''

    from MP.lib.axes_label import __deco_fld__ as deco_fld
    from MP.lib.mpl_lib import wide_fig as wf
    fig=wf(nw=3,nh=2,iarange=True,w0=0.2,ws=0.6,w1=0.2,
           uw=3.0,uh=3.0)

    dat=read_lines(fn,xy_trans=xy_trans)
    n_lines    =dat.shape[1]
    n_positions=dat.shape[2]
    n_frames   =dat.shape[3]

    new_dat=np.zeros((7,n_lines,n_positions,n_frames))
    for i in range(3):
        new_dat[i,:,:,:] = dat[i,:,:,:]
    new_dat[3,:,:,:] = +dat[1,:,:,:]+dat[2,:,:,:]
    new_dat[3,:,:,:] = new_dat[3,:,:,:] * (-1)

    exx_diff=np.diff(new_dat[1,:,:,:])
    eyy_diff=np.diff(new_dat[2,:,:,:])
    ezz_diff=np.diff(new_dat[3,:,:,:])

    new_dat[:,:,:,0] = np.nan
    new_dat[4,:,:,1:]=exx_diff[::]
    new_dat[5,:,:,1:]=eyy_diff[::]
    new_dat[6,:,:,1:]=ezz_diff[::]

    ## fig_lines=wf(nw=2,nh=3,iarange=True)
    for i in range(n_lines):
        c = []
        dat_line=dat[:,i,:,:]
        i0, i1 = find_lost_of_correlation(dat_line)
        if i0==None or i1==None: pass
        else:
            x   = dat[0,i,:,i0]
            exx = dat[1,i,:,i0] ## exx
            eyy = dat[2,i,:,i0] ## eyy
            l, = fig.axes[0].plot(
                x,exx,label='line #%i'%(i))
            c.append(l.get_color())
            fig.axes[1].plot(x,eyy,color=l.get_color())

            fig.axes[0].plot(x[i1],exx[i1],'x',
                             color=l.get_color(),
                             label='line #%i'%(i))

            fig.axes[1].plot(x[i1],eyy[i1],'x',
                             color=l.get_color(),
                             label='line #%i'%(i))

            for j in range(n_adjacent):
                x   = dat[0,i,i1-j,0:i0]
                exx = dat[1,i,i1-j,0:i0]
                eyy = dat[2,i,i1-j,0:i0]
                fig.axes[2].plot(
                    eyy,exx,
                    color=l.get_color())
                # fig.axes[3].plot(
                #     eyy[-1],exx[-1],'x',
                #     color=l.get_color())

                for k in range(n_last_point):
                    x      = new_dat[0,i,i1-j,i0-k]
                    exx    = new_dat[1,i,i1-j,i0-k]
                    eyy    = new_dat[2,i,i1-j,i0-k]
                    ezz    = new_dat[3,i,i1-j,i0-k]
                    exxdot = new_dat[4,i,i1-j,i0-k]
                    eyydot = new_dat[5,i,i1-j,i0-k]
                    ezzdot = new_dat[6,i,i1-j,i0-k]
                    dum= '%9i %9i %9i %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f %6s'%(
                        i,j,k,x,exx,eyy,ezz,
                        exxdot,eyydot,ezzdot,
                        l.get_color())
                    print dum
                    dat_string =dat_string + dum + '\n'

            for j in range(n_last_point):
                x      = new_dat[0,i,:,i0-j]
                exxdot = new_dat[4,i,:,i0-j]
                eyydot = new_dat[5,i,:,i0-j]
                ezzdot = new_dat[6,i,:,i0-j]

                fig.axes[3].plot(
                    x,exxdot)
                fig.axes[4].plot(
                    x,eyydot,l.get_color())
                fig.axes[5].plot(
                    x,abs(ezzdot),l.get_color())

    fig.axes[3].set_yscale('log')
    fig.axes[4].set_yscale('log')
    fig.axes[5].set_yscale('log')

    xylab(fig.axes[0],iopt=0,xy='xx',xyp='X')
    xylab(fig.axes[1],iopt=0,xy='yy',xyp='X')
    fig.axes[2].set_ylabel(r'$\bar{E}_\mathrm{RD}$')
    fig.axes[2].set_xlabel(r'$\bar{E}_{\mathrm{TD}}$')

    # deco_fld(fig_lines.axes[0],iopt=2)
    # deco_fld(fig_lines.axes[1],iopt=2)

    xylab(fig.axes[3],iopt=1,xy='xx',xyp='X')
    xylab(fig.axes[4],iopt=1,xy='yy',xyp='X')
    xylab(fig.axes[5],iopt=2,xy='zz',xyp='X')

    return dat_string

def main_disc(fn='BB/2011NOV10/node_disc_36.csv',
              dist=[5,10,15,20],
              fang=90.,ixy_flip=False,ioff=True,
              ntot_col=14):
    """
    fn       = file name: " /dum/[strainpath]/[date].csv "
    dist     = distances from the critical point
    fang     = fracture angle
    ixy_flip = flag if x-y cooridnate should be switched
    ioff     = interactive mode switch for MPL
    ntot_col = total number of column belonging to a dat block

    ## below is the column information for a block
    icol:     0 1 2 3   4   5   6  7  8  9      10     11     12
    variable: x,y,z,exx,eyy,exy,e1,e2,vm,exxdot,eyydot,exydot,e33
    """
    import os
    from MP.lib.axes_label import __deco_fld__ as deco_fld
    import matplotlib.pyplot as plt
    from MP.lib.mpl_lib import wide_fig as wf
    from MP.lib.mpl_lib import rm_all_lab as ral

    ## output file name
    strain_path = fn.split(os.sep)[1]
    date = fn.split(os.sep)[2].split('.')[0]
    fn_out = '%s_%s'%(strain_path,date)

    ## MPL interactive mode switch
    if ioff: plt.ioff()
    else: plt.ion()

    ## Read data from VIC-3d software
    dat=read_disc(fn=fn,ndat_col=ntot_col,
                  ixy_flip=ixy_flip)

    ## through thickness strain rate minimum
    fig=wf(nw=3,nh=4,w0=0.0,ws=1.0,w1=0.0,
           h0=0.0,hs=1.0,h1=0.0,uw=2.3,uh=2.3,
           left=0.2,up=0.1,down=0.2,iarange=True)
    fig1=wf(nw=5,nh=1,left=0.1,w0=0.2,ws=0.6,w1=0.2,
            down=0.15,iarange=True)
    n_image=len(fig.axes)
    IND=None; dum=None


    for i in range(len(fig.axes)):
        ax=fig.axes[i]
        plot_xy_disc_e33dot(dat,istp=-n_image-i,icol=12,
                            ax=ax,
                            fang=fang,
                            dist_bb=dist,
                            dist_rib=[1,2,3,4])
        return
                        
                        


    # """ Core of the 'point' method """
    # ## Strain color-map for the last 12 images
    # ## Also, find coordinate
    # for i in range(len(fig.axes)):
    #     ax=fig.axes[i]
    #     dum = plot_xy_2dc(
    #         dat,istp=-n_image+i,ax=ax,icol=12,dist=dist,
    #         fang=fang,is_max=False,IND=dum)
    #     if i==n_image-1: IND=dum[::]

    # dst=[0]
    # for i in range(len(dist)):
    #     dst.append(dist[i]); dst.append(-dist[i])

    # ls = ['k-','r-','r--','b-','b--','g-','g--',
    #       'm-','m--','c-','c--']
    # EYY,EXX,LS=[],[],[]
    # for i in range(len(dst)):
    #     ind = IND[i]
    #     exx=dat[:,ind,3]
    #     eyy=dat[:,ind,4]
    #     e33=dat[:,ind,12]
    #     exxdot=dat[:,ind,9]
    #     eyydot=dat[:,ind,10]
    #     exydot=dat[:,ind,11]
    #     ezzdot = - exxdot - eyydot
    #     ndat=len(eyy); index = np.arange(ndat)
    #     ndat0 = int(ndat*0.9)
    #     fig1.axes[0].plot(index,-e33,            ls[i],label=dst[i])
    #     fig1.axes[1].plot(index[ndat0:],-e33[ndat0:],   ls[i],label=dst[i])
    #     fig1.axes[2].plot(index,-ezzdot,         ls[i],label=dst[i])
    #     fig1.axes[3].plot(index[ndat0:],-ezzdot[ndat0:],ls[i],label=dst[i])

    #     fig1.axes[4].plot(eyy,exx,ls[i])## RD//y TD//x
    #     EYY.append(eyy)
    #     EXX.append(exx)
    #     LS.append(ls[i])
    # """ Core of the 'point' method """

    # ## Deco figures
    # fig1.axes[0].legend(loc='best',ncol=2,fontsize=7,fancybox=True).\
    #     get_frame().set_alpha(0.5)

    # ## deco_fld(fig1.axes[2],iopt=2)
    # for i in range(4):
    #     fig1.axes[i].set_xlabel('index')

    # fig1.axes[0].set_ylabel(r'$\bar{E}_{33}$')
    # fig1.axes[1].set_ylabel(r'$\bar{E}_{33}$')
    # fig1.axes[2].set_ylabel(r'$\dot{\bar{E}}_{33}$')
    # fig1.axes[3].set_ylabel(r'$\dot{\bar{E}}_{33}$')

    # fig1.axes[2].set_yscale('log')
    # fig1.axes[3].set_ylim(0.,)

    # fig.axes[0].set_title(r'$\dot{\bar{E}}_{33}$')
    # ral(fig.axes)
    # fig.savefig('%s_e33dot.png'%fn_out)
    # fig1.savefig('%s_e33dot_e33_time1.png'%fn_out)

    # return EYY,EXX,LS

def find_coordinate_index(
    dat_disc,ind_image=-1,xy=[0,0],is_max=True,
    icol=11,dist_filter=2.):
    """
    Find index that cooresponds
    to the given 'xy' coordinate
    """
    from MP import ssort
    sort=ssort.shellSort

    x0,y0 = xy
    d     = dat[istp,:,:]
    x,y   = dat[istp,:,0],dat[istp,:,1]
    dist  = []; ind_ref = []
    for i in range(len(x)):
        if np.isnan(x0) or np.isnan(y0): pass
        else:
            d=np.sqrt((x0-x)**2 + (y0-y)**2)
            if d<dist_filter:
                dist.append(d)
                ind_ref[i]
    if len(ind_ref)==0: return -1
    dist_sorted, indc = sort(dist)
    return in_ref[indc[0]]

def plot_xy_disc_e33dot(dat,istp=-1,
                        ax=None,
                        icol=11,
                        fang=90,
                        dist_bb=[2.5,5,10], ## backbone
                        dist_rib=[1,2,3,4]
                        ):
                        
    """
    Extract DIC data from 'disc' area
    """
    from MP.lib.mpl_lib import wide_fig as wf
    from MP.lib.mpl_lib import rm_all_lab as ral

    import matplotlib as mpl
    from MP.lib import mpl_lib
    from MP import ssort


    sort=ssort.shellSort

    d=dat[istp,:,:]
    dz = np.log10(abs(d[:,icol]))
    for i in range(len(dz)):
        if dz[i]==np.inf or dz[i]==-np.inf:
            dz[i] = np.nan
    mn=min(dz);mx=max(dz)
    norm=mpl.colors.Normalize(vmin=mn, vmax=mx)
    cmap, m = mpl_lib.norm_cmap(
        mn=mn,mx=mx,cm_name='brg')
    n=d.shape[0]


    ## find values in the back-bone branch shape
    ids = []; xs = []; ys = []; zs = []
    for i in range(n):
        x,y,z=d[i,0],d[i,1],-d[i,icol]
        if np.isnan(z): pass
        else:
            zs.append(d[i,icol])
            ids.append(i)
            xs.append(d[i,0])
            ys.append(d[i,1])
        c = m.to_rgba(np.log10(z))
        ax.plot(x,y,'o',ms=1.5,
                color=c,mfc=c,mec='None')

    val, ind = sort(zs)
    imn = ids[ind[0]]
    ax.plot(d[imn,0],d[imn,1],'g.') ## minimum

    ## draw a backbone grid
    x0,y0=d[imn,0],d[imn,1] # center of backbone

    db=np.array(dist_bb); dr=np.array(dist_rib)
    DB=[]; DR=[]
    for i in range(len(db)):
        DB.append(-db[-i-1])
    DB.append(0)
    for i in range(len(db)):
        DB.append(db[i])
    for i in range(len(dr)):
        DR.append(-dr[-i-1])
    DR.append(0)
    for i in range(len(dr)):
        DR.append(dr[i])

    db=DB[::];dr=DR[::]
    xy_db_dr = np.zeros((len(db),len(dr),2))
    for i in range(len(db)):
        for j in range(len(dr)):
            xy_db_dr[i,j,:] = db[i], dr[j]

    ## rotate xy by fang
    for i in range(len(db)):
        for j in range(len(dr)):
            y,x = xy_db_dr[i,j,:]
            xy_db_dr[i,j,:] = rot_xy([x,y],fang)

    ## translate
    xy_db_dr[:,:,0] =     xy_db_dr[:,:,0] + x0 
    xy_db_dr[:,:,1] =     xy_db_dr[:,:,1] + y0

    fig=wf(nw=len(db),nh=len(dr),w0=0.0,ws=1.0,w1=0.0,
           h0=0.0,hs=1.0,h1=0.0,uw=2.3,uh=2.3,
           left=0.2,up=0.1,down=0.2,iarange=True)
    fig1=wf(nw=1,nh=1,w0=0.0,ws=1.0,w1=0.0,
            h0=0.0,hs=1.0,h1=0.0,uw=2.3,uh=2.3,
            left=0.2,up=0.1,down=0.2,iarange=True)

    iax=0
    sym=['k^','k+','kd','kx','k<','k>','k*',
         '','']
    for i in range(len(db)):
        indices_along_ribs = []
        for j in range(len(dr)):
            x,y = xy_db_dr[i,j,:]
            ind = find_nearest([x,y], xs,ys)
            indices_along_ribs.append(ind)

        epsx = []; epsy = []; epst_dot = []
        for k in range(len(indices_along_ribs)):
            ind = indices_along_ribs[k]
            x = d[ind,0]
            y = d[ind,1]
            ax.plot(x,y,sym[i],ms=3)

            ex=dat[:,ind,9]
            ey=dat[:,ind,10]
            et=dat[:,ind,12]

            ex_d=dat[:,ind,9]
            ey_d=dat[:,ind,10]
            ez_d=-ex_d-ey_d
            epst_dot.append(ez_d)
            fig.axes[iax].plot(ez_d)
            iax=iax + 1

        epst_dot=np.array(epst_dot).T
        D=[]
        for k in range(len(epst_dot)):
            D.append(np.mean(epst_dot[k]))
        fig1.axes[0].plot(D,label=db[i])
    fig1.axes[0].legend(loc='best',ncol=2,fontsize=7)

def plot_xy_2dc(dat,istp=-1,ax=None,icol=11,
                is_max=True,fang=90,
                dist=[2.5,5,10],IND=None):
    """
    Extract DIC data from 'disc' area
    """
    import matplotlib as mpl
    from MP.lib import mpl_lib
    from MP import ssort
    sort=ssort.shellSort

    d=dat[istp,:,:]
    dz = np.log10(abs(d[:,icol]))
    for i in range(len(dz)):
        if dz[i]==np.inf or dz[i]==-np.inf:
            dz[i] = np.nan
    mn=min(dz);mx=max(dz)

    norm=mpl.colors.Normalize(vmin=mn, vmax=mx)
    cmap, m = mpl_lib.norm_cmap(
        mn=mn,mx=mx,cm_name='brg')
    n=d.shape[0]

    ids = []; xs = []; ys = []; zs = []
    ## index seeking if necessary
    for i in range(n):
        x=d[i,0]
        y=d[i,1]
        z=abs(d[i,icol])
        if np.isnan(z): pass
        else:
            zs.append(d[i,icol])
            ids.append(i)
            xs.append(d[i,0])
            ys.append(d[i,1])
        c = m.to_rgba(np.log10(z))
        ax.plot(x,y,'o',ms=1.5,
                color=c,mfc=c,mec='None')
    if IND==None:
        val, ind = sort(zs)
        ## max
        imx=ids[ind[-1]]
        ax.plot(d[imx,0], d[imx,1],'r.')
        ## min
        imn=ids[ind[0]]
        ax.plot(d[imn,0], d[imn,1],'g.')

        ## draw a line [-10, 10] perpendicular to the min
        ## find locations close to [-10, -5, 0, 5, 10]
        IND = np.zeros((len(dist)*2+1),dtype='int')
        if is_max:
            x=d[imx,0];y=d[imx,1]
            IND[0] = imx
        else:
            x=d[imn,0];y=d[imn,1]
            IND[0] = imn
        for i in range(len(dist)):
            xy0, xy1 = find_xy(
                r=dist[i],fang=fang,
                x=x,y=y)
        ## find coordinate index that is close to those numbers.
        ## calculate distance from this point
            ind0 = find_nearest(xy0,xs,ys)
            ind1 = find_nearest(xy1,xs,ys)
            IND[i*2+1]=ind0; IND[i*2+2]=ind1
            ax.plot(xs[ind0],ys[ind0],'kx')
            ax.plot(xs[ind1],ys[ind1],'k+')
        return IND
    else:
        for i in range(len(IND)):
            ind0 = IND[i]
            if i==0:
                sym='g.'
            elif np.mod(i,2)==1:
                sym='kx'
            elif np.mod(i,2)==0:
                sym='k+'
            ax.plot(d[ind0,0],d[ind0,1],sym)
            #ax.plot(d[ind0],ys[ind0],sym)
        return IND

def rot_xy_o(xy,o,th):
    xy=np.array(xy)
    o =np.array(o)
    XY = xy - o
    XY = rot_xy(XY,th)
    return XY + o
    

def rot_xy(xy,th):
    c=np.cos(th*np.pi/180.)
    s=np.sin(th*np.pi/180.)
    rot = np.array([[c,-s],[s,c]])
    return np.dot(rot,xy)

def find_xy(r,fang,x,y):
    fang=np.pi*fang/180. + np.pi/2.
    x1=x+r*np.cos(fang);    y1=y+r*np.sin(fang)
    x2=x-r*np.cos(fang);    y2=y-r*np.sin(fang)
    return [x1,y1],[x2,y2]

def find_nearest(xy0,xs,ys):
    from MP import ssort
    sort=ssort.shellSort
    ds=calc_dist(xy0,xs,ys)
    val,ind=sort(ds)
    return ind[0]

def calc_dist(xy0,xs,ys):
    ds=[]
    for i in range(len(xs)):
        xy1=[xs[i],ys[i]]
        ds.append(__calc_dist__(xy0,xy1))
    return ds

def __calc_dist__(xy0,xy1):
    x0,y0=xy0; x1,y1=xy1
    return np.sqrt((x0-x1)**2 + (y0-y1)**2)

def ex02(
    ## non-smooth
    fns=['Bsteel_non_smooth/BB/2011NOV10.csv',
         'Bsteel_non_smooth/BB/2012JUL11.csv',
         'Bsteel_non_smooth/PSRD/2011NOV17.csv',
         'Bsteel_non_smooth/PSRD/2012JUL12.csv',
         'Bsteel_non_smooth/PSRD/2012JUL13.csv',
         'Bsteel_non_smooth/PSTD/2011NOV17.csv',
         'Bsteel_non_smooth/PSTD/2012JUL17.csv'],
    fangs=[90,90,0,0,0,90,90],
    dist=[3,5,7],
    ## dist=[5,10,15,20],
    ixy_flip=[False,False,False,False,
              False,False,True],
    ntot_col=[14,14,14,14,14,14,14],
    ioff=True):

    from MP.lib.mpl_lib import wide_fig as wf
    import matplotlib.pyplot as plt
    from MP.lib.axes_label import __deco_fld__ as deco_fld

    fig=wf(nw=2,nh=1)
    for i in range(len(fns)):
        try:
            EYY,EXX,LS = main_disc(fn=fns[i],
                              dist=dist,
                              fang=fangs[i],
                              ixy_flip=ixy_flip[i],
                              ntot_col=ntot_col[i],
                              ioff=ioff)
            for j in range(len(EYY)):
                fig.axes[0].plot(EYY[j],EXX[j],LS[j])

                ind=-1
                while (np.isnan(EXX[j][ind])
                       or
                       np.isnan(EYY[j][ind])):
                    ind=ind-1
                fig.axes[1].plot(EYY[j][ind],EXX[j][ind],
                                 'o',mfc='None',mec=LS[j][0])
        except ValueError:
            print 'error occured in %s'%fns[i]
    deco_fld(fig.axes[0],iopt=2)
    deco_fld(fig.axes[1],iopt=2)
    fig.axes[0].grid('on');fig.axes[1].grid('on')
    plt.draw()
    fig.savefig('EXP_FLD.pdf'); fig.savefig('EXP_FLD.png')

def ex01(fns=['BB/2012JUL11/epsilon_xlines.csv',
              'PSRD/2012JUL12/epsilon_ylines.csv',
              'PSTD/2012JUL17/epsilon.csv'],
         xy_trans=[False,True,False]):

    from MP.lib.mpl_lib import wide_fig as wf
    from StringIO import StringIO
    print '%9s %9s %9s %9s %9s %9s %9s %9s %9s %9s %6s'%(
        'lines','adjacent',
        'tframes','Loc','exx','eyy','ezz','exxdot',
        'eyydot','ezzdot','color')
    dat=''
    for i in range(len(fns)):
        dstring = main(
            fn=fns[i],
            n_adjacent=3,
            n_last_point=2,
            xy_trans=xy_trans[i])
        print ' --'
        dat=dat+dstring

    dat=np.loadtxt(StringIO(dat),skiprows=1,dtype='string').T
    exx=map(float,dat[4])
    eyy=map(float,dat[5])
    fig=wf(nw=1,nh=1)
    for i in range(len(dat[4])):
        fig.axes[0].plot(dat[4][i],dat[5][i],'x',
                         color=dat[10][i])

    np.savetxt('FLD.out',np.array([exx,eyy]).T)

