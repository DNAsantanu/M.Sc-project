import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.colors import LogNorm
import matplotlib.ticker as ticker
from matplotlib.ticker import FormatStrFormatter
from matplotlib import rc

#### INPUTS ########################


#### MLE ###################

r=9.448626e+03	# Mpc
rprime=1.773767e+01	# Mpc/MHz
deltanu=0.1	# MHz

############################


#### Bin along kparallel ###

kmin=0.025
kmax=0.525
mode='lin'

############################


#### Plot Heatmap ##########

nuc=166.2 # MHz
Dia=23.850471219185987 # m

############################


####################################

#### Horizon #######

def f(x):
    a=r*x/(rprime*nuc)
    return a


#### First null #######

def g(x):
    b=(r/(rprime*nuc))*1.22*300.*x/(nuc*Dia)
    return b


def errcov(ell,N_c,Nrea,filename,Output):

    ######### Estimate the mean #######################################

    mean=np.zeros([N_c,ell])

    for ii in range (Nrea):

        real=np.load(filename+str(ii+1)+'_cl_sort.npy')
#		l=real[0,1:]
        real=real[1:,:]
        mean+=real[:N_c,1:ell+1]

    mean/=Nrea

    ######### Estimate the error covariance ###########################

    cov=np.zeros([N_c,N_c,ell])

    for ii in range (ell):

        for jj in range (Nrea):
            real=np.load(filename+str(jj+1)+'_cl_sort.npy')
            real=real[1:,:]
            diff=real[:N_c,ii+1]-mean[:,ii]

            for kk in range (N_c):
                for ll in range (N_c):
                    cov[kk,ll,ii]+=(diff[kk]*diff[ll])

    cov/=(Nrea-1)

    np.save(Output,cov)

    return cov


def MLE(cl,var,N_E,window_type,OUTNAME):

    ell=cl[0,1:]
    cl=cl[1:,:]
    cl[:,0]*=deltanu

#	deltanu=cl[1,0]/1.e6	# MHz
    N_c=len(cl)
    nbin_kper=len(cl[0])-1
    print("deltanu=%le N_c=%d nbin=%d"%(deltanu,N_c,nbin_kper))

    def window(typ):
        if(typ=='BN'):
            A=[3.6358e-1, 4.891775e-1, 1.3659e-1, 1.06411e-2, 0., 0., 0.]
        if(typ=='BH4'):
            A=[3.5875e-1, 4.8829e-1, 1.4128e-1, 1.1680e-2, 0., 0., 0.]
        if(typ=='MS5'):
            A=[3.2321e-1, 4.7149e-1, 1.7553e-1, 2.8496e-2, 1.2613e-3, 0., 0.]
        if(typ=='MS6'):
            A=[2.9355e-1, 4.5193e-1, 2.0141e-1, 4.7926e-2, 5.0261e-3, 1.3755e-4, 0.]
        if(typ=='MS7'):
            A=[2.7122e-1, 4.3344e-1, 2.1800e-1, 6.5785e-2, 1.07618e-2, 7.7001e-4, 1.3680e-5]
        w=np.zeros(N_c,dtype='double')
        for ii in range (0,len(w),1):
            if(typ=='none'):
                w[ii]=1.
            else:
                for jj in range (0,len(A),1):
                    w[ii]+=((-1.)**jj*A[jj]*np.cos(2.*(ii+N_c-1.)*jj*np.pi/(((2.*N_c)-1)-1.)))
        return(w)


    win=window(window_type)
#	print(win)

    kper=ell/r

    kpara=[]
    power=[]

    for jj in range (0,N_E,1):
        kpara.append(2.*np.pi*jj/(rprime*(2.*N_E-2.)*deltanu))
        power.append([])

    A=np.zeros([N_c,N_E])

    A[0:N_c,0]=1.0/(2.*N_E-2.)

    for ii in range (0,N_c,1):
        A[ii,N_E-1]=np.cos(np.pi*cl[ii,0]/cl[1,0])/(2.*N_E-2.)

    for ii in range (0,N_c,1):
        for jj in range (1,N_E-1,1):
            A[ii,jj]=2.*np.cos(cl[ii,0]*jj*np.pi/(cl[1,0]*(N_E-1.)))/(2.*N_E-2.)

    A=np.matrix(A)
    ADag=A.transpose()

    Noise=np.identity(N_c)
    N=np.matrix(Noise)

    X=np.zeros([N_c,1])
    pk=np.zeros([N_E,1])

    for ii in range (nbin_kper):
        v=0.
        for jj in range (N_c):
            for kk in range (N_c):
                N[jj,kk]=var[jj,kk,ii]
                v+=var[jj,kk,ii]
        N=N/v

        Ninv=np.linalg.inv(N)
        AdagNinv=ADag*Ninv
        AdagNinvAwinv=np.linalg.inv(ADag*Ninv*A)
        fac=AdagNinvAwinv*AdagNinv

        X[:,0]=win*cl[:,ii+1]
        pk=fac*X
        for jj in range (0,N_E,1):
#			print('%e\t%e\t%e'%(2.*np.pi*jj/(rprime*(2.*N_c-2.)*deltanu),ell[ii]/r,r**2.*rprime*deltanu*pk[jj]))
            power[jj].append(r**2.*rprime*deltanu*pk[jj])

    power=np.array(power)

    np.save(OUTNAME+'kper',kper)
    np.save(OUTNAME+'kpara',kpara)
    np.save(OUTNAME+'pk',power)

    return kper, kpara, power


def binkpara(kper,kpara,power,nbin_kpara,OUTNAME):

    #### bin along k_parallel ####

    nbin_kper=len(kper)
    N_E=len(kpara)

    if(mode=='lin'):
        kv=(kmax-kmin)/(1.*nbin_kpara)
#         print(kv)
    if(mode=='log'):
        kv=math.log(kmax/kmin,10)/(1.*nbin_kpara)


    kpara_bin=np.zeros(nbin_kpara)
    wt=np.zeros(nbin_kpara)
    val=np.zeros([nbin_kpara,nbin_kper])

    for ii in range (0,N_E,1):
        wg=1.
        if(mode=='lin'):
            binno=int(math.floor((kpara[ii]-kmin)/kv))
        if(mode=='log'):
            if(kpara[ii]==0.):
                binno=0
            else:
                binno=int(math.floor(math.log(kpara[ii]/kmin,10)/kv))
#         if(binno>0):
#             binno=binno
#         else:
#             binno=0
#         print(binno)
        if(binno>=0 and binno<10):
#             print("inside",binno)
            kpara_bin[binno]+=(wg*kpara[ii])
            wt[binno]+=wg
            for jj in range (0,nbin_kper,1):
                val[binno,jj]+=(wg*power[ii,jj])

    for ii in range (0,nbin_kpara,1):
        if wt[ii]>0:
            val[ii,:]/=wt[ii]
            kpara_bin[ii]/=wt[ii]
#         print(kpara_bin)
    kparallel=[]
    value=[]

    for ii in range (0,nbin_kpara,1):
        if(kpara_bin[ii]!=0.):
            kparallel.append(kpara_bin[ii])
            value.append([])
            for jj in range (0,nbin_kper,1):
                value[len(value)-1].append(val[ii,jj])

    value=np.array(value)

#	print(kparallel)
#	print(value)

#	np.save(OUTNAME+'kpara_bin',kparallel)
#	np.save(OUTNAME+'pk_bin',value)

    return kparallel, value



def plot_heatmap(kper,kparallel,value,OUTNAME):

    nbin_kper=len(kper)
    nbin_kpara=len(kparallel)

    #plt.rc('text',usetex=True)

    X=np.zeros([nbin_kpara,nbin_kper])
    Y=np.zeros([nbin_kpara,nbin_kper])

    for ii in range (0,nbin_kpara,1):
        for jj in range (0,nbin_kper,1):
            X[ii,jj]=kper[jj]

    for ii in range (0,nbin_kpara,1):
        for jj in range (0,nbin_kper,1):
            Y[ii,jj]=kparallel[ii]

    fig=plt.figure()

    z_minimum=[]

    z_minimum = np.append(z_minimum,np.amin(np.abs(value)))    
    z_min=np.amin(abs(z_minimum))*10**6

    z_maximum=[]

    z_maximum = np.append(z_maximum,np.amax(np.abs(value)))   
    z_max = np.amax(abs(z_maximum))*10**6


    c = plt.pcolormesh(X, Y, np.abs(value)*10**6, cmap = 'RdYlGn_r',norm = LogNorm(vmin = z_min, vmax = z_max))

    plt.xscale('log')
    plt.yscale('log')

    plt.tick_params(labelcolor = 'k',size = 4,labelsize = 20,pad = 10,direction = 'in',which = 'both')

    plt.axis([X.min(), X.max(), Y.min(), Y.max()])

    # horizon
    Hori= r*kper/(rprime*nuc)
#	plt.plot(kper,Hori, marker='', ls='--', c='c')  

    # first null
    FN= (r/(rprime*nuc))*1.22*300.*kper/(nuc*Dia)
#	plt.plot(kper,FN, ls=':', c='k') 

    # colorbar axis... 

    cb = plt.colorbar(c)
#   cb.set_label(label = 'P ($k_{\perp}$,$k_{\parallel}$) mK$^2$ MPc$^{3}$',labelpad = 10,fontsize = 18)
    cb.ax.tick_params(labelcolor='k',size = 4, labelsize = 20, pad = 10, direction= 'in', which = 'both')

#     plt.xlabel(r'$k_{\perp}$ MPc$^{-1}$', fontsize = 18,labelpad=10)
#     plt.ylabel(r'$k_{\parallel}$ MPc$^{-1}$', fontsize = 18,labelpad=10)

    plt.show()
#	plt.savefig(OUTNAME)
