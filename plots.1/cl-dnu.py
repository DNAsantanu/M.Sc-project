#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


# In[2]:


def collapse_channels(gv, nc):
    nb = np.intc(Nchan/nc)
    Nu=np.shape(gv)[0]
    cgv=np.zeros((Nu,Nu,nb), dtype=np.complex128)
    for ux in range (Nu):
        for vy in range (Nu):
            bb=0
            for ii in range(nb):
                aa = bb
                bb = aa+nc
                if np.count_nonzero(gv[ux,vy,aa:bb])>0:
                    cgv[ux,vy,ii] = np.sum(gv[ux,vy,aa:bb])/np.count_nonzero(gv[ux,vy,aa:bb])
                else:
                    cgv[ux,vy,ii]=0.
    return cgv


# In[3]:


def cross_corr_dnu(fvcg1,fvcg2,NC):
    # read and collapse channels
    vcg_XX = collapse_channels(fvcg1, NC)
    vcg_YY = collapse_channels(fvcg2, NC)
    
    # dimensions of data
    Nu=np.shape(vcg_XX)[0]
    n_ave=np.shape(vcg_XX)[2]
    print("Nu=%d, Nc=%d"%(Nu,n_ave))
    
    eg=np.zeros(vcg_XX.shape)
    
    for ii in range(Nu):
        for jj in range(Nu):
            for ch1 in range (n_ave):
                for ch2  in range (ch1,n_ave):
                    eg[ii,jj,ch2-ch1]+=(np.real((vcg_XX[ii,jj,ch1]*vcg_YY[ii,jj,ch2].conjugate())+(vcg_YY[ii,jj,ch1]*vcg_XX[ii,jj,ch2].conjugate()))/2.)
            
    return(eg)


# In[4]:


def corr_dnu(vcg,ag):
    # dimensions of data
    Nu=np.shape(vcg)[0]
    n_ave=np.shape(vcg)[2]
    
    eg=np.zeros(vcg.shape)
    
    for ii in range(Nu):
        for jj in range(Nu):
            for ch1 in range (n_ave):
                for ch2  in range (ch1,n_ave):
                    eg[ii,jj,ch2-ch1]+=np.real(vcg[ii,jj,ch1]*vcg[ii,jj,ch2].conjugate())

            eg[ii,jj,:]-=np.real(ag[ii,jj,:])    #  Subtract the self-terms for C_{\ell}
            
    return(eg)


# In[5]:


def bin_cl_dnu(GV,GVconst,OUTNAME):

	##################################

	if(mode=='lin'):
		kv=(1.*Nbin)/(Umax-Umin)
	if(mode=='log'):
		kv=(1.*Nbin)/math.log((Umax/Umin),10)

	#### read input header #####

	Nu=GV.shape[0]
	Nv=GV.shape[1]
	n_ave=GV.shape[2]
	Ng=int((Nu-1)/2)
	Nuv= Nu*Nv

	print("%d %d %d %d"%(Nu,Nv,n_ave,Ng))
	print(GV.shape)

	#### read Mg header #####

	Nu=GVconst.shape[0]
	Nv=GVconst.shape[1]

	print("%d %d"%(Nu,Nv))
	print(GVconst.shape)


	binval=np.zeros((n_ave+1, Nbin+1),dtype=np.float64)
	weight=np.zeros((n_ave, Nbin), dtype=np.float64)
	uval=np.zeros(Nbin, dtype=np.float64)
	num=np.zeros(Nbin, dtype=np.float64)

	for ii in range (0,Nu,1):
		for jj in range (Ng,Nv,1):
			ii1=ii-Ng
			jj1=jj-Ng
			UvalGrid=dU*(ii1*ii1+jj1*jj1)**0.5
			if(UvalGrid>= Umin and UvalGrid< Umax):
				if(mode=='lin'):
					NUGrid=int(math.floor(kv*(UvalGrid-Umin)))
				if(mode=='log'):
					NUGrid=int(math.floor(kv*math.log((UvalGrid/Umin),10)))
				if(NUGrid>0):
					NUGrid=NUGrid
				else:
					NUGrid=0
				for chan in range (0,n_ave,1):
					Ag2=GV[ii,jj,chan]
					Cg2=GVconst[ii,jj,chan]
					wg=GVconst[ii,jj,chan]
					#wg=1.
					if(Cg2>cutoff):
						binval[chan+1,NUGrid+1]+=(wg*Ag2/Cg2)
						weight[chan,NUGrid]+=wg
						if(chan==0):
							uval[NUGrid]+=(UvalGrid*wg)
							num[NUGrid]+=1.

	for ii in range (0,Nbin,1):
		for chan in range (0,n_ave,1):
			if(weight[chan,ii]>0.):
				binval[chan+1,ii+1]/=weight[chan,ii]
			else:
				binval[chan+1,ii+1]=1.e-99

	#### Write the output #####

	fpout=open(OUTNAME+'Power_2D.dat','w')

	for ii in range (0,Nbin,1):
		if(weight[0,ii]>0.):
			fpout.write('%d\t%e\t%e\n'%(int((2.*np.pi*uval[ii])/weight[0,ii]),binval[0+1,ii+1],num[ii]))
			binval[0,ii+1]=int((2.*np.pi*uval[ii])/weight[0,ii])
		else:
			fpout.write('%e\t%e\t%e\n'%(1.e-99,binval[0+1,ii+1],num[ii]))
			binval[0,ii+1]=1.e-99

	fpout.close()

	for ii in range (0,n_ave,1):
#		binval[ii+1,0]=ii*CDELT[2]
		binval[ii+1,0]=ii
		
	binval=np.array(binval)

	np.save(OUTNAME+'cl_dnu_py',binval)

	return 0


# In[6]:


def bin_nua_nub_to_dnu(GV,ell,OUTNAME):
    
    Nbin=GV.shape[2]
    n_ave=GV.shape[0]
    
    binval=np.zeros((n_ave+1, Nbin+1),dtype=np.float64)
    
    binval[0,1:]=ell
    binval[1:,0]=np.arange(0,n_ave,1)
    
    for ii in range (Nbin):
        for chan in range (n_ave):
#            if np.count_nonzero(np.diag(GV[:,:,ii],k=chan))>0:
#                binval[chan+1,ii+1]=np.sum(np.diag(GV[:,:,ii],k=chan))/np.count_nonzero(np.diag(GV[:,:,ii],k=chan))
#            else:
#                binval[chan+1,ii+1]=1.e-99
            binval[chan+1,ii+1]=np.sum(np.diag(GV[:,:,ii],k=chan))/(n_ave-chan)
    
    binval=np.array(binval)
    np.save(OUTNAME+'cl_nua_nub_dnu_py',binval)
    
    return 0


# # INPUTS

# In[7]:


#### INPUTS ################

Umax=477.0 # in lambda
Umin=80. # in lambda
Nbin=6
mode='lin'  #  'lin' or 'log'
cutoff=0.0
dU=4.57 # in lambda
Nchan=151

NCol= 1 # Channels to collapse
Nim = 1 # No. of realizations of Images
Nmg = 10 # No. of realizations of Mg

#############################


# # Data : $C_{\ell}(\Delta\nu)$

# In[8]:


# Correlate data

XX=np.load('Output/vcg_XX.npy')
XXsq=np.load('Output/vcg_sq_XX.npy')

Eg_XX=corr_dnu(XX,XXsq)


# Correlate Mg

Eg_mg_XX=np.zeros(Eg_XX.shape)

for ii in range (Nmg):
    XX=np.load('Output/Mg/vcg_mg_XX_'+str(ii+1)+'.npy')
    XXsq=np.load('Output/Mg/vcg_sq_mg_XX_'+str(ii+1)+'.npy')

    Eg_mg_XX+=corr_dnu(XX,XXsq)

    
# Avegrage Mg

Eg_mg_XX /= Nmg


# Bin data

OUTNAME='Output/MAPS/XX_'
bin_cl_dnu(Eg_XX,Eg_mg_XX,OUTNAME)


# In[ ]:




