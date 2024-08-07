"""import numpy as np
import multiprocessing as mp # multiprocessing
import matplotlib.pyplot as plt
import sys

def read_vcg(filename, nc, NC, xn):
    headersize = (nc+3)*np.float64().itemsize
#     aa = np.fromfile(filename, dtype='float64', count=headersize, offset=0)
    aa = np.fromfile(filename, dtype='float64', offset=0)
#     FWHM, dUg, f, nua = aa[0], aa[1], aa[2], aa[3:nc+3] # header
    FWHM, dUg, f, nua, bb = aa[0], aa[1], aa[2], aa[3:nc+3], aa[nc+3:]
    hinfo = (FWHM, dUg, f, nua, bb)
    aa = np.fromfile(filename, dtype='complex128', count=-1, offset=headersize)
#     Nu = int(np.sqrt(len(aa)/(3*nc)))
    Nu = int(np.sqrt(bb.shape[0]/(6*nc)))
    ee = 4*Nu*Nu*nc
    GV = np.reshape(bb[:ee],(2,Nu,Nu,nc,2))
    AG = np.reshape(bb[ee:],(Nu,Nu,nc,2))

    GVRR = GV[0,:,:,:,0] + 1j*GV[0,:,:,:,1]
#     GVLL = GV[1,:,:,:,0] + 1j*GV[1,:,:,:,1]

    AGRR = AG[:,:,:,0]
#     AGLL = AG[:,:,:,1]
#     NN = Nu*Nu*nc
#     GR = collapse(np.reshape(aa[:NN], (Nu, Nu, nc))[xn:Nu-xn, xn:Nu-xn, :], NC)
#     GL = collapse(np.reshape(aa[NN:2*NN], (Nu, Nu, nc))[xn:Nu-xn, xn:Nu-xn, :], NC)
#     AA = collapse(np.reshape(aa[2*NN:3*NN], (Nu, Nu, nc))[xn:Nu-xn, xn:Nu-xn, :], NC)
#     return hinfo, GR, GL, np.real(AA), np.imag(AA)
    return hinfo, GVRR, AGRR, GV, AG, bb, ee, Nu 


# hinfo, GVRR, AGRR, GV, AG, bb, ee, Nu  = read_vcg('/home/santanu-das/ska_data/combined_0_150', 151, 1, 0)
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


# Create a new figure
fig = plt.figure(figsize=(7,5))

for i in range(150):
    plt.clf()  # Clear previous frame
    plt.imshow(np.angle(GVRR[:,:,i]), cmap='viridis')  # Plot absolute values with viridis colormap
    plt.pause(0.01)
plt.show()