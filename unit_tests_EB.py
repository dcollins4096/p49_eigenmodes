
import os
import sys
import pdb
#sys.path.append('cmbtools')
#sys.path.append('../cmbtools_nofits')
sys.path.append('/Users/dcollins/local-other-2018-01-05/cmbtools_nofits')
if os.path.isdir('code') :
    sys.path.append('code')
import cmbtools
#import pyfits
import astropy.io.fits as pyfits
from pylab import *
import re
import glob
if 'frame' not in dir():
    frame = 0
def SpiralQU():
    simname = 'spiral'
    xmax = 5 *np.pi/180
    ymax = 5 *np.pi/180
    a = 1.0 * np.pi/180
    b = 0.15
    thetamin = -np.pi/2
    thetamax = 3*np.pi/2
    Nsteps = 300
    r = np.linspace(a*np.exp(b*thetamin),a*np.exp(b*thetamax),Nsteps)
    theta = np.log(r/a)/b if (b!=0.0) else np.linspace(thetamin,thetamax,Nsteps)
    angtan = np.pi/2 if (b==0.0) else np.arctan(1/b) # angle between tangent line and radial direction
    angpol = theta + angtan - np.pi/2
    xc = 0.5*xmax
    yc = 0.5*ymax
    x = r*np.cos(theta) + xc
    y = r*np.sin(theta) + yc
    r0 = 0.1 * pi/180
    u = r0*np.cos(angpol)
    v = r0*np.sin(angpol)
    psi = angpol - np.pi/2
    N = np.array([512,512],dtype = int32)
    size2d = np.array([xmax,ymax])
    Delta = size2d/N
    Deltal = cmbtools.Delta2l(Delta,N)
    lmax = 10000
    lbins = np.linspace(0,lmax,50)
    lcent = lbins[:-1] + np.diff(lbins)/2.
    Q =np.zeros(N)
    U =np.zeros(N)
    # identify coordinates of places to fill in
    s = 1.0
    sres = 20
    for i in range(Nsteps) :
        crossx = u[i]*np.linspace(-s,s,sres)
        crossy = v[i]*np.linspace(-s,s,sres)
        cox = crossy/5
        coy = -crossx/5

        regx = (x[i] + [[ x0+x1 for x0 in crossx] for x1 in cox]).flatten()
        regy = (y[i] + [[ y0+y1 for y0 in crossy] for y1 in coy]).flatten()

        pixx = np.array(regx/Delta[1],dtype=int)
        pixy = np.array(regy/Delta[0],dtype=int)

        for p in zip(pixy,pixx):
            Q[p] =np.cos(2*psi[i])
            U[p] =np.sin(2*psi[i])

    from scipy.ndimage.filters import gaussian_filter
    Q = gaussian_filter(Q,1.0/60.*pi/180/Delta[0],mode='wrap')
    U = gaussian_filter(U,1.0/60.*pi/180/Delta[0],mode='wrap')
    return {'Q':Q,'U':U}

class QU_TO_EB():
    def __init__(self,Q,U,EB_KLUDGE=1):
        self.Q=Q
        self.U=U
        N = array(shape(Q),dtype = int32)
        xsize = 5 * pi / 180
        size2d = array([xsize,xsize])
        Delta = size2d/N

        Deltal = cmbtools.Delta2l(Delta,N)

        self.Qharm = cmbtools.map2harm(Q,Delta)
        self.Uharm = cmbtools.map2harm(U,Delta)

        self.Eharm, self.Bharm = cmbtools.QU2EB(self.Qharm,EB_KLUDGE*self.Uharm,Deltal)

        self.E = cmbtools.harm2map(self.Eharm,Delta)
        self.B = cmbtools.harm2map(self.Bharm,Delta)

        lmax = Deltal[0]*N[0]
        lbins = linspace(0,lmax,100)
        lcent = lbins[:-1] + diff(lbins)/2.

        self.ClEE = cmbtools.harm2cl(self.Eharm,Deltal,lbins)
        self.ClBB = cmbtools.harm2cl(self.Bharm,Deltal,lbins)

def PlotQUEB(Q,U,E,B, outname = 'QUEB_plot.png'):
    plt.clf()
    fig, axes = plt.subplots(2,2,sharex=True,sharey=True)
    args={'interpolation':'nearest','origin':'lower'}
    axes[0][0].imshow(Q,**args)
    axes[0][0].set_title('Q')
    axes[1][0].imshow(U,**args)
    axes[1][0].set_title("U")
    oot=axes[0][1].imshow(E,**args)
    fig.colorbar(oot)
    axes[0][1].set_title("E")
    axes[1][1].imshow(QUEB.B,**args)
    axes[1][1].set_title("B")
    fig.savefig(outname)
    print(outname )
#Spiral=SpiralQU()
#QUEB=QU_TO_EB(Spiral['Q'],Spiral['U'])
#PlotQUEB(Spiral['Q'],Spiral['U'],QUEB.E,QUEB.B)
