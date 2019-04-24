
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
def SpiralQU(Nzones=512):
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
    N = np.array([Nzones,Nzones],dtype = int32)
    size2d = np.array([xmax,ymax])
    Delta = size2d/N
    Deltal = cmbtools.Delta2l(Delta,N)
    lmax = 10000
    lbins = np.linspace(0,lmax,50)
    lcent = lbins[:-1] + np.diff(lbins)/2.
    Q =np.zeros(N)
    U =np.zeros(N)
    Hx =np.zeros(N)
    Hy =np.zeros(N)
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
            Hx[p] = np.cos(psi[i])
            Hy[p] = np.sin(psi[i])

    from scipy.ndimage.filters import gaussian_filter
    Q = gaussian_filter(Q,1.0/60.*pi/180/Delta[0],mode='wrap')
    U = gaussian_filter(U,1.0/60.*pi/180/Delta[0],mode='wrap')
    return {'Q':Q,'U':U, 'Hx':Hx,'Hy':Hy}

def QU2EB_take2(self,Q,U,EB_KLUDGE=1):
    if not Q.flags['C_CONTIGUOUS']:
        Q = np.ascontiguousarray(Q)
    if not U.flags['C_CONTIGUOUS']:
        U = np.ascontiguousarray(U)
    N = array(shape(Q),dtype = int32)
    xsize = 5 * pi / 180
    size2d = array([xsize,xsize])
    Delta = size2d/N

    Deltal = cmbtools.Delta2l(Delta,N)

    Qharm = cmbtools.map2harm(Q,Delta)
    Uharm = cmbtools.map2harm(U,Delta)

    Eharm, Bharm = cmbtools.QU2EB(Qharm,EB_KLUDGE*Uharm,Deltal)

    E = cmbtools.harm2map(Eharm,Delta)
    B = cmbtools.harm2map(Bharm,Delta)

    lmax = Deltal[0]*N[0]
    lbins = linspace(0,lmax,100)
    lcent = lbins[:-1] + diff(lbins)/2.
    return {'Q':Q,'U':U,'E':E,'B':B,'Eh':Eh,'Bh':Bh}

    #ClEE = cmbtools.harm2cl(self.Eharm,Deltal,lbins)
    #ClBB = cmbtools.harm2cl(self.Bharm,Deltal,lbins)

class QU_TO_EB():
    def __init__(self,Q,U,EB_KLUDGE=1):
        if not Q.flags['C_CONTIGUOUS']:
            Q = np.ascontiguousarray(Q)
        if not U.flags['C_CONTIGUOUS']:
            U = np.ascontiguousarray(U)
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
Spiral=SpiralQU(Nzones=64)
QUEB=QU_TO_EB(Spiral['Q'],Spiral['U'])
PlotQUEB(Spiral['Q'],Spiral['U'],QUEB.E,QUEB.B,outname='Spiral_64.png')
print("butts")
def plot_HxHy(Spiral):
    plt.close('all')
    args={'interpolation':'nearest','origin':'lower'}
    fig, axes = plt.subplots(2,1,sharex=True,sharey=True)
    Y,X = np.mgrid[0:1:64j,0:1:64j]
    args={'interpolation':'nearest','origin':'lower','extent':[0,1,0,1]}
    axes[0].imshow(Spiral['Hx'],**args)
    axes[0].set_title('Hx')
    axes[1].imshow(Spiral['Hy'],**args)
    axes[1].set_title('Hy')
    axes[1].streamplot(X,Y,Spiral['Hx'],Spiral['Hy'],color='k')
    fig.savefig('HxHy.png')
def ToEnzo(Spiral,prefix="../spiral_test/test_", axis=[2]):
    import enzo_write
    reload(enzo_write)
    enzo_field_list = ['d','vx','vy','vz','hx','hy','hz','p']
    map_to_label ={'d':'density','vx':'x-velocity','vy':'y-velocity','vz':'z-velocity',
                   'hx':'Bx','hy':'By','hz':'Bz','e_specific':'TotalEnergy','p':'GasPressure'}
    Hx = Spiral['Hx']
    Hy = Spiral['Hy']
    size = Hx.shape[0]
    data={}
    u0 = {'d':1,'vx':0,'vy':0,'vz':0,'hx':0,'hy':0,'hz':0,'e_specific':1, 'p':1.}
    for field in enzo_field_list:
        data[field]=np.zeros([size]*3) + u0[field]
    #data['hx'][:,:,size//2]=Hx
    #data['hy'][:,:,size//2]=Hy

    if 2 in axis:
        data['hx'][:,:,size//2]=Hx.swapaxes(0,1)
        data['hy'][:,:,size//2]=Hy.swapaxes(0,1)

    if 0 in axis:
        data['hy'][0,:,:]=Hx.swapaxes(0,1)
        data['hz'][0,:,:]=Hy.swapaxes(0,1)

    if 1 in axis:
        data['hx'][:,0,:]=Hy#.swapaxes(0,1)
        data['hz'][:,0,:]=Hx#.swapaxes(0,1)
    for field in enzo_field_list:
        this_filename = "%s%s_%d.h5"%(prefix,map_to_label[field], size)
        enzo_write.dump_h5(data[field],this_filename)
    return data
#data=ToEnzo(Spiral, prefix="../spiral_test_x1/test_", axis=[0])
#data=ToEnzo(Spiral, prefix="../spiral_test_y1/test_", axis=[1])
#data=ToEnzo(Spiral, prefix="../spiral_test_z1/test_", axis=[2])
def plot_HxHyHz_proj(data):
    plt.close('all')
    args={'interpolation':'nearest','origin':'lower'}
    fig, axes = plt.subplots(3,3,sharex=True,sharey=True)

    for ax in [0,1,2]:# [0,1,2]:
        the_x = ['hy','hz','hx'][ax]
        the_y = ['hz','hx','hy'][ax]
        the_z = ['hx','hy','hz'][ax]
        this_hx = np.sum(data[the_x], axis=ax)
        this_hy = np.sum(data[the_y],axis=ax)
        this_hz = np.sum(data[the_z],axis=ax)
        Nx = this_hx.shape[0]
        this_hx.shape = (Nx,Nx)
        this_hy.shape = (Nx,Nx)
        this_hz.shape = (Nx,Nx)
        Y,X = np.mgrid[0:1:Nx*1j,0:1:Nx*1j]
        args={'interpolation':'nearest','origin':'lower','extent':[0,1,0,1]}
        axes[ax][0].imshow(this_hx,**args)
        axes[ax][0].set_title(the_x + "proj %d"%ax)

        axes[ax][1].imshow(this_hy,**args)
        axes[ax][1].set_title(the_y)
        axes[ax][1].streamplot(X,Y,this_hx,this_hy,color='k')

        axes[ax][2].imshow(this_hz,**args)
        axes[ax][2].set_title(the_z)
        fig.savefig('HxHyHz_proj.png')
def plot_HxHy_proj(data):
    plt.close('all')
    args={'interpolation':'nearest','origin':'lower'}
    fig, axes = plt.subplots(2,1,sharex=True,sharey=True)
    this_hx =np.sum(data['hx'],axis=2)
    this_hy = np.sum(data['hy'],axis=2)
    Nx = this_hx.shape[0]*1j
    Y,X = np.mgrid[0:1:Nx,0:1:Nx]
    args={'interpolation':'nearest','origin':'lower','extent':[0,1,0,1]}
    axes[0].imshow(this_hx,**args)
    axes[0].set_title('Hx')
    axes[1].imshow(this_hy,**args)
    axes[1].set_title('Hy')
    axes[1].streamplot(X,Y,this_hx,this_hy,color='k')
    fig.savefig('HxHy_proj.png')
#plot_HxHyHz_proj(data)
#import h5py
if 1:
    bxname = 'test_Bx_64.h5'
    bx = h5py.File('../spiral_test/%s'%bxname)[bxname].value.swapaxes(0,2)
    byname = 'test_By_64.h5'
    by = h5py.File('../spiral_test/%s'%byname)[byname].value.swapaxes(0,2)
    
if 0:
    fptr = h5py.File('../OT/DD0000/data0000.cpu0000')['Grid00000001']
    bx = fptr['Bx'].value.swapaxes(0,2)
    by = fptr['By'].value.swapaxes(0,2)
    #fptr.close()
#plot_HxHy_proj({'hx':bx,'hy':by})








