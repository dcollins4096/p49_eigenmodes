#!/usr/bin/env python
import numpy as N
# try:
#     import scipy.fftpack as fpack
# except ImportError:
#     import numpy.fft as fpack
import numpy.fft as fpack

__version__ = "$Id$"

def kvectors(uu,deltay=0.,Lx=1.,dx=0.,dy=0.,dz=0.):
    """construct wavevectors for a given shape input array.
    
    deltay -- shear distance from pencil code. if given, create
    sheared (Eularian) wave numbers. for Lagrangian k, just use deltay=0.
    
    if dx = dy = dz = 0., then construct dimensionless wavevectors (default).
    """
    sz = uu.shape
    nz = sz[0]
    ny = sz[1]
    nx = sz[2]
    
    if nx%2 == 0:
        #even
        xmax_index  = nx/2
        xmax_index1 = nx/2
    else:
        #odd
        xmax_index = (nx-1)/2
        xmax_index1 = xmax_index + 1 
    if ny%2 == 0:
        #even
        ymax_index  = ny/2
        ymax_index1 = ny/2
    else:
        ymax_index = (ny-1)/2
        ymax_index1 = ymax_index + 1
    if nz > 7:
        if nz%2 == 0:
            #even
            zmax_index  = nz/2
            zmax_index1 = nz/2
        else:
            zmax_index = (nz-1)/2
            zmax_index1 = zmax_index + 1

        kz_,ky_,kx_ = N.ogrid[-zmax_index:zmax_index1,
                              -ymax_index:ymax_index1,
                              -xmax_index:xmax_index1]
        kz = 2*N.pi*fpack.fftshift(kz_)
        ky = 2*N.pi*fpack.fftshift(ky_)
        kx = 2*N.pi*fpack.fftshift(kx_)
        if dx == dy == dz == 0.:
            pass
        else:
            kx /= dx * float(nx) 
            ky /= dy * float(ny) 
            kz /= dz * float(nz)
    else:
        # 2D case
        ky_,kx_ = N.ogrid[-ymax_index:ymax_index1,
                          -xmax_index:xmax_index1]
        ky = fpack.fftshift(ky_)
        kx = fpack.fftshift(kx_)
        kz = 0
        if dx == dy == dz == 0.:
            pass
        else:
            kx /= dx * float(ny) 
            ky /= dy * float(nx) 
            
    #add shear
    if ky.size > kx.size:
        diff = ky.size - kx.size 
        kx = kx + deltay/Lx*ky[:,diff/2:-diff/2,:].flatten()
    else:
        kx = kx + deltay/Lx*ky.flatten()

    return [kz,ky,kx]

def fft_shear(uu,deltay,x,dx,dy,dz,ireverse=False):
    """take a fourier transform in a shearing frame.

    uu is a 3D scalar field
    
    use fourier interpolation methods to apply the shear.
    """
    Lx = x[-1] - x[0]
    kz,ky,kx = kvectors(uu,deltay,Lx,dx,dy,dz)

    # construct x dependent "shear" wavenumber. must be a 3D array in
    # order to multiply it out with ky. use ogrid to create dummy y
    # and z directions.
    deltay_x = -deltay/Lx * (x-(x[0]+Lx/2))

    z_,y_,x_ = N.ogrid[0:uu.shape[0],
                       0:uu.shape[1],
                       0:uu.shape[2]]

    if (not ireverse):
        # do y direction first so shear can be applied via fourier
        # interpolation
        uu_hat = fpack.fft(uu,axis=1)
        uu_hat *= N.exp(1j*ky*(0.*z_+0.*y_+deltay_x))
        
        # x transform, now that shear is applied and domain is periodic
        uu_hat = fpack.fft(uu_hat,axis=2)
        
        # z transform
        uu_hat = fpack.fft(uu_hat,axis=0)

    else:
        #shift back
        #uu *= N.exp(-1j*ky*(0.*z_+0.*y_+deltay_x))
        #uu_hat = fpack.fftn(uu)
        uu_hat = fpack.ifft(uu,axis=2)
        uu_hat = fpack.ifft(uu_hat,axis=0)
        uu_hat *= N.exp(-1j*ky*(0.*z_+0.*y_+deltay_x))
        uu_hat = fpack.ifft(uu_hat,axis=1)
        
    return uu_hat

def ifft_shear(uu,deltay,x,dx,dy,dz):
    return fft_shear(uu,deltay,x,dx,dy,dz,ireverse=True)

def power_spec(data,dx=0.,dy=0.,dz=0.,deltay=0.,x=0,
               ampl_spec=False, check_parseval=False):
    """compute a power spectrum. if ampl_spec is true, then return the amplitude spectrum, 

    """
    if deltay != 0:
        fft = fft_shear(data,deltay,x,dx,dy,dz)
    else:
        fft = fpack.ifftn(data)
    if ampl_spec:
        pwr = N.abs(fft)
    else:
        pwr = (fft*fft.conj())
    if check_parseval:
        print("real sum: %8.5e" % (data**2).sum())
        print("k space sum: %8.5e" %(pwr.sum()*data.size))
    return pwr.real
    

def power_spec_1D(data,dx=0.,dy=0.,dz=0.,deltay=0.,x=0, ampl_spec=False):
    """output 1D spectra kx,ky,kz.
    
    """
    power = power_spec(data,dx,dy,dz,deltay,x,ampl_spec=ampl_spec)
    nz,ny,nx = power.shape
    # calculate a "two-sided" PSD
    power_x  = 2.*power.mean(axis=0).mean(axis=0)[0:nx/2]
    power_y  = 2.*power.mean(axis=2).mean(axis=0)[0:ny/2]
    power_z  = 2.*power.mean(axis=2).mean(axis=1)[0:nz/2]

    return [power_z,power_y,power_x]
