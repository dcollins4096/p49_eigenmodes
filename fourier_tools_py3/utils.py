import numpy as na
import numpy.fft as fpack
import pencil as pc
from .fourier_filter import FourierFilter
from .fft_shear import fft_shear, ifft_shear

def pc_en_spec(var,magnetic=False):
    """shell averaged energy power spectrum on pencil code data
    """
    if magnetic:
        data = (pc.curl(var.f[4:7],var.dx,var.dy,var.dz))[:,3:-3,3:-3,3:-3]
    else:
        data = var.f[0:3,3:-3,3:-3,3:-3]
    uuhat = []
    for i in range(3):
        if hasattr(var,'deltay'):
            uuhat.append(fft_shear(data[i],var.deltay,var.x[3:-3],var.dx,var.dy,var.dz)/data.size)
        else:
            uuhat.append(fpack.fftn(data[i])/data[i].size)
    power = 0.5*(na.array([na.abs(b)**2 for b in uuhat])).sum(axis=0)
    ff = FourierFilter(power)

    return na.array([power[ff.get_shell(bin)].sum() for bin in range(ff.nx)])

def pc_hel_spec(var,magnetic=False):
    """shell averaged helicity power spectrum on pencil code data
    """
    if magnetic:
        # bb is field; aa is vector potential
        bb = (pc.curl(var.f[4:7],var.dx,var.dy,var.dz))[:,3:-3,3:-3,3:-3]
        aa = var.f[4:7,3:-3,3:-3,3:-3]
    else:
        # bb is vorticity; aa is velocity
        aa = var.f[0:3,3:-3,3:-3,3:-3]
        bb = (pc.curl(var.f[0:3],var.dx,var.dy,var.dz))[:,3:-3,3:-3,3:-3]

    helhat = []
    enhat = []
    for i in range(3):
        if hasattr(var,'deltay'):
            aahat = fft_shear(aa[i],var.deltay,var.x[3:-3],var.dx,var.dy,var.dz)/aa.size
            bbhat = fft_shear(bb[i],var.deltay,var.x[3:-3],var.dx,var.dy,var.dz)/bb.size
        else:
            aahat = fpack.fftn(aa[i])/aa[i].size
            bbhat = fpack.fftn(bb[i])/bb[i].size
        if magnetic:
            enhat.append(bbhat)
        else:
            enhat.append(uuhat)
        helhat.append(aahat*bbhat.conj() + aahat.conj()*bbhat)

    hpower = 0.5*(na.array(helhat)).sum(axis=0)
    epower = 0.5*(na.array([na.abs(b)**2 for b in enhat])).sum(axis=0)
    ff = FourierFilter(hpower) #they both have the same shape, so this is fine
    k = ff.get_shell_k(dx=var.dx)    
    hspec = na.array([hpower[ff.get_shell(bin)].sum() for bin in range(ff.nx)])
    espec = na.array([epower[ff.get_shell(bin)].sum() for bin in range(ff.nx)])

    return k, hspec, espec

def pc_calc_spec(data):
    """compute a shell averaged power spectrum for an arbitrary field.

    """
    power = 0.5 * na.abs(fpack.fftn(data)/data.size)**2
    ff = FourierFilter(power)

    return na.array([power[ff.get_shell(bin)].sum() for bin in range(ff.nx)])


def pc_calc_Ak(var, k):
    """calculate real space filtered vector potentials containing only 

    k-dk/2 <= k < k+dk/2

    components in fourier space.

    inputs
    ------
     var -- pencil code var object
     k   -- integer wavenumber in units of k1 = 2*pi/L where L is the *smallest* 
            length of the box

    """
    aa = var.f[4:7]

    aa_k = na.zeros_like(aa)
    for i in range(3):
        if hasattr(var,'deltay'):
            ff = FourierFilter(aa[i],fft_func=fft_shear,ifft_func=ifft_shear, deltay=var.deltay, x=var.x,dx=var.dx,dy=var.dy,dz=var.dz)
        else:
            ff = FourierFilter(aa[i])
        aa_k[i] = ff.filter_field(k)
    
    return aa_k

def pc_calc_helk(var,k):
    """return the helicity density at k, computed in real space

    H_k = < A_k * B_k> V/dk

    and 

    M_k  = < B_k**2 > V/dk (magnetic energy spectrum)

    since the latter comes for free

    see Brandenburg & Subramanian (2005) Phys Rep pp30-31 for derivation
    """

    aa_k = pc_calc_Ak(var,k)
    bb_k = pc.curl(aa_k,var.dx,var.dy,var.dz)
    
    hel = (pc.dot(aa_k,bb_k))[3:-3,3:-3,3:-3]
    mag  = (pc.dot2(bb_k))[3:-3,3:-3,3:-3]
    vol = (var.x[-3]-var.x[3])*(var.y[-3]-var.y[3])*(var.z[-3]-var.z[3])

    # NB: this assumes cubic zones
    dk = FourierFilter(aa_k[0]).dk * 2*na.pi/var.dx
    if k == 0:
        k = 1

    # i think it should be k, and not dk; but i'm not sure yet
    return hel.mean()/k, mag.mean()/(2.*k)
    #return hel.mean()*vol/dk, mag.mean()*vol/dk

def pc_calc_helk_all(var):
    """return the helicity densities at all k

    """

    ff = FourierFilter(var.f[0,3:-3,3:-3,3:-3])
    hspec = []
    mspec = []

    for i in range(ff.nx):
        print("k = %i" % i)
        h,m = pc_calc_helk(var,i)
        hspec.append(h)
        mspec.append(m)

    hspec = na.array(hspec)
    mspec = na.array(mspec)
    k = ff.get_shell_k(dx=var.dx)
    data = na.array([k,hspec,mspec])

    return data
