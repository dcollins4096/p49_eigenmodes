import numpy as na
import numpy.fft as fpack

class FourierFilter(object):
    """FourierFilter takes a dataset and filters it in the Fourier domain. 

    Internally, the bins are calculated in wavenumber space in units
    of the nyquist wavenumber. However, when asking for a shell in
    fourier space, one identifies it by the *bin number*, which is an
    integer between zero and and the number of grid points in the
    shortest dimension of the box. Physically, of course, this means
    that long wavenumbers are missing, not short ones: the nyquist
    wave number in a box that is not uniform in length is the same in
    all directions. This will work fine AS LONG AS ALL DIMENSIONS HAVE
    EQUALLY SPACED GRIDS.

    """
    def __init__(self, data, fft_func=fpack.fftn,ifft_func=fpack.ifftn, **kwargs):
        self.data = data
        self.debug = False
        self._init_k()
        self._data_hat = []
        self._extra_fft_args = kwargs
        self._fftn = fft_func
        self._ifftn = ifft_func

    def _init_k(self):
        """this is done entirely in dimensionless wavenumber/frequency
        space, with nyquist frequency as 1.

        """
        self._kk = []
        l = len(self.data.shape)
        for i,dim in enumerate(self.data.shape):
            sl = i*(1,)+(dim,)+(l-i-1)*(1,)
            k = fpack.fftfreq(dim)
            k.resize(sl)
            self._kk.append(k)

        self.kradius = na.sqrt(na.sum((k**2 for k in self._kk)))

        shortind = (na.array([a.size for a in self._kk])).argmin()
        self.bins = self._kk[shortind].flatten()
        self.nx = self.bins.size/2
        self.dk = self.bins[1] - self.bins[0]

    def get_shell_k(self, dx=0):
        """returns all wavenumbers in the shell bins.
        
        if dx is set to the grid spacing, returned wavenumbers will be
        in code units

        """
        if dx != 0:
            k_scale = (2.*na.pi/dx)
        else:
            k_scale = 1.

        return self.bins[0:self.nx] * k_scale

    def get_shell(self,nk):
        """returns an boolean index array of dimensions equal to that
        of the data. when used to index fft'd data, will return only
        modes that fall in the shell.

        nk: the integer bin number, 0 <= nk <= Nbins, where Nbins is
        half the number of grid points in the SHORTEST direction of
        the box (Nbins corresponds to the nyquist wavenumber)
        """
        # this assumes equal width bins!

        bin_center = nk*self.dk
        bins = self.bins[0:self.nx]

        lower = bin_center - self.dk/2.
        upper = bin_center + self.dk/2.
        if self.debug:
            print lower, upper, bin_center
        return ((self.kradius >= lower) & (self.kradius <= upper))

    @property
    def fft(self):
        if len(self._data_hat) == 0:
            self._data_hat = self._fftn(self.data,**self._extra_fft_args)
        
        return self._data_hat
    
    def ifft(self, data_hat):
        return self._ifftn(data_hat,**self._extra_fft_args)

    def hi_pass(self, kthresh):
        pass
    
    def lo_pass(self, kthresh):
        pass

    def filter_field(self, k, dk=1.):
        shell = self.get_shell(k)
        hat = (self.fft).copy()
        hat[~shell] = 0
        return (self._ifftn(hat, **self._extra_fft_args)).real

if __name__ == "__main__":
    import pencil as pc
    import pylab as P
    import os
    import cPickle
    from utils import pc_en_spec
    datadir = os.path.expanduser('~/vc/pencil-code/samples/helical-MHDturb/data/')
    t,pcpow = pc.read_power('power_kin.dat',datadir=datadir)
    var = pc.read_var(datadir=datadir)
    power2 = 0.5*abs(fpack.fftn(pc.dot2(var.f[0:3,3:-3,3:-3,3:-3]))/var.f[0,3:-3,3:-3,3:-3].size)
    spec = pc_en_spec(var)
    q = FourierFilter(power2)
    spec2 = na.array([power2[q.get_shell(bin)].sum() for bin in range(q.nx)])
    P.loglog(pcpow[-1,:],'rx',markersize=10)
    P.loglog(spec)
    P.loglog(spec2)
    P.show()
