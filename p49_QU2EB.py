import os
import sys
import pdb
#sys.path.append('cmbtools')
#sys.path = ['../cmbtools_nofits']+sys.path
frbname = 'frbs'
#if os.path.isdir('code') :
#    sys.path.append('code')

import cmbtools
#import pyfits
import astropy.io.fits as pyfits
from scipy.optimize import leastsq
from pylab import *
import pylab
import re
import glob

def linear_chi2(par,x,y) :
    m = par[0]
    y0 = par[1]

    xcent = (x[0]+x[-1])/2

    ymodel = m*(x-xcent)+y0
#    ymodel = m*(x)+y0
    
    return(y - ymodel)



def slopes_powers(car,frame, n0=1,p=1, fitlmin=1e3, fitlmax=8e3, plot_format="pdf"):
    #base = '%s/FRBs/DD%04d_Cl%s_n0-%04d_p-%d.dat'%(car.directory,frame, "%s", n0, p)
    #output_base = '%s_DD%04d_n0-%04d_p-%d_%s.%s'%(car.outname,frame, n0, p,"%s","%s")
    base = '%s/%s/DD%04d_Cl%s.dat'%(car.directory,frbname,frame, "%s")
    output_base = '%s_DD%04d_%s.%s'%(car.outname,frame,"%s","%s")
    shortname = car.name
    shortname_n0p = "%s_n0-%04d_p-%d"%(shortname, n0, p)

#fitlmin = 5e3
#fitlmax = 3e4


    ClE = {}
    ClB = {}

    for ax in ['x','y','z'] :
        ell, ClE[ax], ClB[ax] = pylab.loadtxt(base % ax,usecols=[0,1,2],unpack=True)
        
    ellmask = (fitlmin < ell)*(ell < fitlmax)

    Eslope = {}
    Bslope = {}

    Eamp = {}
    Bamp = {}

    import pdb
    #pdb.set_trace()
    for ax in ['x','y','z'] :
        mb = np.array([0,0])
        ellmask_finite = np.logical_and(ellmask, ClE[ax] > 0)
        ellmask_finite = np.logical_and(ellmask_finite, ClB[ax] > 0)


        if ellmask_finite.sum() == 0:
            make_plots=False
            print( "Warning: No positive values of harmonics in target ell range")
            Bslope[ax]=np.nan
            Bamp[ax]=np.nan
            Eslope[ax]=np.nan
            Eamp[ax]=np.nan
            continue
        make_plots=True
        x =np.log10(ell[ellmask_finite])

        y = np.log10(ClE[ax][ellmask_finite])
        res = leastsq(linear_chi2, mb, args=(x,y) )
        Eslope[ax] = res[0][0]
        Eamp[ax] = pow(10,res[0][1])
        
        y = np.log10(ClB[ax][ellmask_finite])
        res = leastsq(linear_chi2, mb, args=(x,y) )
        Bslope[ax] = res[0][0]
        Bamp[ax] = pow(10,res[0][1])

        
    if make_plots:
        coldict = {'x':'r', 'y':'g', 'z':'b'}
        #values_dict = {}
        ##pickle_fname = 'p49_all_values.pickle'
        ##if glob.glob(pickle_fname) != []:
        ##    values_dict = fPickle.load(pickle_fname)
        #    #[shortname][frame][fit range, Eslope(3) Bslope(3) Eamp(3) Bamp(3)]
        #if not values_dict.has_key(shortname_n0p):
        #    values_dict[shortname_n0p]={}
        #if not values_dict[shortname_n0p].has_key(frame):
        #    values_dict[shortname_n0p][frame]={'fit_range':[fitlmin, fitlmax]}
        #    values_dict[shortname_n0p][frame]['Eslope'] = Eslope
        #    values_dict[shortname_n0p][frame]['Bslope'] = Bslope
        #    values_dict[shortname_n0p][frame]['Eamp'] = Eamp
        #    values_dict[shortname_n0p][frame]['Bamp'] = Bamp
        ##fPickle.dump(values_dict,pickle_fname)
            

# plots
        lfit = pow(10,(x[0]+x[-1])/2)

        fig = plt.figure()
        axis = fig.add_subplot(111)
        axis.plot(ell,Eamp['x']*(ell/lfit)**(-2.43), c=[0.5,0.5,0.5], label='%0.2f'%(-2.43))
        for ax in ['x','y','z'] :
            axis.plot(ell,ClE[ax],'-',label='E%s %0.2f'%(ax, Eslope[ax]), c=coldict[ax])
            axis.plot(ell,ClB[ax],'--',label='B%s %0.2f'%(ax, Bslope[ax]), c=coldict[ax])
            axis.plot(ell,Eamp[ax]*(ell/lfit)**Eslope[ax],'k:')
            axis.plot(ell,Bamp[ax]*(ell/lfit)**Bslope[ax],'k:')

            print('%s Eslope %f Bslope %f Eamp/Bamp %f' % (ax,Eslope[ax],Bslope[ax],Eamp[ax]/Bamp[ax]))
            
        axis.legend()
        axis.loglog()
        axis.set_xlabel(r'$\ell$')
        axis.axvspan(fitlmax,pylab.gca().get_xlim()[1],facecolor='k',alpha=0.1,ls=None)
        axis.axvspan(pylab.gca().get_xlim()[0],fitlmin,facecolor='k',alpha=0.1,ls=None)
        fig.savefig(output_base%("EE_BB",plot_format))

        fig = plt.figure()
        axis=fig.add_subplot(111)
        axis.plot(ell, np.zeros_like(ell)+2.0, c=[0.5]*3)
        for ax in ['x','y','z'] :
            axis.plot(ell,ClE[ax]/ClB[ax],'-',label=ax, c=coldict[ax])

        axis.set_title('ClE/ClB')
        axis.legend()
        axis.loglog()
        axis.set_xlabel(r'$\ell$')
        fig.savefig(output_base%("Ratio",plot_format))



        fig = plt.figure()
        axis=fig.add_subplot(111)
        ax = 'y'
        axis.plot(ell,ClE[ax],'-',label='E%s'%ax, c=coldict[ax])
        axis.plot(ell,ClB[ax],'--',label='B%s'%ax, c=coldict[ax])
        axis.plot(ell,Eamp[ax]*(ell/lfit)**Eslope[ax],'k:')
        axis.plot(ell,Bamp[ax]*(ell/lfit)**Bslope[ax],'k:')

        axis.set_xlim(500,1.3e4)

#legend()
        axis.loglog()

        axis.axvspan(pylab.gca().get_xlim()[0],fitlmin,facecolor='k',alpha=0.1,ls=None)
        axis.axvspan(fitlmax,pylab.gca().get_xlim()[1],facecolor='k',alpha=0.1,ls=None)

#fig.savefig('EBfit.pdf',bbox_inches='tight')
        axis.set_xlabel(r'$\ell$')
        fig.savefig(output_base%("EBfit",plot_format))

    return {'Eslope':Eslope,'Eamp':Eamp,'Bslope':Bslope,'Bamp':Bamp}

#if len(sys.argv) > 1:
#    rootdir = sys.argv[1]
#    frame = int(sys.argv[2])
#    Qlist = glob.glob(rootdir+'/DD%04d_Q[xyz]*.fits'%frame)
#    print Qlist
#else:
#    rootdir = './output'
#    Qlist = glob.glob(rootdir+'/*/DD*_Q[xyz]*.fits')
def EBfromQU(Q,U, return_quharm=False):

    if not Q.flags['C_CONTIGUOUS']:
        Q = np.ascontiguousarray(Q)
    if not U.flags['C_CONTIGUOUS']:
        U = np.ascontiguousarray(U)
    N = array(shape(Q),dtype = int32)
    xsize = 5 * pi / 180
    size2d = array([xsize,xsize])
    Delta = size2d/N

    print("N = ",N)
    print("Delta = ",Delta)

    Deltal = cmbtools.Delta2l(Delta,N)

    Qharm = cmbtools.map2harm(Q,Delta)
    Uharm = cmbtools.map2harm(U,Delta)
    Tharm = cmbtools.map2harm(Density,Delta)

    Eharm, Bharm = cmbtools.QU2EB(Qharm,Uharm,Deltal)

    E = cmbtools.harm2map(Eharm,Delta)
    B = cmbtools.harm2map(Bharm,Delta)
    if return_quharm:
        output = {'E':E,'B':B,'Eh':Eharm,'Bh':Bharm, 'Qh':Qharm, 'Uh':Uharm, 'Deltal':Deltal,'Delta':Delta,'Tharm':Tharm}
    else:
        output = E,B,Eharm,Bharm, Tharm
    return output
def QU2EB(rootdir,frame):
    Qlist = glob.glob(rootdir+'/DD%04d_Q[xyz]*.fits'%frame)
    Ulist = []
    for Qfile in Qlist:
        mo = re.match('(.*/DD[0-9]{4}_)Q([xyz].*)(.fits)',Qfile)
        Ufile = mo.group(1)+'U'+mo.group(2)+'.fits'
        Ulist.append(Ufile)

    QUlist = zip(Qlist,Ulist)

    for Qfile, Ufile in QUlist :

        Q = array(pyfits.open(Qfile)[0].data,dtype=double)
        U = array(pyfits.open(Ufile)[0].data,dtype=double)

        N = array(shape(Q),dtype = int32)
        xsize = 5 * pi / 180
        size2d = array([xsize,xsize])
        Delta = size2d/N
        Deltal = cmbtools.Delta2l(Delta,N)

        E,B, Eharm, Bharm, Tharm = EBfromQU(Q,U)

        lmax = Deltal[0]*N[0]
        lbins = linspace(0,lmax,100)
        lcent = lbins[:-1] + diff(lbins)/2.
        
        ClEE = cmbtools.harm2cl(Eharm,Deltal,lbins)
        ClBB = cmbtools.harm2cl(Bharm,Deltal,lbins)
        ClTE = cmbtools.harm2clcross_samegrid(Eharm, Tharm,Deltal,lbins)

        # write the output near the input
        mo = re.match('(.*/DD[0-9]{4}_)Q([xyz].*)(.fits)',Qfile)
        outroot = mo.group(1)
        outsuf = mo.group(2)

        Efile = outroot+'E'+outsuf+'.fits'
        Bfile = outroot+'B'+outsuf+'.fits'
        Clfile = outroot+'Cl'+outsuf+'.dat'

        print( Efile, Bfile, Clfile)

        hdu = pyfits.PrimaryHDU(E)
        hdulist = pyfits.HDUList([hdu])
        hdulist.writeto(Efile,overwrite=True)

        hdu = pyfits.PrimaryHDU(B)
        hdulist = pyfits.HDUList([hdu])
        hdulist.writeto(Bfile,overwrite=True)

        savetxt(Clfile, zip(lcent,ClEE,ClBB))

'''
# Plot Q/U
figure(figsize=[20,6])
subplot(1,4,1)
title('Q')
imshow(Q,interpolation='nearest')
colorbar()
subplot(1,4,2)
title('U')
imshow(U,interpolation='nearest')
colorbar()

# Plot E/B
subplot(1,4,3)
title('E')
imshow(E,interpolation='nearest')
colorbar()
subplot(1,4,4)
imshow(B,interpolation='nearest')
title('B')
colorbar()

#plot Cl
figure()
plot(lcent,ClEE,label='EE')
plot(lcent,ClBB,label='BB')
loglog()
legend()

show()
'''
'''
N = array([128,128],dtype = int32)
xsize = 1 * pi / 180
size2d = array([xsize,xsize])

Delta = size2d/N

print("N = ",N)
print("Delta = ",Delta)

Deltal = cmbtools.Delta2l(Delta,N)

print("Deltal = ",Deltal)

border = int(N[0]/10)
mask = zeros(N)
mask[border:-border,border:-border] = 1.0
#mask = ones(N)
fmask = sum(mask)/prod(N)

lmax = 10000
lbins = linspace(0,lmax,50)
lcent = lbins[:-1] + diff(lbins)/2.

# Make some E/B

E = gaussian_filter(normal(size=N),2,mode='wrap')
#B = zeros(N)
B = normal(size=N)
B = gaussian_filter(B,4,mode='wrap') - gaussian_filter(B,2,mode='wrap')


#E = zeros(N)#array([[cos(3*2*pi*(i+j)/N[1]) for i in range(N[1]) ] for j in range(N[0])])
#E = array([[cos(3*2*pi*(i+2*j)/N[1]) for i in range(N[1]) ] for j in range(N[0])])
# + sin(2*pi*i/N[1])

#B = array([[cos(3*2*pi*(i+j)/N[1]) for i in range(N[1]) ] for j in range(N[0])]);


# Now do the analysis

Eharm = cmbtools.map2harm(E,Delta)
Bharm = cmbtools.map2harm(B,Delta)

ClEE = cmbtools.harm2cl(Eharm,Deltal,lbins)
ClBB = cmbtools.harm2cl(Bharm,Deltal,lbins)

Qharm, Uharm = cmbtools.EB2QU(Eharm,Bharm,Deltal)

Q = cmbtools.harm2map(Qharm,Delta) * mask
U = cmbtools.harm2map(Uharm,Delta) * mask

Qharm2 = cmbtools.map2harm(Q,Delta)
Uharm2 = cmbtools.map2harm(U,Delta)

Eharm2, Bharm2 = cmbtools.QU2EB(Qharm2,Uharm2,Deltal)

#ClB = cmbtools.harm2cl(Bharm2)
#ClE = cmbtools.harm2cl(Eharm2)
ClEEout = cmbtools.harm2cl(Eharm2,Deltal,lbins)
ClBBout = cmbtools.harm2cl(Bharm2,Deltal,lbins)


E2 = cmbtools.harm2map(Eharm2,Delta)
B2 = cmbtools.harm2map(Bharm2,Delta)


# Plot E/B
figure()
subplot(1,2,1)
title('E')
imshow(E,interpolation='nearest')
colorbar()
subplot(1,2,2)
imshow(B,interpolation='nearest')
title('B')
colorbar()

# Plot Q/U
figure()
subplot(1,2,1)
title('Q')
imshow(Q,interpolation='nearest')
colorbar()
subplot(1,2,2)
title('U')
imshow(U,interpolation='nearest')
colorbar()

# Plot Q/U harmonics
figure()
subplot(2,2,1)
imshow(real(Qharm2), interpolation='nearest')
colorbar()
title('Qharm2')
subplot(2,2,2)
title('Uharm2')
imshow(real(Uharm2), interpolation='nearest')
colorbar()
subplot(2,2,3)
imshow(imag(Qharm2), interpolation='nearest')
colorbar()
subplot(2,2,4)
imshow(imag(Uharm2), interpolation='nearest')
colorbar()


# Plot output E/B harmonics
figure()
subplot(2,2,1)
imshow(real(Eharm2), interpolation='nearest')
colorbar()
title('Eharm2')
subplot(2,2,2)
title('Bharm2')
imshow(real(Bharm2), interpolation='nearest')
colorbar()
subplot(2,2,3)
imshow(imag(Eharm2), interpolation='nearest')
colorbar()
subplot(2,2,4)
imshow(imag(Bharm2), interpolation='nearest')
colorbar()

# Plot output E/B
figure()
subplot(1,2,1)
title('E2')
imshow(E2,interpolation='nearest')
colorbar()
subplot(1,2,2)
title('B2')
imshow(B2,interpolation='nearest')
colorbar()


# Plot input/output spectrum
figure()
plot(lcent,ClEE,label='EE')
plot(lcent,ClEEout/fmask,':',lw=6,label='EEout')
plot(lcent,ClBB,label='BB')
plot(lcent,ClBBout/fmask,':',lw=6,label='BBout')
legend()

show()
'''
