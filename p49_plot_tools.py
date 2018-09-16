from go_lite_pyt3 import *
import matplotlib.colors as colors
import p49_eigen
reload(p49_eigen)
import fourier_tools_py3.fourier_filter as Filter
reload(Filter)
from p49_print_tools import *
##
field_list_prim = ['d','vx','vy','vz','hx','hy','hz','p']
field_list_enzo = ['d','vx','vy','vz','hx','hy','hz','e']
field_list_cons = ['d','px','py','pz','hx','hy','hz','e']

wave_list  =['f-', 'a-','s-','c','f+','a+','s+']
debug = 1
def chomp(directory, HydroMethod = 4, eigen='rb96'):
    map_to_label ={'d':'density','vx':'x-velocity','vy':'y-velocity',
                   'vz':'z-velocity',
                   'hx':'Bx','hy':'By','hz':'Bz','p':'GasPressure'}
    if eigen == 'rj95':
        map_to_label ={'d':'density','vx':'x-velocity','vy':'y-velocity',
                       'vz':'z-velocity',
                       'hx':'Bx','hy':'By','hz':'Bz','e':'TotalEnergy'}
    stuff = {'cubes':p49_eigen.fieldthing(),'means':p49_eigen.fieldthing()}
    for field in map_to_label:
        setname = "%s_16.h5"%map_to_label[field]
        fname = "%s/%s"%(directory,setname)
        if debug > 0:
            print("read %s"%fname)
        if not os.path.exists(fname):
            print("oh no i'm going to fail. No file %s"%fname)
        fptr = h5py.File(fname,'r')
        cube = fptr[setname][:]
        if HydroMethod == 6 and field in ['hx','hy','hz']:
            cube = cube[cubic_slice]
        else:
            cubic_slice = [slice(cube.shape[0]),slice(cube.shape[1]),slice(cube.shape[2])]

            
        stuff['cubes'][field]=cube.swapaxes(0,2)
        stuff['means'][field] = np.mean(stuff['cubes'][field])
    these_means = stuff['means']
    these_ffts  = p49_eigen.get_ffts(stuff['cubes'], these_means)
    kall,mwut=p49_eigen.rotate_back(these_ffts, these_means)
    stuff['kall']=kall
    stuff['wut']=mwut
    return stuff

def print_wave_content(mwut=None,stuff=None):

    if mwut is None:
        these_means = stuff['means']
        these_ffts  = p49_eigen.get_ffts(stuff['cubes'], these_means)
        kall,mwut=p49_eigen.rotate_back(these_ffts, these_means)
    plt.clf()
    ymax=0
    for wave in wave_list:
        this_fft = mwut.wave_content[wave]
        power = this_fft*this_fft.conj()
        ff = Filter.FourierFilter(power)
        power_1d = np.array([power[ff.get_shell(bin)].sum() for bin in range(ff.nx)])
        k2 = np.array([ff.get_shell(bin).sum() for bin in range(ff.nx)])
        power_1d=power_1d/k2
        ymax=max([ymax,power_1d.max()])
        #print("ampl wave %3s this_max %0.2e all maxx %0.2e"%(wave, power_1d.max(), ymax))

        vvv = power_1d.max()
        vvv = mwut.wave_content[wave].max()
        vvv = np.abs(mwut.wave_content[wave]).max()
        if vvv > 1e-12:
            max_str = "%5.2e"%vvv
        else:
            max_str = "%5s"%"-"
        print("ampl wave %3s this_max %s "%(wave, max_str))
##
def plot_wave_mag(mwut=None,stuff=None,output_name="thig.png", nbins=None):

    real=True
    if mwut is None:
        these_means = stuff['means']
        these_ffts  = p49_eigen.get_ffts(stuff['cubes'], these_means, real=real)
        kall,mwut=p49_eigen.rotate_back(these_ffts, these_means,real=real)
    kmag = np.sqrt(kall[0,...]**2+kall[1,...]**2+kall[2,...]**2)
    if nbins is None:
        nbins = kmag.shape[0]
    dbin = (kmag.max()-kmag.min())/nbins
    bins = np.linspace(kmag.min()-0.5*dbin,kmag.max()+0.5*dbin,nbins+1)
    bin_center = 0.5*(bins[1:]+bins[:-1])
    dig=np.digitize(kmag,bins)
    ymax=0
    

    plt.clf()
    p1d={}
    p1d['dig']=dig
    p1d['bins']=bins
    p1d['bin_center']=bin_center
    p1d['kmag']=kmag
    
    print("==== wave power ====")
    for wave in  ['f-', 'a-','s-','c','f+','a+','s+']:
        this_fft = mwut.wave_content[wave]
        p1d['ug']=this_fft
        power_complex = (this_fft*this_fft.conj())
        power=power_complex.real
        power_1d = np.array([ power[ dig==this_bin].sum() for this_bin in range(1,nbins+1)])
        k2 = np.array([ (dig==this_bin).sum() for this_bin in range(1,nbins+1)])
        power_1d=power_1d/k2
        vvv=power_1d.max()
        if vvv > 1e-12:
            max_str = "%5.2e"%vvv
        else:
            max_str = "%5s"%"-"
        print("ampl wave %3s this_max %s"%(wave, max_str))
        plt.plot(bin_center,power_1d,label=wave)
        p1d[wave]=power_1d
        nz1 = nonzero(power_1d)
        ymax=max([ymax,power_1d.max()])
    #plt.ylim([-1e-10,ymax])
    plt.yscale('log')
    plt.legend(loc=0)
    plt.savefig(output_name)
    print(output_name)
    return p1d
def plot_wave_mag_fourtool(mwut=None,stuff=None,output_name="thig.png"):
    #kinda busted

    real=True
    if mwut is None:
        these_means = stuff['means']
        these_ffts  = p49_eigen.get_ffts(stuff['cubes'], these_means, real=real)
        kall,mwut=p49_eigen.rotate_back(these_ffts, these_means,real=real)
    kmag = np.sqrt(kall[0,...]**2+kall[1,...]**2+kall[2,...]**2)
    plt.clf()
    ymax=0
    p1d={}
    for wave in ['f+']:# ['f-', 'a-','s-','c','f+','a+','s+']:
        this_fft = mwut.wave_content[wave]
        #this_fft = np.ones_like(this_fft)
        power = (this_fft*this_fft.conj()).real
        ff = Filter.FourierFilter(power)
        print(ff.nx, kall[0,...])
        power_1d = np.array([power[ff.get_shell(bin)].sum() for bin in range(ff.nx)])
        k2 = np.array([ff.get_shell(bin).sum() for bin in range(ff.nx)])
        power_1d=power_1d/k2
        ymax=max([ymax,power_1d.max()])
        #print("ampl wave %3s this_max %0.2e all maxx %0.2e"%(wave, power_1d.max(), ymax))
        print("ampl wave %3s this_max %0.2e"%(wave, power_1d.max()))
        plt.plot(range(ff.nx), power_1d,label=wave)
        p1d[wave]=power_1d
        nz1 = nonzero(power_1d)
        print("NZ1 %3s %s"%(wave,str(nz1.size)))
    plt.yscale('log')
    plt.ylim([1e-10,ymax])
    plt.legend(loc=0)
    plt.savefig(output_name)
    print(output_name)
    #plt.clf()
    #plt.plot(kall,marker='*')
    #plt.savefig('moo.png')
    p1d['kmag']=kmag
    return p1d


def plot_k_rad(wut=None, stuff=None,prefix="THING", this_format='png'):
    if wut is None:
        these_means = stuff['means']
        these_ffts  = p49_eigen.get_ffts(stuff['cubes'], these_means)
        kall,wut=p49_eigen.rotate_back(these_ffts, these_means)
    for wave in ['f-', 'a-','s-','c','f+','a+','s+']:
        fig = plt.figure(figsize=(8,8)) # Notice the equal aspect ratio
        fig.suptitle('%s %s'%(prefix,wave))
        #ax = [fig.add_subplot(1,1,i+1) for i in range(6)]
        ax = [fig.add_subplot(1,1,1,projection='polar')]

        for a in ax:
            a.set_xticklabels([])
            a.set_yticklabels([])
            a.set_aspect('equal')


        this_fft = wut.wave_content[wave]
        all_angle = np.angle(this_fft)
        flag = np.abs(this_fft) > 1e-9
        this_kmag = kmag[flag]
        this_angle = all_angle[flag]
        oname = '%s_%s_rtheta.%s'%(prefix, wave, this_format)
        ax[0].scatter(this_angle, this_kmag)
        for a in ax:
            a.set_rmax(16)
        fig.savefig(oname)
        print(oname)

def plot_var_fields(stuff=None,prefix="A_PLOT",this_format='png',
                   do = 'slice'):
    field_list=field_list_prim
    Nx = len('xyz')
    Xinch = Nx*4
    Ny = len(field_list)
    Yinch = Ny*4

    fig = plt.figure(figsize=(Xinch,Yinch)) 
    ax = [fig.add_subplot(Ny,Nx,i+1) for i in range(Ny*Nx)]

    for a in ax:
        a.set_xticklabels([])
        a.set_yticklabels([])
        a.set_aspect('equal')

    fig.subplots_adjust(wspace=0, hspace=0)
    this_sym = colors.SymLogNorm(linthresh=1e-5, linscale=0.1, 
                                 vmin=-1e-3, vmax=1e-3)
    img_args={ 'origin':'lower','interpolation':'nearest'}#,\
            #'norm': this_sym}
    oots=[None]*len(ax)

    for nf,field in enumerate(field_list):
        this_field = stuff['cubes'][field]
        nplt = nar([0,1,2])+nf*3
        for dim in [0,1,2]:
            if do == 'std_proj':
                this_proj=np.std(this_field,axis=dim)
            if do == 'slice':
                this_slice = [slice(None),slice(None),slice(None)]
                this_slice[dim]=this_field.shape[dim]//2
                my_dims = list(this_field.shape)
                my_dims.pop(dim)
                this_proj=this_field[this_slice]
                this_proj.shape=my_dims
            #print("var %3s ax %d range(%0.2e %0.2e)"%(field,dim,
            #                        this_proj.min(),this_proj.max()))
            oots[nf] = ax[nplt[dim]].imshow(this_proj,**img_args)
        ax[nplt[0]].set_ylabel(field)
    outname = "%s_%s.%s"%(prefix,do,this_format)
    fig.savefig(outname)
    if debug > 0:
        print("Wrote %s"%outname)

def default_proj_func(array,dim):
    return np.sum(np.abs(array),axis=dim)
default_imshow={ 'origin':'lower','interpolation':'nearest'}
def six_plot(array, oname='moo.png', title='', img_args=default_imshow, func=default_proj_func):
    fig = plt.figure(figsize=(12,8)) # Notice the equal aspect ratio
    fig.suptitle(title)
    ax = [fig.add_subplot(2,3,i+1) for i in range(6)]

    for a in ax:
        a.set_xticklabels([])
        a.set_yticklabels([])
        a.set_aspect('equal')

    fig.subplots_adjust(wspace=0, hspace=0)

    oots=[None]*len(ax)
    px = array.shape[0]
    def labs(tax,ii,jj,kk):
        tax.text(px-3,0,ii)
        tax.text(px-2,0,kk)
        tax.text(px-2,1,jj)
    oots[0] = ax[0].imshow(func(array.real,dim=0),**img_args)
    ax[0].set_ylabel('Real')
    labs(ax[0],'z','y','x')
    oots[1] = ax[1].imshow(func(array.real,dim=1),**img_args)
    labs(ax[1],'z','x','y')
    oots[2] = ax[2].imshow(func(array.real,dim=2),**img_args)
    labs(ax[2],'x','y','z')
    oots[3] = ax[3].imshow(func(array.imag,dim=0),**img_args)
    ax[3].set_ylabel('Imag')
    oots[4] = ax[4].imshow(func(array.imag,dim=1),**img_args)
    oots[5] = ax[5].imshow(func(array.imag,dim=2),**img_args)
    cb=fig.colorbar(oots[5],ax=ax)
    cb.cmap.set_under('w')
    #for pl in oots:
    #    pl.set_clim(-1e-3,1e-3)
    #fig.colorbar(pl)
    #fig.subplots_adjust(wspace=0, hspace=0)
    fig.savefig(oname)
    if debug > 0:
        print('wrote %s'%oname)

def plot_k_proj(wut=None,stuff=None,prefix="THING", this_format='png', func=default_proj_func):
    if wut is None:
        these_means = stuff['means']
        these_ffts  = p49_eigen.get_ffts(stuff['cubes'], these_means)
        kall,wut=p49_eigen.rotate_back(these_ffts, these_means)
    setup = False
    for wave in ['f-', 'a-','s-','c','f+','a+','s+']:
        cube =wut.wave_content[wave] 
        for cube_part in [cube.real,cube.imag]:
            for dim in [0,1,2]:
                this_fft = func(cube_part,dim=dim)
                nonzero = np.logical_and(np.abs(this_fft)>0, np.isnan(this_fft)==False)
                nonzero = np.logical_and(nonzero, np.isinf(this_fft)==False)
                if nonzero.sum() < 1:
                    continue
                this_fft=this_fft[nonzero]
                if not setup:
                    all_max=max([this_fft.real.max()])
                    nonzero_min = min([this_fft.real.min()])
                    setup=True
                else:
                    all_max=max([all_max,this_fft.real.max()])
                    nonzero_min = min([nonzero_min,this_fft.real.min()])
                #print("range ",all_max, nonzero_min)

    for wave in ['f-', 'a-','s-','c','f+','a+','s+']:
        this_fft = wut.wave_content[wave]
        this_title = '%s %s'%(prefix,wave)
        img_args={ 'origin':'lower','interpolation':'nearest'}
        img_args['norm'] = colors.Normalize(vmin=nonzero_min,vmax=all_max)
                #'norm': colors.SymLogNorm(linthresh=1e-5, linscale=0.1, vmin=-1e-3, vmax=1e-3)}
                #'norm': colors.SymLogNorm(linthresh=1e-5, linscale=0.1, vmin=-1e-3, vmax=1e-3)}
        oname = '%s_%s_yz.%s'%(prefix, wave, this_format)
        six_plot(this_fft,oname=oname,title=this_title,img_args=img_args, func=func)

def do_stuff(print_wave=False,plot_fields=False,k_mag=False,k_proj=False,
             outdir=".", stuff=None,k_func=default_proj_func):

    od=outdir
    if print_wave:
        print("=== wave content ===")
        print_wave_content(stuff=stuff)

    if plot_fields:
        print("=== plot fields ===")
        prefix = "%s%s"%(od,"fields")

        plot_var_fields(stuff=stuff, prefix=prefix, do='slice')

    if k_proj:
        print("=== k proj === ")
        prefix = "%s%s"%(od,"K")
        plot_k_proj(stuff=stuff, prefix=prefix, func=k_func)

    if k_mag:
        print("=== k mag ===")
        prefix = "%sMag.png"%od
        p1d=plot_wave_mag(stuff=stuff,output_name=prefix)
        #print_wave_content(stuff=stuff)

