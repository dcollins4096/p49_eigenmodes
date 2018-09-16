
from go_lite_pyt3 import *
import enzo_write
reload(enzo_write)
import p49_eigen
reload(p49_eigen)
import p49_plot_tools
reload(p49_plot_tools)
from p49_print_tools import *
plt.close('all')

#size = 128
size = 32
sim = 'rx10'
sim = 'rx06b_s-'

blah=False
if 1:
    #setup FFT init
    blah=False
    #Target woodshed
    #k3 = p49_eigen.make_k_freqs(size)
    ks = p49_eigen.make_k_freqs_and_int(size)
    k3 = ks['k_freq']
    kint = ks['k_int']
    kint_norm = np.sqrt( kint[0,...]**2 + kint[1,...]**2 + kint[2,...]**2 )
    k_norm = np.sqrt( k3[0,...]**2 + k3[1,...]**2 + k3[2,...]**2 )
    ok = k_norm > 0


amplitude=1e-3
ampl=np.zeros_like(k_norm)*(1+0j)
old_b   = nar([1.0, 1.41421, 0.5])
sim_dir = "../p49_run"
out_prefix_analysis = "%s_"%sim
out_dir_analysis="../p49_plot/"+sim
out_name = "%s/%s"%(out_dir_analysis,out_prefix_analysis)
if sim=='rx01':
    #single wave.  Works, man.
    ampl[2,0,2]=amplitude
    out_dir_s = sim_dir+'/rx01_202'

if sim=='rx02':
    #reverse wave.  Works.
    ampl[-2,0,2]=amplitude
    out_dir_s = sim_dir+'/rx02_-202'

if sim=='rx03':
    #single wave, different angle.  Works fine.
    ampl[2,2,0]=amplitude
    out_dir_s = sim_dir+'/rx03_220'

if sim=='rx04':
    #single wave, more angle.  Generates even odd.  fascinating.
    #Seems to get corrugated in Vz and Bz first.
    ampl[2,2,1]=amplitude
    out_dir_s = sim_dir+'/rx04_221'

if sim=='rx05':
    #single wave, more angle, rotated.  Hx gets corrugated first?
    ampl[2,1,2]=amplitude
    out_dir_s = sim_dir+'/rx05_212'

if sim=='rx06':
    #single wave, more even-odd hunting.  Even-odd in y? Yes, probably.
    old_b   = nar([1.0, 0.5, 1.41421])
    ampl[2,1,2]=amplitude
    out_dir_s = sim_dir+'/rx06_212_b2'

if sim=='rx06b_s-':
    #single wave, more even-odd hunting.  Even-odd in y? Yes, probably.
    old_b   = nar([1.0, 0.5, 1.41421])
    ampl[2,1,2]=amplitude
    wave='s-'
    out_dir_s = sim_dir+'/rx06b_212_s-'

if sim=='rx07':
    if size != 128:
        print("SIZE ERROR THIS ONE NEEDS TO BE LARGE")
        print("SIZE ERROR THIS ONE NEEDS TO BE LARGE")
        print("SIZE ERROR THIS ONE NEEDS TO BE LARGE")
        print("SIZE ERROR THIS ONE NEEDS TO BE LARGE")
    #single wave, more even-odd hunting.  
    #WITH RiemannSolver=1, this is a horrid even-odd. 
    #RiemannSolver =6 works nicely.
    wave='f-'
    old_b   = nar([1.0, 0.5, 1.41421])
    ampl[2,1,2]=amplitude
    out_dir_s = sim_dir+'/rx07_212_128'
    out_dir_s = sim_dir+'/rx07b_212_128_HLLD'
    out_dir_s = sim_dir+'/rx07b_212_128_HLLD_check'
else:
    if size != 32:
        print("SIZE ERROR needs to be 32")
        print("SIZE ERROR needs to be 32")
        print("SIZE ERROR needs to be 32")
        print("SIZE ERROR needs to be 32")


if sim=='r701':
    #like r701.  Good.
    ampl=np.zeros_like(k_norm)*(1+0j)
    wave='f-'
    ts_temp = p49_eigen.waves(hx=1.0,hy=1.41421,hz=0.5,p=0.6,this_wave=wave, form='rb96')
    ratio = ts_temp.speeds['cf']/ts_temp.speeds['aa']
    ampl[1,0,0]=1e-6*ratio
    out_dir_s = sim_dir+'/r701c_32'

if sim=='rx09':
    #rx09: two waves.
    theta = np.pi*2*np.random.random(ampl.shape)
    theta = 0.
    phase = np.exp(theta*1j)
    mask=ok
    #mask = np.logical_and(mask,k3[0,...]>0)
    #mask = np.logical_and(mask,k3[1,...]>0)
    #mask = np.logical_and(mask,k3[2,...]>0)

    mask[:]=False
    mask[2,1,2]=True
    #mask[1,:,7]=True
    mask[5,1,7]=True
    phase=1.0#phase[mask]
    ampl[mask] = k_norm[mask]**(-5./3)*phase
    ampl[mask] *= 1e-3/ampl[mask].max()
    wave='s-'
    out_dir_s= sim_dir+'/rx09_two_wves_212_517_s-'


if sim=='rx08':
    #rx08: ensemble of waves. Fixed 
    theta = 0. #np.pi*2*np.random.random(ampl.shape)
    phase = np.exp(theta*1j)
    mask=ok

    phase=1.0
    mask = np.logical_and(mask, kint[2,...]>0)
    ampl[mask] = k_norm[mask]**(-5./3)*phase
    ampl[mask] *= 1e-3/ampl[mask].max()
    out_dir_s= sim_dir+'/rx08_k-53'

if sim=='rx08_seg_test':
    #this is actually rx06, enzo threw a seg-fault.  It works
    old_b   = nar([1.0, 0.5, 1.41421])
    ampl[2,1,2]=amplitude
    out_dir_s= sim_dir+'/rx08_k-53'

if sim=='rx10':
    #rx10: ensemble of waves, random phase.
    theta = np.pi*2*np.random.random(ampl.shape)
    phase = np.exp(theta*1j)
    mask=ok

    mask = np.logical_and(mask, kint[2,...]>0)
    ampl[mask] = k_norm[mask]**(-5./3)*phase[mask]
    ampl[mask] *= 1e-3/ampl[mask].max()
    out_dir_s= sim_dir+'/rx10_k-53_phase'


if 0:
    #test the constituents
    import matplotlib.colors as colors
    img_args={ 'origin':'lower','interpolation':'nearest'}
    img_args['norm'] = colors.Normalize(vmin=0,vmax=np.abs(ampl).max())
    p49_plot_tools.six_plot(ampl,oname = 'AAA/ampl.png',img_args=img_args)
    p49_plot_tools.six_plot(k_norm,oname = 'AAA/k_norm.png')
    p49_plot_tools.six_plot(kint_norm,oname = 'AAA/k_int.png')

if 0:
    ampl=np.ones_like(k_norm)*(amplitude+0j)
    this_mask = np.ones_like(ampl,dtype=bool)
    this_mask[-5,7,1]=True
    blank_mask = this_mask == False
    ampl[blank_mask]=0.


if 0:
    #works.
    directory_good = sim_dir+'/test_good'
    kint_good = nar([[2,1],[0,0],[2,0]])
    ampl_good = nar([amplitude,0])
    blah_good=True
if 0:
    #works.
    kint_good = nar([[-2,1],[2,0],[1,0]])
    ampl_good = nar([amplitude,0])
    blah_good=True

if 'wave' not in dir():
    wave = 's+'

if 'unlazy' not in dir():
    unlazy=True
if 'ts1' not in dir() or unlazy:
    WRITE=True

    if 'ampl_good' in dir():
        ts1 = p49_eigen.waves(hx=old_b[0],hy=old_b[1],hz=old_b[2],p=0.6,
                              this_wave=wave, form='rb96', HydroMethod=4)
        ts1.rot_write(pert_shape='fft',base_size=nar([size]*3),
                      pert=ampl_good,directory=directory_good,
                      k_rot=kint_good,
                      wave=wave, start=True,write=WRITE, blah=blah_good)

    ts2 = p49_eigen.waves(hx=old_b[0],hy=old_b[1],hz=old_b[2],p=0.6,
                          this_wave=wave, form='rb96', HydroMethod=4)
    print("wrtie the what?")
    ts2.rot_write(pert_shape='fft',base_size=nar([size]*3),
                  pert=ampl,directory=out_dir_s,
                  k_rot=kint,
                  wave=wave, start=True,write=WRITE, blah=blah)

in_dir = out_dir_s
sn = p49_plot_tools.chomp(in_dir)
def maxis(inarr,dim):
    arr=np.zeros_like(inarr)
    ok =  inarr>1e-13
    arr[ok] = inarr[ok]
    out = np.log(np.max(np.abs(arr),axis=dim))

    #out = np.sum(np.abs(inarr),axis=dim)
    #out[ out==0] == -1
    return out

analysis={'print_wave':True,
          'plot_fields':1,
          'k_mag':1,
          'k_func':maxis,
          'k_proj':1}
p49_plot_tools.do_stuff(stuff=sn,outdir=out_name,**analysis)
