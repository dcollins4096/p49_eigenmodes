
from starter import *
import unit_tests
def test_k_multiple(state=None,WRITE=False, directory ='./Run/' ):
    if state is None:
        state = p49_eigen.waves(hx=1.0,hy=1.41421,hz=0.5,p=0.6,this_wave="f-", HydroMethod=4)
    size = 16
    ks = p49_eigen.make_k_freqs_and_int(size)
    k3 = ks['k_freq']
    kint = ks['k_int']
    k_norm    = np.sqrt( k3[0,...]**2 + k3[1,...]**2 + k3[2,...]**2 )
    ampl  = np.zeros_like(k_norm)*(1+0j)
    theta = np.pi*2*np.random.random(ampl.shape)
    phase = np.exp(theta*1j)
    mask  = np.ones_like(kint[2,...],dtype='bool')
    mask  = np.logical_and(mask, kint[2,...]==1)
    mask  = np.logical_and(mask, kint[1,...]==1)
    mask  = np.logical_and(mask, kint[0,...]==1)
    ampl[mask] = k_norm[mask]**(-5./3)*phase[mask]
    ampl[mask] *= 1e-3/ampl[mask].max()
    
    wave='f-'
    state.perturb(pert_shape='fft',base_size=nar([size]*3),
                  ampl=ampl,directory=directory,
                  k_rot=kint,
                      wave=wave, start=True,write=WRITE, single_wave_selector=False)
    mask  = np.ones_like(kint[2,...],dtype='bool')
    mask  = np.logical_and(mask, kint[2,...]==1)
    mask  = np.logical_and(mask, kint[1,...]==2)
    mask  = np.logical_and(mask, kint[0,...]==3)
    ampl[mask] = k_norm[mask]**(-5./3)*phase[mask]
    ampl[mask] *= 1e-3/ampl[mask].max()
    
    wave='s+'
    state.perturb(pert_shape='fft',base_size=nar([size]*3),
                  ampl=ampl,directory=directory,
                  k_rot=kint,
                      wave=wave, start=False,write=WRITE, single_wave_selector=False)
    these_ffts  = p49_eigen.get_ffts(state.cubes, means={})
    return state, these_ffts
#state, mask=test_k_multiple(directory = './Run3', WRITE=True)
ic_stuff=unit_tests.test_read_ics(ic_dir="./Run3",plot_dir="./Vis3/")
for frame in [0,1]:
    read_stuff=unit_tests.test_read_outputs(data_dir="./Run3",plot_dir="./Vis3/",frame=frame)
