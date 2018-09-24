
from go_lite_pyt3 import *
import yt
import enzo_write
reload(enzo_write)
import p49_eigen
reload(p49_eigen)
import p49_plot_tools
reload(p49_plot_tools)
from p49_print_tools import *
import unit_tests
reload(unit_tests)
class miscellaneous_things():
    def __init__(self):
        self.wave_list=['f-', 'a-','s-','c','f+','a+','s+']
        self.field_list = ['d','px','py','pz','hx','hy','hz','p']
        self.form_s = " %5s"
        self.form_f = " %5.2f"
t=miscellaneous_things()
def test_k_array_simple(state=None,WRITE=False, directory ='./Run/' ):
    if state is None:
        state = p49_eigen.waves(hx=1.0,hy=1.41421,hz=0.5,p=0.6,this_wave="f-")
    size = 32
    ks = p49_eigen.make_k_freqs_and_int(size)
    k3 = ks['k_freq']
    kint = ks['k_int']
    kint_norm = np.sqrt( kint[0,...]**2 + kint[1,...]**2 + kint[2,...]**2 )
    k_norm    = np.sqrt( k3[0,...]**2 + k3[1,...]**2 + k3[2,...]**2 )
    ampl  = np.zeros_like(k_norm)*(1+0j)
    theta = np.pi*2*np.random.random(ampl.shape)
    phase = np.exp(theta*1j)
    mask  = np.ones_like(kint[2,...],dtype='bool')
    mask  = np.logical_and(mask, kint[2,...]>0)
    mask  = np.logical_and(mask, kint[1,...]>0)
    mask  = np.logical_and(mask, kint[0,...]>0)

    mask  = np.logical_and(mask, kint[0,...] == 2)
    mask  = np.logical_and(mask, kint[1,...] == 1)
    mask  = np.logical_and(mask, kint[2,...] == 1)
    ampl[mask] = k_norm[mask]**(-5./3)*phase[mask]
    ampl[mask] *= 1e-3/ampl[mask].max()
    
    wave='f-'
    state.perturb(pert_shape='fft',base_size=nar([size]*3),
                  ampl=ampl,directory=directory,
                  k_rot=kint,
                      wave=wave, start=True,write=WRITE, single_wave_selector=False)
    these_ffts  = p49_eigen.get_ffts(state.cubes, means={})
    return state, these_ffts
def test_k_array_powers(state=None,WRITE=False, directory ='./Run/' ):
    if state is None:
        state = p49_eigen.waves(hx=1.0,hy=1.41421,hz=0.5,p=0.6,this_wave="f-")
    size = 32
    ks = p49_eigen.make_k_freqs_and_int(size)
    k3 = ks['k_freq']
    kint = ks['k_int']
    kint_norm = np.sqrt( kint[0,...]**2 + kint[1,...]**2 + kint[2,...]**2 )
    k_norm    = np.sqrt( k3[0,...]**2 + k3[1,...]**2 + k3[2,...]**2 )
    ampl  = np.zeros_like(k_norm)*(1+0j)
    theta = np.pi*2*np.random.random(ampl.shape)
    phase = np.exp(theta*1j)
    mask  = np.ones_like(kint[2,...],dtype='bool')
    mask  = np.logical_and(mask, kint[2,...]>0)
    ampl[mask] = k_norm[mask]**(-5./3)*phase[mask]
    ampl[mask] *= 1e-3/ampl[mask].max()
    
    wave='f-'
    state.perturb(pert_shape='fft',base_size=nar([size]*3),
                  ampl=ampl,directory=directory,
                  k_rot=kint,
                      wave=wave, start=True,write=WRITE, single_wave_selector=False)
    these_ffts  = p49_eigen.get_ffts(state.cubes, means={})
    return state, these_ffts
state, mask=test_k_array_powers(directory = './Run2', WRITE=True)
unit_tests.test_read_ics(ic_dir="./Run2",plot_dir="./Vis2/")
#unit_tests.test_read_outputs(data_dir="./Run2",plot_dir="./Vis2")
