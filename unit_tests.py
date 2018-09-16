
from go_lite_pyt3 import *
import enzo_write
reload(enzo_write)
import p49_eigen
import p49_plot_tools
reload(p49_plot_tools)
from p49_print_tools import *

class miscellaneous_things():
    def __init__(self):
        self.wave_list=['f-', 'a-','s-','c','f+','a+','s+']
        self.field_list = ['d','px','py','pz','hx','hy','hz','p']
        self.form_s = " %5s"
        self.form_f = " %5.2f"

def test_mean_state(state=None):
    """one vector along \hat{x}, one along \hat{z}, left-going fast wave."""
    
    reload(p49_eigen)
    t=miscellaneous_things()

    #The mean background state.  rho=1, \vec{v}=0,0,0, gamma=5./3
    #This populates a set of left- and right-eigen vectors, and
    #the average state, *quan*
    if state is None:
        state = p49_eigen.waves(hx=1.0,hy=1.41421,hz=0.5,p=0.6,this_wave="f-")

    print("=== Mean ===")
    for q in state.quan:
        print("%15s %5.2f"%(q,state.quan[q]))

    quans = list(state.right['f-'])
    
    print("=== Right Eigenvectors ===")
    print(t.form_s%" "+t.form_s*len(quans)%tuple(quans))
    for wave in t.wave_list:
        out = t.form_s%wave
        for q in quans:
            out += t.form_f%state.right[wave][q]
        print(out)
    print("=== Left Eigenvectors ===")
    print(t.form_s%" "+t.form_s*len(quans)%tuple(quans))
    for wave in t.wave_list:
        out = t.form_s%wave
        for q in quans:
            out += t.form_f%state.left[wave][q]
        print(out)
    print("=== Right_{ik} Left_{kj} = \delta_{ij} === ")
    state.check_orthonormality()
    return state

#state = test_mean_state()
def test_two(state=None):

    if state is None:
        state = p49_eigen.waves(hx=1.0,hy=1.41421,hz=0.5,p=0.6,this_wave="f-")
    t=miscellaneous_things()
    reload(p49_eigen)
    #two waves, one in \hat{x}, one in \hat{z}
    k_test = nar([[1.,0.],[0.,0.],[0.,1]])

    #amplitudes for each of the waves.
    #This ratio makes things consistent between two eigen systems, ignore it.
    ratio = state.speeds['cf']/state.speeds['aa'] 
    ampl = nar([1e-6*ratio,2.3]) 

    #do the perturbation. 
    #write the result in *directory*
    state.perturb(pert_shape='fft',base_size=nar([16]*3),ampl=ampl,directory=".",
                          wave="f-",k_rot=k_test,write=False,single_wave_selector=True)

    
    for q in t.field_list:
        s = state.rot[q]
        print(t.form_s%q + t.form_f*len(s)%tuple(s))
    return state
state = test_two()



