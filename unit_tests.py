"""Here are several test functions for the eigen-vector project.
Read.  Then horse around.
Not everything works all that great.  You might need to fix stuff.

There are two formulations of the eigen system, "rj95" and "rb96".  
We'll use rb96 always.  There are minor differences.

In several places the field needs to be written with the right label, so Enzo can read it.

px,py,pz is momentum.
hx,hy,hz is magnetic field.
There's a tool called "fieldthing" that takes care of keeping momentum and density consistent.
Access is kind of like a dictionary, but I haven't fully fleshed out the object.



"""
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

def test_two(state=None, write=False):

    if state is None:
        state = p49_eigen.waves(hx=1.0,hy=1.41421,hz=0.5,p=0.6,this_wave="f-")
    t=miscellaneous_things()
    reload(p49_eigen)
    #two waves, one in \hat{x}, one in \hat{z}
    k_test = nar([[1.,0.],[1.,0.],[1.,1]])

    #amplitudes for each of the waves.
    #This ratio makes things consistent between two eigen systems, ignore it.
    ratio = state.speeds['cf']/state.speeds['aa'] 
    ampl = nar([1e-6*ratio,0]) 

    #do the perturbation. 
    #write the result in *directory*
    state.perturb(pert_shape='fft',base_size=nar([16]*3),ampl=ampl,directory="./Runs",
                          wave="f-",k_rot=k_test,write=write,single_wave_selector=True)

    

    nstates = len(state.rot['d'])
    print(t.form_s%" " + t.form_s*nstates%tuple(["k%d"%i for i in range(nstates)]))
    print("A hat (along B0 vector)")
    for d in range(3):
        print(t.form_s%"a%d"%d + t.form_f*nstates%tuple(state.a_unit[d,...]))
    print("B hat (perp. to B0 and \hat{k})")
    for d in range(3):
        print(t.form_s%"a%d"%d + t.form_f*nstates%tuple(state.b_unit[d,...]))
    print("C hat (in B0 k plane, perp to B0)")
    for d in range(3):
        print(t.form_s%"a%d"%d + t.form_f*nstates%tuple(state.c_unit[d,...]))

    print("Perturbations")
    print(t.form_s%" " + t.form_s*nstates%tuple(["k%d"%i for i in range(nstates)]))
    for q in t.field_list:
        s = state.rot[q]
        print(t.form_s%q + t.form_f*len(s)%tuple(s))
    print("Net arrays (mean, std)")
    for q in t.field_list:
        print( t.form_s%q+ p49_plot_tools.mean_min_max(state.cubes[q]))
    return state

def test_read_ics(ic_dir="./Runs", plot_dir="./Plots/"):
    """read in initial conditions from *ic_dir*.
    Make some plots in *plot_dir*
    p49_plot_tools.read_initial_conditions consumes the initial data files in *ic_dir*
    p49_plot_tools.do_stuff makes several figures, based on the arguments.
    Here we use a fun trick of double-asterisk to pass in a bunch of keyword args."""
    data = p49_plot_tools.read_initial_conditions(ic_dir)
    analysis={'print_wave':True,  #prints the wave content to the screen.
              'plot_fields':True, #makes slices (or something) of density, velocity, etc.
              'k_func':p49_plot_tools.maxis,     #The projection function for the above.
              'k_mag':True,       #Makes power spectra for each wave.
              'k_proj':True}      #makes an inscrutable plot of K-space
    p49_plot_tools.do_stuff(stuff=data,outdir=plot_dir,**analysis)
state = test_mean_state()
state = test_two(write=True)
test_read_ics()

def test_temp_fft(ic_dir="./Runs", plot_dir="./Plots/"):
    """read in initial conditions from *ic_dir*.
    Make some plots in *plot_dir*
    p49_plot_tools.chomp consumes the initial data files in *ic_dir*
    p49_plot_tools.do_stuff makes several figures, based on the arguments.
    Here we use a fun trick of double-asterisk to pass in a bunch of keyword args."""
    data = p49_plot_tools.chomp(ic_dir)
    ffts = p49_eigen.get_ffts(data['cubes'])
    return data, ffts

ics = test_temp_fft()


