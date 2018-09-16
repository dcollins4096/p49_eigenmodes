
from go_lite_pyt3 import *
import enzo_write
reload(enzo_write)
import p49_eigen
reload(p49_eigen)
import p49_plot_tools
reload(p49_plot_tools)
from p49_print_tools import *


def prnt1(thing):
    for field in ['vx','vy','vz','hx','hy','hz']:
        print("%4s %5.2f"%(field,thing[field]))
def prnt2(thing):
    for field in ['vx','vy','vz','hx','hy','hz']:
        print("%4s %5.2f %5.2f"%(field,thing[field][0],thing[field][1]))
def prnt1b(thing):
    for field in ['d','p','px','py','pz','hx','hy','hz']:
        print("%4s %5.2f"%(field,thing[field]))
def prnt2b(thing):
    for field in ['d','p','vx','vy','vz','hx','hy','hz']:
        print("%4s %5.2f %5.2f"%(field,thing[field][0],thing[field][1]))

reload(p49_eigen)

if 1:
    #r701 r701_rb_fft_f-
    directory = '/Users/dcollins/scratch/Paper49b_play/Eigen/r701_rb96_fft_f-'
    name = 'r701'
    wave='f-'
    ts701 = p49_eigen.waves(hx=1.0,hy=1.41421,hz=0.5,p=0.6,this_wave=wave, form='rb96')
    #ts701 = p49_eigen.waves(hx=1.0,hy=1.0,hz=1.0,p=0.6,this_wave=wave, form='rb96')
    if 0:
        #straight forward, works.
        k_test = nar([[1.,0.],[0.,0.],[0.,1]])
    if 1:
        k_test = nar([[1.,0.],[1.,0.],[1.,1]])
    #print('kludge; using questionable vector')
    #k_test = nar([[3.,3.],[2.,5.],[0.,1]])
    ratio = ts701.speeds['cf']/ts701.speeds['aa']
    ampl = nar([1e-6*ratio,2.3])
    ampl = nar([1,0.2])
    ampl = nar([1,0.])
    ts701.rot_write(pert_shape='fft',base_size=nar([16]*3),pert=ampl,directory=directory,
                          wave=wave,k_rot=k_test,write=False)
    ret701 = p49_eigen.waves(form=ts701.form,**ts701.quan)

if 0:
    #eigen vector check
    #ret701.check_orthonormality(k_test,ts701.rot)
    ts701.check_orthonormality()
if 0:
    #rotation check; dummy is dumb, rotten rotates, stillrotten is dummy back again.
    dummy={'vx':0.1,'vy':0.2,'vz':0.3,'hx':101.,'hy':1e3,'hz':1e4,
           'px':0.1,'py':0.2,'pz':0.3,'d':1.}
    rotten = ts701.rotate_to_xyz(dummy)
    rotten['d']=1.
    stillrotten = ts701.rotate_to_abc(rotten)
    stillrotten['d']=1
    print("=== original ===")
    prnt1(dummy)
    print("=== two rotated ===")
    prnt2(rotten)
    print("=== rotated back ===")
    prnt2(stillrotten)
 
if 0:
    #rot is the eigen vectors along k_test, rotated to physics.
    #Successfully extracts just the f- wave.
    ret701 = p49_eigen.waves(form=ts701.form,**ts701.quan)
    ret701.fields_to_wave_frame(k_test,ts701.rot)
    
    print("=== eigen: truth ===")
    prnt1b(ts701.right[wave])
    print("=== real space: along k ===")
    prnt2b(ts701.rot)
    print("=== wave space: these should match ===")
    prnt2b(ret701.wave_frame)


    #ret701.project_to_waves(k_test, ts701.rot)
    zero = {}
    for field in ts701.rot.stuff:
        zero[field] = 0.
    ret701.project_to_waves(k_test, ts701.rot, means=zero)
    print("=== wave content: direct (only f-) ===")
    wave_content_plain = ret701.to_waves(ret701.right['f-'],ret701.left)
    for wave in wave_content_plain:
        print("%5s "%wave, wave_content_plain[wave])

    print("=== wave content: ===")
    for wave in ret701.wave_content:
        print("%5s "%wave, ret701.wave_content[wave])
    # 1.) start with rotated field; t701.rot, real space.
    # 2.) given k, rotate back, project.



if 0:
    #Work In Progress:  Two waves.
    #Successfully extracts just the f- wave.
    ts22 = p49_eigen.waves(hx=1.0,hy=1.41421,hz=0.5,p=0.6, form='rb96', HydroMethod=4)
    ts22 = p49_eigen.waves(form=ts701.form,**ts701.quan)
    ratio = 1.0# ts22.speeds['cf']/ts22.speeds['aa']

    ampl = nar([1e-6*ratio,0])
    size = nar([4]*3)
    ts22.rot_write(pert_shape='fft',base_size=size,pert=ampl,directory=directory,
                          wave='f-',k_rot=nar([[1.,0.],[0.,0.],[0.,1]]),start=True)
    #ts22.rot_write(pert_shape='fft',base_size=size,pert=ampl,directory=directory,
    #                      wave='s+',k_rot=nar([[3,4],[2,4],[0,1]]),start=False)
    ts22.fields_to_wave_frame(k_test,ts701.rot)
    
    print("=== eigen: truth ===")
    prnt1b(ts701.right[wave])
    print("=== real space: along k ===")
    prnt2b(ts701.rot)
    print("=== wave space: these should match ===")
    prnt2b(ret701.wave_frame)


    #ret701.project_to_waves(k_test, ts701.rot)
    zero = {}
    for frame in ts701.rot:
        zero[frame] = 0.
    ret701.project_to_waves(k_test, ts701.rot, means=zero)
    print("=== wave content: ===")
    for wave in ret701.wave_content:
        print("%5s "%wave, ret701.wave_content[wave])
    # 1.) start with rotated field; t701.rot, real space.
    # 2.) given k, rotate back, project.


if 0:
    field_list = ['d','vx','vy','vz','hx','hy','hz','p']
    def nz(field):
        nz = np.abs(field) > 1e-12
        return field[nz]
    def wnz(field):
        nz = np.abs(field) > 1e-12
        return np.where(nz)


    if 0:
        #really stupid test.
        field_list = ['d','vx','vy','vz','hx','hy','hz','p']
        these_cubes={}
        for field in field_list:
            these_cubes
        fobase_k = np.zeros([4,4,4])*1j
    if 0:
        #great test. Works.  Two waves.
        wave='f-'
        tsfft = p49_eigen.waves(hx=1.0,hy=1.0,hz=1.0,p=0.6,this_wave=wave, form='rb96')
        k_test = nar([[1.,0.],[0.,0.],[0.,1]])
        ampl = nar([1,1.])
        nn=32
        tsfft.rot_write(pert_shape='fft',base_size=nar([nn]*3),pert=ampl,directory='',
                              wave=wave,k_rot=k_test,write=False)
        k_test = nar([[1.,1.],[0.,0.],[0.,1]])
        tsfft.rot_write(pert_shape='fft',base_size=nar([nn]*3),pert=ampl,directory='',
                              wave='a+',k_rot=k_test,write=False)
        these_ffts = p49_eigen.get_ffts(tsfft.temp_cubes, tsfft.temp_means)
        these_means = tsfft.temp_means
    if 1:
        #More modes.  IN THE WORKS.
        wave='f-'
        write=True
        name = 'rB01'
        directory = '/Users/dcollins/scratch/Paper49b_play/Eigen/rB01_rb_several'
        tsfft = p49_eigen.waves(hx=1.0,hy=1.41421,hz=0.5,p=0.6,this_wave=wave, form='rb96', HydroMethod=4)
        k_test = nar([[2.,1.],[3.,0.],[0.,1]])
        ratio = 1.0# ts.speeds['cf']/ts.speeds['aa']
        ampl = nar([1e-6*ratio,0])
        #ts.rot_write(pert_shape='fft',base_size=nar([16]*3),pert=ampl,directory=directory,
        #                      wave=wave,k_rot=k_test)
        tsfft.rot_write(pert_shape='fft',base_size=nar([32]*3),pert=ampl,directory=directory,
                              wave='s+',k_rot=nar([[1,4],[0,4],[0,1]]), write=write)
        these_ffts = p49_eigen.get_ffts(tsfft.temp_cubes, tsfft.temp_means)
        these_means = tsfft.temp_means
    if 0:
        #also works.
        directory = '/Users/dcollins/scratch/Paper49b_play/Eigen/y701_rb96_fft_f-_play'
        name = 'y701'
        wave='f-'
        #tsy701 = p49_eigen.waves(hx=1.0,hy=1.0,hz=1.0,p=0.6,this_wave=wave, form='rb96')
        tsy701 = p49_eigen.waves(hx=1.0,hy=1.41421,hz=0.5,p=0.6,this_wave=wave, form='rb96')
        k_test = nar([[1.,1.],[0.,0.],[0.,1]])
        ratio = tsy701.speeds['cf']/tsy701.speeds['aa']
        ampl = nar([1e-6*ratio,0])
        tsy701.rot_write(pert_shape='fft',base_size=nar([16]*3),pert=ampl,directory=directory,
                              wave=wave,k_rot=k_test, write=True)
        these_ffts = p49_eigen.get_ffts(tsy701.temp_cubes, tsy701.temp_means)
        these_means = tsy701.temp_means
    if 0:
        #yay! Load the data, k=100
        frame = 0
        directory = '/Users/dcollins/scratch/Paper49b_play/Eigen/y701_rb96_fft_f-_play'
        ds = yt.load("%s/DD%04d/data%04d"%(directory,frame,frame))
        stuff = p49_eigen.get_cubes_cg(ds)
        these_means = stuff['means']
        these_ffts  = p49_eigen.get_ffts(stuff['cubes'], these_means)
    if 0:
        #harder run; k=110.  Broken
        wave='f-'
        directory = '/Users/dcollins/scratch/Paper49b_play/Eigen/r801_rj95_110_f-'
        ts801 = p49_eigen.waves(hx=1.0,hy=1.41421,hz=0.5,p=0.6,this_wave=wave, form='rb96')
        k_test = nar([[1.,1.],[1.,0.],[0.,0]])
        ratio = ts801.speeds['cf']/ts801.speeds['aa']
        ampl= nar([1e-6*ratio,0])
        #ampl = nar([1.,0.])
        ts801.rot_write(pert_shape='fft',base_size=nar([16]*3),pert=ampl,directory=directory,
                              wave=wave,k_rot=k_test,write=True)
        tsfft = ts801
        these_ffts = p49_eigen.get_ffts(tsfft.temp_cubes, tsfft.temp_means)
        these_means = tsfft.temp_means
        ###
    if 0:
        #harder test: read rotated
        frame = 0
        directory = '/Users/dcollins/scratch/Paper49b_play/Eigen/r801_rj95_110_f-'
        ds = yt.load("%s/DD%04d/data%04d"%(directory,frame,frame))
        stuff = p49_eigen.get_cubes_cg(ds)
        these_means = stuff['means']
        these_ffts  = p49_eigen.get_ffts(stuff['cubes'], these_means)
        ###
    if 0:
        #Rotated, rb96.  Works by itself.
        #rA01 
        directory = '/Users/dcollins/scratch/Paper49b_play/Eigen/rA01_rb96_110_f-'
        name = 'rA01'
        wave='f-'
        tsA01 = p49_eigen.waves(hx=1.0,hy=1.41421,hz=0.5,p=0.6,this_wave=wave, form='rb96',\
                                HydroMethod = 4)
        k_test = nar([[1.,1.],[1.,0.],[0.,0]])
        ratio =  tsA01.speeds['cf']/tsA01.speeds['aa']
        ampl = nar([1e-6*ratio,0])
        tsA01.rot_write(pert_shape='fft',base_size=nar([16]*3),pert=ampl,directory=directory,
                              wave=wave,k_rot=k_test, write=True)
        these_ffts = p49_eigen.get_ffts(tsA01.temp_cubes, tsA01.temp_means)
        these_means = tsA01.temp_means
        tsfft = tsA01
    if 0:
        #Rotate, rb96, life: IN PROCESS.
        frame = 50
        directory = '/Users/dcollins/scratch/Paper49b_play/Eigen/rA01_rb96_110_f-'
        ds = yt.load("%s/DD%04d/data%04d"%(directory,frame,frame))
        stuff = p49_eigen.get_cubes_cg(ds)
        these_means = stuff['means']
        these_ffts  = p49_eigen.get_ffts(stuff['cubes'], these_means)

    print_fields = False
    print_waves = True
    kall,wut=p49_eigen.rotate_back(these_ffts, these_means)
    fl =  np.zeros_like(wut.wave_frame['d']).astype('bool')
    if print_fields:
        for field in  wut.wave_frame:
            print(" ===== %s ===="%field)
            thisthing =  wut.wave_frame[field]
            thisthing =  wut.dumb[field]
            #print("  eigen    %s"%str(tsfft.right['f-'][field]))
            #print("  rot      %s"%str(tsfft.rot[field]))
            #print("all_hat  %3s %s"%(field, nz(tsfft.all_hats[field])))
            #aaa = these_ffts[field] #is good.
            #print("also fft input  k %3s %s"%(field, str(nz(aaa).size)))
            print("this wave frame n k %3s %s"%(field, str(nz(thisthing).size)))
            print("this wave frame   a %3s %s"%(field, str(nz(thisthing))))
            print("this wave frame   k %3s %s"%(field, str(wnz(thisthing))))
    if print_waves:
        for wave in wut.wave_content:
            thisthing = wut.wave_content[wave]
            bang_or_not = ""
            if ( np.abs(thisthing)>1e-12).sum() > 0:
                bang_or_not = "!!!"*8 + " meann %0.2e max %0.2e"%(np.mean(np.abs(thisthing)),np.abs(thisthing).max())
            print("=== Wave %s %s"%(wave, bang_or_not))
            s1 = str(nz(thisthing.real).size)
            s2 = str(nz(thisthing.imag).size)
            print("wave real nz %s imag %s"%(s1,s2))
        

