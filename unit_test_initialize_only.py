#from unit_tests import *
#from unit_tests_2 import *
from starter import *
from yt.funcs import mylog
reload(enzo_write)
reload(p49_eigen)
reload(p49_plot_tools)
mylog.setLevel(40)
def ic_and_output_one(data_dir="./Runs",plot_dir="./Plots",frame=0):
    
    ds_fname = "%s/DD%04d/data%04d"%(data_dir,frame,frame)
    prefix="RunTest"
    ds = yt.load(ds_fname)
    read_stuff = p49_eigen.get_cubes_cg(ds)
    ic_stuff = p49_plot_tools.read_initial_conditions(data_dir)
    return read_stuff, ic_stuff
read_stuff,ic_stuff=ic_and_output_one(data_dir="./Run2")
exec(open('tmp.py').read(),globals(),locals())
