from go import *
if '/home/dcollins/repos/p49c/p49_eigenmodes' not in sys.path:
    sys.path.append('/home/dcollins/repos/p49c/p49_eigenmodes')
if './p49c' not in sys.path:
    sys.path.append('./p49c')

import p49_eigen
reload(p49_eigen)
import p49_plot_tools
reload(p49_plot_tools)
import p49_QU2EB
import turb_quan
reload(turb_quan)
import plots
reload(plots)

from p49_print_tools import *
#p49_eigenmodes/unit_tests_EB.py
#import unit_tests_EB
#reload(unit_tests_EB)
def project_and_fft(array, line_of_sight):
    a_flat = np.sum(array,axis=line_of_sight)
    #a_flat = array[8,:,:]
    a_flat.shape = (array.shape[0],array.shape[1])
    a_flat=np.ascontiguousarray(a_flat)
    a_fft =  np.fft.rfftn( a_flat)
    #a_fft[a_fft < 1e-8]=0
    return a_fft
def one_wave(size, k_unit_Q=None,k_unit_U=None):
    x,y = np.mgrid[0:1:1./size, 0:1:1./size]
    if k_unit_Q is None: k_unit_Q = np.array([0,1])*2*np.pi
    Q_flat =   (np.exp(1j*(k_unit_Q[0]*x+k_unit_Q[1]*y))).imag
    if k_unit_U is None: k_unit_U = np.array([0,1])*2*np.pi
    U_flat =   (np.exp(1j*(k_unit_U[0]*x+k_unit_U[1]*y))).imag
    Q=Q_flat
    U=U_flat
    #Q = np.zeros([size,size])
    #U = np.zeros([size,size])
    #dx = 1./size
    #x_list = np.arange(0,1,dx)
    #y_list = np.arange(0,1,dx)
    #for i,x in enumerate(x_list):
    #    for j,y in enumerate(y_list):
    #        Q[i,j] = (np.exp(1j*(k_unit_Q[0]*x+k_unit_Q[1]*y))).imag
    #        U[i,j] = (np.exp(1j*(k_unit_U[0]*x+k_unit_U[1]*y))).imag
    return {'Q':Q, 'U':U}

class extents():
    def __init__(self):
        self.minmax=[]
    def __call__(self,array):
        if len(self.minmax):
            self.minmax[0] = min([array.min(),self.minmax[0]])
            self.minmax[1] = max([array.max(),self.minmax[1]])
        else:
            self.minmax=[array.min(),array.max()]
    def __getitem__(self,index):
        return self.minmax[index]

class k_maker():
    def __init__(self,size):
        self.size=size
    def setup(self,size):
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
        self.k_norm=k_norm
        self.k3=k3
        self.kint=kint
        self.ampl=ampl
        self.mask = mask
        self.phase=phase
class modes_k001(k_maker):
    def __call__(self,size):
        self.setup(size)
        self.mask  = np.logical_and(self.mask, self.kint[0,...]==0)
        self.mask  = np.logical_and(self.mask, self.kint[1,...]==0) #
        self.mask  = np.logical_and(self.mask, self.kint[2,...]==1) #
        self.ampl[self.mask] =  self.k_norm[self.mask]#**(-5./3)*phase[self.mask]
        self.ampl[self.mask] *= 1e-3/self.ampl[self.mask].max()

class modes_k011(k_maker):
    def __call__(self,size):
        self.setup(size)
        self.mask  = np.logical_and(self.mask, self.kint[0,...]==0)
        self.mask  = np.logical_and(self.mask, self.kint[1,...]==1) #
        self.mask  = np.logical_and(self.mask, self.kint[2,...]==1) #
        self.ampl[self.mask] =  self.k_norm[self.mask]#**(-5./3)*phase[self.mask]
        self.ampl[self.mask] *= 1e-3/self.ampl[self.mask].max()

class modes_k032(k_maker):
    def __call__(self,size):
        self.setup(size)
        self.mask  = np.logical_and(self.mask, self.kint[0,...]==0)
        self.mask  = np.logical_and(self.mask, self.kint[1,...]==3) #
        self.mask  = np.logical_and(self.mask, self.kint[2,...]==2) #
        self.ampl[self.mask] =  self.k_norm[self.mask]#**(-5./3)*phase[self.mask]
        self.ampl[self.mask] *= 1e-3/self.ampl[self.mask].max()



def modes(h0=1,theta=0,state=None,WRITE=False, directory ='./Runs/', wave='f-', size=32, mode_tool=modes_k001):
    if state is None:
        r2=np.sqrt(2)
        state = p49_eigen.waves(hx=h0*np.cos(theta),hy=h0*np.sin(theta)/r2,hz=h0*np.sin(theta)/r2,p=0.6,this_wave=wave, HydroMethod=4)
        #state = p49_eigen.waves(hx=h0*np.cos(theta),hy=0,hz=h0*np.sin(theta),p=0.6,this_wave=wave, HydroMethod=4)
        
    
    mode_instance = mode_tool(size)
    mode_instance(size)
    ampl, kint = mode_instance.ampl, mode_instance.kint
    state.perturb(pert_shape='fft',base_size=nar([size]*3),
                  ampl=ampl,directory=directory,
                  k_rot=kint,
                      wave=wave, start=True,write=WRITE, single_wave_selector=False)

    these_ffts  = p49_eigen.get_ffts(state.cubes, means={})
    state.kstuff={'ampl':ampl,'mask':mode_instance.mask,'k':mode_instance.k3}
    return state, these_ffts
def state_to_qu(state,data=None, line_of_sight=0):
    field_horizontal = {0:'hy',1:'hz',2:'hx'}[line_of_sight]
    field_vertical   = {0:'hz',1:'hx',2:'hy'}[line_of_sight]
    axis_vertical   = {0:'hy',1:'hx',2:'hx'}[line_of_sight]
    axis_horizontal = {0:'hz',1:'hz',2:'hy'}[line_of_sight]
    axis_labels=[axis_horizontal,axis_vertical]

    if data is None:
        data=state.cubes
    epsilon = data['d']
    B_sq = data['hx']**2.0 + data['hy']**2.0 + data['hz']**2.0
    size=B_sq.shape[0]
    #epsilon = n
    Q_local =  ( epsilon * (((data[field_horizontal])**2.0) - ((data[field_vertical])**2.0))/B_sq )
    #Q_local =  (epsilon * (((data[field_horizontal])**2.0) - ((data[field_vertical])**2.0)) )#/B_sq )
    #n = data['density']
    #B_sq = data['Bx']**2.0 + data['By']**2.0 + data['Bz']**2.0    
    #epsilon = n
    U_local = (2.0 * epsilon * ((data[field_horizontal]) * (data[field_vertical]))/B_sq)
    U_local[B_sq == 0 ] = 0
    U_flat = np.sum(U_local,axis=line_of_sight)/size
    U_flat.shape = (U_local.shape[0],U_local.shape[1])
    Q_flat = np.sum(Q_local,axis=line_of_sight)/size
    Q_flat.shape = (Q_local.shape[0],U_local.shape[1])

    dhat=project_and_fft(data['d'], line_of_sight)
    b2ha=project_and_fft(B_sq, line_of_sight)
    hhhat=project_and_fft(data[field_horizontal], line_of_sight)
    hvhat=project_and_fft(data[field_vertical], line_of_sight)
    hats = {'dhat':dhat,'b2ha':b2ha,'hhhat':hhhat,'hvhat':hvhat}

    turb_quan.plotter(np.abs(dhat), 
                      np.abs(b2ha), 
                      np.abs(hhhat), 
                      np.abs(hvhat), 
                      'p49c_plots/test_field_hat.png',
                      norm='positive', 
                      labs=['nhat','h2hat','Hhhat','Hvhat'], 
                      axis_labels=axis_labels )


    
    #print("KLUDGE QU")
    #qu = one_wave(size,k_unit_Q=np.array([0,2])*2*np.pi,k_unit_U=np.array([0,2])*2*np.pi)
    #Q_flat=qu['Q']
    #U_flat=qu['U']
    #print("KLUDGE Q")
    #y,x = np.mgrid[0:1:1./size, 0:1:1./size]
    #k_unit = nar([8.,14])*2*np.pi
    #Q_flat =   (np.exp(1j*(k_unit[0]*x+k_unit[1]*y))).imag
    output = {'Q':Q_flat,'U':U_flat}
    output['proj']=[dhat,b2ha,hhhat,hvhat]
    output['Q_local']=Q_local
    output['U_local']=U_local
    output['n']=epsilon
    output['H2']=B_sq
    output['Hh']=data[field_horizontal]
    output['Hv']=data[field_vertical]
    output['axis_labels']=axis_labels
    output['h_pos']=[field_horizontal,field_vertical]
    output['hats']=hats
    return output
def wave_to_EB(wave, theta=np.pi/4,size=16,mode_tool=modes_k001):
    state, fft = modes(theta=theta, wave=wave,h0=1, size=size, mode_tool=mode_tool)
    qu = state_to_qu(state, line_of_sight=0)
    #E,B,Eh,Bh,Qh,Uh = p49_QU2EB.EBfromQU(qu['Q'],qu['U'], return_quharm=True)
    EB = p49_QU2EB.EBfromQU(qu['Q'],qu['U'], return_quharm=True)
    turb_quan.plotter(qu['Q'],qu['U'],EB['E'],EB['B'],
                      'p49c_plots/test_EB_%s_%02f.png'%(wave,theta),norm='ind',
                      axis_labels=qu['axis_labels'] )
    flat = {}
    for f in ['n','H2','Hh','Hv']:
        flat[f] = np.sum(qu[f],axis=0)
        flat[f].shape = (qu[f].shape[0],qu[f].shape[1])
    turb_quan.plotter2([flat['n'],flat['H2'],flat['Hh'], flat['Hv']],
                        'p49c_plots/test_Physics_%s_%0.2f.png'%(wave,theta),
                      labs=['n','B','Hh','Hv'],norm='ind',
                      axis_labels=qu['axis_labels'])
    plots.field_with_streams(flat['n'],flat['Hh'],flat['Hv'],'p49c_plots/with_streams.png')
    outname = 'p49c_plots/test_QhUhEhBh_%s_%0.2f.png'%(wave,theta)
    turb_quan.plotter(np.abs(EB['Qh']), 
                      np.abs(EB['Uh']), 
                      np.abs(EB['Eh']),
                      np.abs(EB['Bh']),
                      outname,labs=['Qh','Uh','Eh','Bh'],norm='positive', axis_labels=qu['axis_labels'])
    print("Plot %s"%outname)
    output = {'fft':fft,'s':state} #,'Q':qu['Q'],'U':qu['U'],'qu':qu}
    output.update(qu)
    output.update(EB)
    return output

def ratios(thing,suffix=""):
    size = thing['s'].cubes['d'].shape[0]
    ks = p49_eigen.make_k_freqs_2d(size, d=1./size)
    thing['psi_ks']=ks
    psi = np.arctan2(ks[1,...],ks[0,...])
    psi=np.ascontiguousarray( psi[:,:size//2+1])
    plt.clf()
    mm = extents()
    npx=2
    npy=3
    fig,axes=plt.subplots(npx,npy)
    plotlist=['Eh','Bh', 'E/B','Qh','Uh']

    for n,f in enumerate(plotlist):
        x=psi.flatten()/np.pi
        if f in ['E/B']:
            Eh = np.ascontiguousarray(np.abs(thing['Eh'])+0)
            Bh = np.ascontiguousarray(np.abs(thing['Bh'])+0)
            ok = Bh>1e-12
            this_field = (Eh[ok]/Bh[ok]).flatten()
            this_field[this_field > 10]=10
            this_field[this_field < .10]=.10
            x=psi[ok].flatten()/np.pi
        else:
            this_field = np.ascontiguousarray(np.abs(thing[f])+0)
        y = np.abs(this_field.flatten())
        x=x[y>0]
        y=y[y>0]
        y[y<1e-12]=1e-12
        marker='.'#['.']*len(y)
        aaa1=n//npy
        aaa2=n%npy
        thisax=axes[aaa1][aaa2]
        if f in ['E/B']:
            e_to_b_ax = thisax
        thisax.set_title(f)
        if len(y) == 0:
            thisax.text(0,0,"No Values")
            continue
        thisax.scatter(x,y,label=f,marker=marker)
        mm(y)
        #plt.yscale('symlog',linthresh=0.1)#1e-18)
        thisax.set_yscale('log')#,linthresh=1e-18)

    def drp(arr):
        a1 = arr.flatten()
        return a1
    ebax=axes[-1][-1]
    ebax.scatter(drp(Eh),drp(Bh))
    ebax.set_xlabel('Eh')
    ebax.set_ylabel('Bh')
    ebax.set_xscale('symlog',linthreshx=1e-13)
    ebax.set_yscale('symlog',linthreshy=1e-13)
    ebax.set_xlim([0,np.round(max([Eh.max(),Bh.max()]),3)])
    ebax.set_ylim([0,np.round(max([Eh.max(),Bh.max()]),3)])

    for n in range(len(plotlist)):
        #axes[n//2][n%2].set_ylim([mm[0],mm[1]])
        axes[n//npy][n%npy].set_ylim(mm[0],mm[1])
    e_to_b_ax.set_ylim(0.1,10)

    outname = 'p49c_plots/p49_psi_quan_%s.png'%suffix
    fig.savefig(outname)
    plt.close(fig)
    print(outname)
    #k3 = ks['k_freq']
    #kint = ks['k_int']
    return psi


def some_fft_things(stuff, line_of_sight):
    field_horizontal = {0:'hy',1:'hz',2:'hx'}[line_of_sight]
    field_vertical   = {0:'hz',1:'hx',2:'hy'}[line_of_sight]
    axis_vertical   = {0:'hy',1:'hx',2:'hx'}[line_of_sight]
    axis_horizontal = {0:'hz',1:'hz',2:'hy'}[line_of_sight]
    axis_labels=[axis_horizontal,axis_vertical]

    data=stuff['s'].cubes
    epsilon = data['d']
    B_sq = data['hx']**2.0 + data['hy']**2.0 + data['hz']**2.0
    size=B_sq.shape[0]
    #epsilon = n
    Q_local = ( epsilon * (((data[field_horizontal])**2.0) - ((data[field_vertical])**2.0))/B_sq )
    #n = data['density']
    #B_sq = data['Bx']**2.0 + data['By']**2.0 + data['Bz']**2.0    
    #epsilon = n
    U_local = (2.0 * epsilon * ((data[field_horizontal]) * (data[field_vertical]))/B_sq)
    U_local[B_sq == 0 ] = 0
    U_flat = np.sum(U_local,axis=line_of_sight)/size
    U_flat.shape = (U_local.shape[0],U_local.shape[1])
    Q_flat = np.sum(Q_local,axis=line_of_sight)/size
    Q_flat.shape = (Q_local.shape[0],U_local.shape[1])
    for field in ['d','hx','hy','hz']:
        print("===== %s ====="%field)
        print(pnz( np.abs(np.fft.rfftn( data[field])),thresh=1e-12))
    for field in ['Q_local','U_local','Q1']:
        print("===== %s ====="%field)
        mmm = {'Q_local':Q_local,
               'U_local':U_local,
               'Q1':(data[field_horizontal])**2.0}
        mmm=mmm[field]
        print(pnz( np.abs(np.fft.rfftn(mmm )),thresh=1e-12))
    for field in ['Qh','Uh','Eh','Bh']:
        print("===== %s ====="%field)
        print(nz( stuff[field], thresh=1e-12))

def from_wave_to_psi(wave,theta,suffix,size=16,mode_tool=modes_k001):
    thing=wave_to_EB(wave,theta,size=size, mode_tool=mode_tool)
    psi=ratios(thing, suffix=suffix)
    #thing['psi']=psi
    return thing

#import k_vecs
#reload(k_vecs)
#k_vecs.kvecs_from_stuff(stuff)
#print(nz(fast['Eh']))
#print(nz(slow['Eh']))
#print(nz(alf['Eh']))
#print(nonzero(fast['Eh']))
#print(nonzero(slow['Eh']))
#print(nonzero(alf['Eh']))
#print(nonzero(fast['Eh']))
#print(nonzero(slow['Eh']))
#print(nonzero(alf['Eh']))
#NEXT: examine fast['Eh']/fast['Bh'] for non-trivial K distribution.
#This will give us the function as a function of Psi.
