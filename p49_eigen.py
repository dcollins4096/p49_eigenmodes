from go_lite_pyt3 import *
import pdb
import numpy as np 
import copy
nar=np.array
import enzo_write
reload(enzo_write)
from p49_print_tools import *
debug = 1
should_i_half = 1.0
def is_iterable(thing):
    if hasattr(thing,'__getitem__') and type(thing) not in [np.float64]:
        return True
    else:
        return False
class fieldthing():
    def __init__(self):
        self.stuff={}
    def get(self,key,default):
        return self.stuff.get(key,default)
    def update(self,other_dict):
        if type(other_dict) is dict:
            self.stuff.update(other_dict)
        if type(other_dict) is type(self):
            self.stuff.update(other_dict.stuff)
    def __getitem__(self,key):
        if key in ['px','py','pz']:
            if key in self.stuff:
                out=self.stuff[key]
            else:
                out = self.stuff[ pv[key]]*self.stuff['d']
                self.stuff[key]=out
        elif key in ['vx','vy','vz','e_specific']:
            if key in self.stuff:
                out = self.stuff[key]
            else:
                out = self.stuff[ vp[key] ]/self.stuff['d']
                self.stuff[key]=out
        else:
            out = self.stuff[key]
        return out
    def __setitem__(self,key,value):
        self.stuff[key]=value
    def __contains__(self,key):
        return key in self.stuff
def wrap_faces(array,field):
    if field not in ['hx','hy','hz']:
        return 
    else:
        if field == 'hx':
            array[-1,:,:]=array[0,:,:]
        elif field == 'hy':
            array[:,-1,:]=array[:,0,:]
        elif field == 'hz':
            array[:,:,-1]=array[:,:,0]
        else:
            print("Wrap Error")
        return 
#def to_ft(fields,mean={}):
#    field_list = ['d','vx','vy','vz','hx','hy','hz','p']
#    for field in 
def make_k_freqs_2d(nk):
    k_freq = np.zeros([2,nk,nk])
    k1=np.fft.fftfreq(nk,d=1)
    #kx, ky, kz = np.meshgrid(k1,k1,k1)
    #k_freq[0,...]=kx
    #k_freq[1,...]=ky
    #k_freq[2,...]=kz
    x = np.repeat(k1,nk)
    x.shape = (nk,nk)
    y = np.repeat(k1,nk)
    y.shape = (nk,nk)
    y=y.swapaxes(0,1)
    k_freq[0,...]=x
    k_freq[1,...]=y
    return k_freq
def make_k_freqs(nk,real=True, d=1):

    ny = nk
    nz = nk
    if real:
        nz = nk//2+1

    k_freq = np.zeros([3,nk,ny,nz])
    k1=np.fft.fftfreq(nk,d=d)
    #kx, ky, kz = np.meshgrid(k1,k1,k1)
    #k_freq[0,...]=kx
    #k_freq[1,...]=ky
    #k_freq[2,...]=kz
    x = np.repeat(k1,nk*nk)
    x.shape = (nk,nk,nk)
    y = np.repeat(k1,nk*nk)
    y.shape = (nk,nk,nk)
    y=y.swapaxes(0,1)
    z = np.repeat(k1,nk*nk)
    z.shape = (nk,nk,nk)
    z=z.swapaxes(0,2)
    if real:
        x = x[:,:,:nz]
        y = y[:,:,:nz]
        z = z[:,:,:nz]
    k_freq[0,...]=x
    k_freq[1,...]=y
    k_freq[2,...]=z
    return k_freq
def make_k_freqs_and_int(nk,real=True):
    k_freq = make_k_freqs(nk,real=real)
    k_int = make_k_freqs(nk,real=real,d=1./nk)
    return {'k_freq':k_freq,'k_int':k_int}
def make_waveset(data):
    #Start with XYZ: data['cubes']
    #hat(XYZ): FFT of real data.  data['ffts'], with get_ffts
    #hat(ABC): Rotated FFT data. data['waveset'].wave_frame
    #hat(FSA): Eigen vectors in FFT space.  data['waveset'].wave_content
    #
    #Also created: waveset.hat_system, contains eigen vectors
    #              waveset.a_unit, b_unit, c_unit: unit vectors.
    # Call order:
    #        rotate_back; makes K vectors, calls;;
    #           this_system = waves(mean)
    #           this_system.project_to_waves(k_all, ffts, means={})
    #              this_system.fields_to_wave_frame(k_all_in,fields=ffts)
    #                +this_system.wave_frame (another system)
    #                +this_system.a_unit, b_unit, c_unit, B0hat
    #                +this_system.hat_system (rotated)
    #                +this_system.wave_content = self.to_waves; to wave frame.
    #                       (self.wave_frame * self.hat_system.left)
    #waveset.wave_frame = hat(ABC)
    #waveset.hat_system = hat(ABC) eigen system
    #waveset.wave_content = hat(FSA) wave system.
    #waveset.a_unit, b_unit, c_unit
    data['ffts']  = get_ffts(data['cubes'], data['means'])
    kall,waveset=rotate_back(data['ffts'], data['means'])
    data['kall']=kall
    data['waveset']=waveset
    return data
def get_cubes_cg(ds,mean={}):
    print(ds['DomainLeftEdge'].astype('float'))
    cg = ds.covering_grid(0,ds['DomainLeftEdge'].astype('float'),ds['TopGridDimensions'])
    field_list = ['d','vx','vy','vz','hx','hy','hz','p']
    map_to_label ={'d':'density','vx':'x-velocity','vy':'y-velocity','vz':'z-velocity',
                   'hx':'Bx','hy':'By','hz':'Bz','e':'TotalEnergy','p':'pressure'}
    cubes=fieldthing()
    means=fieldthing()
    for field in field_list:
        cubes[field] = cg[map_to_label[field] ].v
        means[field] = np.mean( cubes[field])
    for field in ['px','py','pz']:
        means[field] = cubes[field].mean()
    #then get momentum.
    means['Gamma']=ds['Gamma']
    data={'cubes':cubes, 'means':means}
    data = make_waveset(data)
    return data

def get_ffts(cubes, means={}, real=True):
    field_list = ['d','px','py','pz','hx','hy','hz','p']
    if real:
        the_fft = np.fft.rfftn
    else:
        the_fft = np.fft.fftn
    if means == 'get':
        means={}
        for field in field_list:
            means[field] = np.mean(cubes[field])

    ffts = fieldthing()

    kludge={}
    for field in field_list:
        this_cube = cubes[field]-means.get(field,0)
        #if field in ['vx','vy','vz']:
        #    this_cube *= cubes['d']
        #kludge[field]=this_cube

        ffts[field] = the_fft(this_cube)/(np.size(this_cube))/should_i_half #half?

    return ffts

def rotate_back(ffts,means, real=True):

    the_means = means
    if hasattr(the_means,'stuff'):
        the_means=means.stuff
    this_system = waves(form='rb96',**the_means)
    size = nar(ffts['d'].shape)
    #k_all = np.mgrid[:size[0],:size[1],:size[2]]
    k_all = make_k_freqs(size[0], real=real)
    
    #rotate from xyz to abc
    this_system.project_to_waves(k_all, ffts, means={})

    return k_all, this_system



#def to_ampl(ds, mean={}):
#
#    ds.

pv={'px':'vx','py':'vy','pz':'vz'}
vp={'vx':'px','vy':'py','vz':'pz','e_specific':'e'}


class waves():
    def pull_rot(self,n):
        out = {}
        for field in self.rot:
            out[field]=self.rot[field][n]
        return out
    def compute_unit_vectors(self,k_all_in):
        rank = 3
        k_all = k_all_in.astype('float') #just to be sure.
        B0 = self.b0
        B0_mag = (B0**2).sum()**0.5
        h_unit = B0/B0_mag
        self.h_unit=h_unit

        #compute unit vectors
        a_unit = np.zeros_like(k_all).astype('float')
        b_unit = np.zeros_like(a_unit)
        
        knorm = (k_all[0,...]**2+k_all[1,...]**2+k_all[2,...]**2)**0.5
        ok=knorm>0
        for dim in range(rank):
            a_unit[dim,...][ok] = k_all[dim,...][ok]/knorm[ok]
        self.a_unit=a_unit

        c_unit = np.zeros_like(a_unit)
        c_unit[0,...] = a_unit[1,...]*h_unit[2,...]-a_unit[2,...]*h_unit[1,...]
        c_unit[1,...] = a_unit[2,...]*h_unit[0,...]-a_unit[0,...]*h_unit[2,...]
        c_unit[2,...] = a_unit[0,...]*h_unit[1,...]-a_unit[1,...]*h_unit[0,...]
        c_unit_norm  = np.sqrt(c_unit[0,...]**2+c_unit[1,...]**2+c_unit[2,...]**2)


        ok = c_unit_norm > 0
        for dim in range(rank):
            c_unit[dim,...][ok] = c_unit[dim,...][ok]/c_unit_norm[ok]
        self.c_unit_norm=c_unit_norm
        self.c_unit=c_unit

        #if (ok==False).any():
        #    print("Warning: Degenerate system, k||B0. ||c(kpq)||=0 for some kpq Please check corner cases.")

        b_unit = np.zeros_like(a_unit)
        b_unit[0,...] = c_unit[1,...]*a_unit[2,...]-c_unit[2,...]*a_unit[1,...]
        b_unit[1,...] = c_unit[2,...]*a_unit[0,...]-c_unit[0,...]*a_unit[2,...]
        b_unit[2,...] = c_unit[0,...]*a_unit[1,...]-c_unit[1,...]*a_unit[0,...]
        b_unit_norm  = np.sqrt(b_unit[0,...]**2+b_unit[1,...]**2+b_unit[2,...]**2)

        #ok=b_unit_norm>0
        #for dim in range(rank):
        #    b_unit[dim,...][ok] = b_unit[dim,...][ok]/ b_unit_norm[ok]
        #if (ok==False).any():
        #    print("Warning: Degenerate system, k||B0. ||b(kpq)||=0 for some kpq Please check corner cases.")
        self.b_unit=b_unit



        #CORRECT FOR k||B0.  Not done yet.
        #along = np.logical_and( np.abs(a_unit[0,...]-b_unit[0]) < dx[0],\
        #                        np.abs(a_unit[1,...]-b_unit[1]) < dx[1])
        #along = np.logical_and( np.abs(a_unit[2,...]-b_unit[2]) < dx[2],\
        #                        along)

        B0hat = np.zeros_like(a_unit)
        B0hat[0,...] = self.b0[0]*a_unit[0,...]+self.b0[1]*a_unit[1,...]+self.b0[2]*a_unit[2,...]
        B0hat[1,...] = self.b0[0]*b_unit[0,...]+self.b0[1]*b_unit[1,...]+self.b0[2]*b_unit[2,...]
        B0hat[2,...] = self.b0[0]*c_unit[0,...]+self.b0[1]*c_unit[1,...]+self.b0[2]*c_unit[2,...]
        self.B0hat=B0hat

    def fields_to_wave_frame(self,k_all_in, fields):
        """compute eigen vectors for each k_all_in.
        Creates:
                self.wave_frame
                self.a_unit, b_unit, c_unit, B0hat
                hat_system (rotated)
                """

        self.compute_unit_vectors(k_all_in)

        scalars={}
        for field in ['d','p','e', 'px','py','pz']:
            scalars[field]=np.zeros_like(self.B0hat[0,...])+self.quan[field]

        hat_system = waves(hx=self.B0hat[0,...], hy=self.B0hat[1,...], hz=self.B0hat[2,...],
                           Gamma=self.Gamma,form=self.form,**scalars)
        self.hat_system = hat_system
        self.wave_frame=fieldthing()
        def newer():
            return np.zeros_like(self.a_unit[0,...])
        for f in ['px','py','pz','hx','hy','hz']:
            self.wave_frame[f]=np.zeros_like(self.a_unit[0,...])

        self.wave_frame['d'] = fields['d']
        if 'e' in fields:
            self.wave_frame['e'] =fields['e']
        if 'p' in fields:
            self.wave_frame['p'] =fields['p']
        self.wave_frame.update(self.rotate_to_abc(fields))

    #def check_orthonormality(self,k_all_in,fields, means=None):
    def check_orthonormality(self, right=None,left=None):
        field_list = ['d','px','py','pz','hx','hy','hz','p']
        if left is None:
            left = self.left
        if right is None:
            right = self.right
        #field_list = ['hx','hy','hz','p']
        dbg=0

        form = " %8s"
        formn = " %8.2e"
        tensor = {}
        wave_list_1=['f-', 'a-','s-','c','f+','a+','s+']
        wave_list_2=['f-', 'a-','s-','c','f+','a+','s+']
        for wave_1 in wave_list_1:
            tensor[wave_1]={}
            for wave_2 in wave_list_2:
                tensor[wave_1][wave_2]=0.0
                for field in field_list:
                    this_part = right[wave_1][field]*left[wave_2][field]
                    tensor[wave_1][wave_2] += this_part
                    if dbg > 1:
                        print("%s %s %s tot= %s L %s R %s val %s"%(form%wave_1,form%wave_2,form%field,\
                                 formn%tensor[wave_1][wave_2], \
                                 formn%left[wave_2][field],\
                                 formn%right[wave_1][field],\
                                 formn%this_part))

        if dbg > 1:
            print("alph_f**2 a**2 %0.2e"%(self.alph_f**2*self.aa**2))
        print( form%" "+form*len(wave_list_2)%tuple(wave_list_2))
        for wave_1 in wave_list_1:
            this_string = form%wave_1
            for wave_2 in wave_list_2:
                val = tensor[wave_1][wave_2]
                if val != 0.0 and np.abs(val) > 2e-16:
                    this_string += formn%tensor[wave_1][wave_2]
                else:
                    this_string += form%"."
            print(this_string)

        if dbg > 1:
            squares2= self.cf**2*self.alph_f**2+self.cs**2*self.alph_s**2
            print("alph_f**2 + alph_s**2 = 1; %0.2e"%(self.alph_f**2+self.alph_s**2))
            print("alph_f**2 cf**2 + alph_s**2 cs**2 = a**2; %0.2e %0.2e %0.2e"%\
                    (squares2, self.aa**2,1-squares2/self.aa**2 ))
            print("betay**2+betaz**2 = 1; %0.2e"%(self.betay**2+self.betaz**2))
        #print(tensor)

    def project_to_waves(self,k_all_in,fields, means=None):
        field_list = ['d','px','py','pz','hx','hy','hz','p']
        self.fields_to_wave_frame(k_all_in,fields)
        if means is None:
            means = self.quan

        self.wave_content=self.to_waves(self.wave_frame,self.hat_system.left,means)

    def to_waves(self,fields,lefts,means={}):
        wave_content = {}
        field_list = ['d','px','py','pz','hx','hy','hz','p']
        #print("Means in project", means)
        for wave in ['f-', 'a-','s-','c','f+','a+','s+']:
            wave_content[wave] = np.zeros_like(fields['d'])
            for field in field_list:
                #self.wave_content[wave] += (fields[field]-means[field])*self.left[wave][field]
                wave_content[wave] += (fields[field]-means.get(field,0))*lefts[wave][field]
        return wave_content

                





    def rotate_to_abc(self,fields):
        """from xyz to abc, where abc is defined by the wave vector and the field."""
        back = fieldthing()
        for f in ['px','py','pz','hx','hy','hz']:
            back[f]=np.zeros_like(self.a_unit)
        back['px']  = self.a_unit[0,...]*fields['px']
        back['px'] += self.a_unit[1,...]*fields['py']
        back['px'] += self.a_unit[2,...]*fields['pz']

        back['py']  = self.b_unit[0,...]*fields['px']
        back['py'] += self.b_unit[1,...]*fields['py']
        back['py'] += self.b_unit[2,...]*fields['pz']

        back['pz']  = self.c_unit[0,...]*fields['px']
        back['pz'] += self.c_unit[1,...]*fields['py']
        back['pz'] += self.c_unit[2,...]*fields['pz']

        back['hx']  = self.a_unit[0,...]*fields['hx']
        back['hx'] += self.a_unit[1,...]*fields['hy']
        back['hx'] += self.a_unit[2,...]*fields['hz']

        back['hy']  = self.b_unit[0,...]*fields['hx']
        back['hy'] += self.b_unit[1,...]*fields['hy']
        back['hy'] += self.b_unit[2,...]*fields['hz']

        back['hz']  = self.c_unit[0,...]*fields['hx']
        back['hz'] += self.c_unit[1,...]*fields['hy']
        back['hz'] += self.c_unit[2,...]*fields['hz']
        return back

    def rotate_to_xyz(self,fields):
        """from hat system, ahat, bhat, chat, along k and H0
        to xyz"""
        rot = fieldthing() 
        for f in ['px','py','pz','hx','hy','hz']:
            rot[f]=np.zeros_like(self.a_unit)
        rot['px']  = self.a_unit[0,...]*fields['px']
        rot['px'] += self.b_unit[0,...]*fields['py']
        rot['px'] += self.c_unit[0,...]*fields['pz']

        rot['py']  = self.a_unit[1,...]*fields['px']
        rot['py'] += self.b_unit[1,...]*fields['py']
        rot['py'] += self.c_unit[1,...]*fields['pz']

        rot['pz']  = self.a_unit[2,...]*fields['px']
        rot['pz'] += self.b_unit[2,...]*fields['py']
        rot['pz'] += self.c_unit[2,...]*fields['pz']

        rot['hx']  = self.a_unit[0,...]*fields['hx']
        rot['hx'] += self.b_unit[0,...]*fields['hy']
        rot['hx'] += self.c_unit[0,...]*fields['hz']

        rot['hy']  = self.a_unit[1,...]*fields['hx']
        rot['hy'] += self.b_unit[1,...]*fields['hy']
        rot['hy'] += self.c_unit[1,...]*fields['hz']

        rot['hz']  = self.a_unit[2,...]*fields['hx']
        rot['hz'] += self.b_unit[2,...]*fields['hy']
        rot['hz'] += self.c_unit[2,...]*fields['hz']
        return rot
    def wave_to_fields(self,k_all_in, wave=None):
        if wave is None:
            wave = self.wave
        self.compute_unit_vectors(k_all_in)

        scalars={}
        for field in ['d','p','e', 'px','py','pz']:
            scalars[field]=np.zeros_like(self.B0hat[0,...])+self.quan[field]

        hat_system = waves(hx=self.B0hat[0,...], hy=self.B0hat[1,...], hz=self.B0hat[2,...],
                           Gamma=self.Gamma,form=self.form,**scalars)
        
        self.rot=fieldthing()
        for f in ['px','py','pz','hx','hy','hz']:
            self.rot[f]=np.zeros_like(self.a_unit)

        self.rot['d'] = hat_system.right[wave]['d']
        if 'e' in hat_system.right[wave]:
            self.rot['e'] = hat_system.right[wave]['e']
        if 'p' in hat_system.right[wave]:
            self.rot['p'] = hat_system.right[wave]['p']
        self.rot.update(self.rotate_to_xyz(hat_system.right[wave]))
        self.hat_system=hat_system

    def __init__(self, d=1.0,vx=0.0,vy=0.0,vz=0.0,hx=1.0,hy=1.41421,hz=0.5,Gamma=1.6666666667,e=None,p=None, write=True,
                 form='rb96', HydroMethod=4, **kwargs):
        #p = 0.6
        #Gamma=1.6666666667
        self.HydroMethod = HydroMethod
        self.Gamma=Gamma
        self.form=form
        self.wave=None
        self.cubes=None
        if p is not None:
            self.egas=p/((Gamma-1)*d)
            egas=self.egas
            e = 0.5*(vx*vx+vy*vy+vz*vz)+0.5*(hx*hx+hy*hy+hz*hz)/d + egas
        self.quan={'d':d,'vx':vx,'vy':vy,'vz':vz,'hx':hx,'hy':hy,'hz':hz,'p':p,'Gamma':Gamma,'egas':egas,'e':e,\
                   'px':vx*d,'py':vy*d,'pz':vz*d}
        for a,b in [['d','Density'],['vx','x-velocity'],['vz','z-velocity'],['vy','y-velocity'],['e','TotalEnergy'], ['p','GasPressure']]:
            self.quan[b]=self.quan[a]
        self.v0 = np.array([vx,vy,vz])
        self.b0 = np.array([hx,hy,hz])
        self.p = p
        self.d = d
        self.left, self.right = self.eigen(d, vx, vy, vz, hx,hy,hz,e, Gamma=Gamma, form = self.form)
        self.speeds = self.speeds(d, vx, vy, vz, hx,hy,hz,p=p, Gamma=Gamma)
    def __getitem__(self,field):
        this_pert = self.right[self.wave]
        #if hasattr(this_pert['hy'],'__getitem__'):
        if field is 'dv':
            return nar([ this_pert['vx'], this_pert['vy'], this_pert['vz']])
        if field is 'db':
            return nar([ 0, this_pert['hy'], this_pert['hz']])
        if field is 'de':
            return nar(this_pert['d'])
        if field is 'drho':
            return nar(this_pert['e'])
        if field is 'dall':
            return this_pert
        if field is 'b0':
            return self.b0
        if field is 'v0':
            return self.v0
        if field is 'all':
            return self.quan

    def perturb(self,k_rot=None,base_size=None,ampl=1e-6,directory=".", write=True, wave=None, pert_shape='fft',
                  start=True, real_fft=True,single_wave_selector=True):


        field_list = ['d','px','py','pz','hx','hy','hz','e']
        enzo_field_list = ['d','vx','vy','vz','hx','hy','hz','e_specific']
        if self.form =='rb96':
            field_list = ['d','px','py','pz','hx','hy','hz','p']
            enzo_field_list = ['d','vx','vy','vz','hx','hy','hz','p']
        map_to_label ={'d':'density','vx':'x-velocity','vy':'y-velocity','vz':'z-velocity',
                       'hx':'Bx','hy':'By','hz':'Bz','e_specific':'TotalEnergy','p':'GasPressure'}
        #print('need to sort when you need all hats')
        self.all_hats = {}
        if self.HydroMethod == 6:
            face_offset = {'hx':nar([1,0,0]),'hy':nar([0,1,0]),'hz':nar([0,0,1])}
        else:
            face_offset={}
        if start:
            #self.cubes_test={}
            self.all_p = {}
            self.temp_right_back={}
            for f in field_list:
                self.all_p[f]=0.
        cube_slice=[slice(0,base_size[0]),slice(0,base_size[1]),slice(0,base_size[2])]

        self.directory=directory
        if not hasattr(self,'hat_box'):
            self.hat_box={}
            print("make hat box")
        if pert_shape == 'fft':
            self.wave_to_fields(k_rot, wave) 
            kint = k_rot.astype(int)
            for f in field_list:
                if real_fft:
                    this_size=np.array([base_size[0],base_size[1],base_size[2]//2+1])
                else:
                    this_size=base_size #np.array([base_size[0],base_size[1],base_size[2]//2+1])
                self.all_hats[f]=np.zeros(this_size)*1j
                if single_wave_selector:
                    self.all_hats[f][kint[0,...],kint[1,...],kint[2,...]] = self.rot[f]*ampl
                else:
                    self.all_hats[f] = self.rot[f]*ampl
                    #print("shapes: rot    %s ampl %s"%(str(self.rot[f].shape),str(ampl.shape)))

                #print(    "N nonzero K %s "%f + str(nz(self.all_hats[f])))
                #pdb.set_trace()
                #print( "put %s %0.2e"%(f,(self.rot[f]*ampl)[0]))
                if not real_fft:
                    #self.all_hats[f][kint[0,...],kint[1,...],-kint[2,...]] = (self.all_hats[f][kint[0,...],kint[1,...],kint[2,...]] ).conj()
                    if kint is not None:
                        self.all_hats[f][-kint[0,...],-kint[1,...],-kint[2,...]] = (self.all_hats[f][kint[0,...],kint[1,...],kint[2,...]] ).conj()
                    else:
                        print("Shoot, i haven't sorted the conjugation yet.  I hope you did.")
#            for f in  field_list: #there is a bug in an earlier numpy.
#                print("=== pre %3s "%f+str(len(nonzero(self.all_hats[f]))))
#                print("=== pre %3s "%f+str(nz(self.all_hats[f])))
#            self.hat_box[wave] = self.all_hats
#            print("wave",wave)
#            for mw in self.hat_box:
#                print(' -- %s -- '%mw)
#                print_nz_2(self.hat_box[mw]['d'])
            for f in  field_list:
                if 'tmp' in dir():
                    del tmp
                if real_fft:
                    tmp=np.fft.irfftn(self.all_hats[f])
                    #np.fft.irfftn(self.all_hats[f])
                else:
                    tmp=np.fft.ifftn(self.all_hats[f])
                #for the bug.
                #print("=== pos %3s "%f+str(len(nonzero(self.all_hats[f]))))
                real_mean = np.mean(np.abs(tmp.real))
                imag_mean = np.mean(np.abs(tmp.imag))
                if (real_mean+imag_mean)>1e-16:
                    if imag_mean/(real_mean+imag_mean) > 1e-9:
                        print("Warning: large imaginary component")
                self.all_p[f]+=tmp.real*tmp.size*should_i_half
        elif pert_shape == 'square_x':
            size = base_size #+face_offset.get(f,0)
            amplitude = np.zeros(size)
            amplitude[:size[0]//2,:,:] = 1
            amplitude[size[0]//2:,:,:] = -1
            for f in field_list:
                self.all_p[f] = amplitude*ampl[0]*self.right[wave][f] 
            #print("WTF CLOWN ", ampl[0]*self.right[wave][f],wave,f)

            #


        #if self.cubes is None:
        #cubes are always generated fresh.
        #Perturbations are addative.

        self.cubes=fieldthing() 
        for field in field_list:
            size = base_size+face_offset.get(field,0)
            if field in self.cubes:
                this_set = self.cubes[field]
                #print("recalw: %s %s %s"%(field,map_to_label[field],str(this_set.shape)))
            else:
                this_set = np.ones(size)*self.quan[field]
                #print("Setup: %s %s %s"%(field,map_to_label[field],str(this_set.shape)))
            #if field in ['vx','py','pz','e']:
            #    if np.abs(np.mean(self.cubes['density']) - 1.0) > 1e-7:
            #        print("YOU MAY NOT WANT TO BE MULTIPLYING ENERGY BY DENSITY INITALLLY.")
            #        #I think actually we don't, but it's what Enzo currently does.
            #    this_set *= self.cubes['density']
            self.cubes[field]=this_set
        for field in field_list:
            size = base_size+face_offset.get(field,0)
            self.cubes[field][cube_slice] += self.all_p[field] 
            #self.temp_right_back[field]=np.fft.fftn(self.cubes[map_to_label[field]][cube_slice] )
            #self.temp_right_back[field]=np.fft.rfftn(self.cubes[map_to_label[field]][cube_slice] )
            #self.cubes_test[field]  =  self.cubes[map_to_label[field]][cube_slice] - self.quan[field]
        for field in enzo_field_list: 
            #self.cubes[map_to_label[field]] = wrap_faces(self.cubes[map_to_label[field]], field)
            if self.HydroMethod == 6:
                wrap_faces(self.cubes[field], field)
            this_filename = "%s/%s_16.h5"%(directory,map_to_label[field])
            #if field in ['px','py','pz','e']:
            #    vel_field = vp[field]
            #    self.cubes[map_to_label[vel_field]] = self.cubes[map_to_label['d']]self.cubes[map_to_label['d']]
            if write:
                enzo_write.dump_h5(self.cubes[field],this_filename)
                if debug >  0:
                    print("wrote "+this_filename + " with shape "+str(self.cubes[field].shape))



    def right_roe_balsara(self,d,vx,vy,vz,hx,hy,hz,eng, EquationOfState = 0, Gamma=5./3, IsothermalSoundSpeed=1):

        #float d, float vx, float vy, float vz, 
        #float hx, float by, float bz, float eng, 
        #float right[][7])

        #From Roe&Balsara 1996, SIAM
        #hx = magnetic field
        #bx = hx/sqrt(rho)


        p = -1
        bp = 0.5*(hx*hx + hy*hy + hz*hz );
        sqrt2 = np.sqrt(2.0);
        sqrt2i= 1.0/sqrt2;
        sqrtD = np.sqrt(1.*d);
        sqrtDi = 1.0/sqrtD;
        sbx  = np.sign(hx);
        og1  = 1.0/(Gamma - 1);

        bx = hx*sqrtDi
        by = hy*sqrtDi
        bz = hz*sqrtDi

        if EquationOfState == 0 :
            p = (Gamma -1.0 ) * (eng - 0.5* d * (vx*vx + vy*vy + vz * vz ) - 0.5*(bx*bx + by*by + bz*bz));
        else:
            p = IsothermalSoundSpeed*IsothermalSoundSpeed * d;

        #compute wave speeds
        if EquationOfState == 0:
            aa = np.sqrt( 1.0*Gamma* p/d )
        else:
            aa = IsothermalSoundSpeed;

        cs =np.sqrt( 0.5*( aa*aa + 2*bp/d -np.sqrt( pow( (aa*aa + 2*bp/d ),2) - 4* aa*aa*bx*bx/d ) ) );
        ca =np.sqrt( bx*bx/d ); 
        cf =np.sqrt( 0.5*( aa*aa + 2*bp/d +np.sqrt( pow( (aa*aa + 2*bp/d ),2) - 4* aa*aa*bx*bx/d ) ) );
        #print("yay 5", Gamma)
        #self.thing1 =  p #aa# (aa*aa )#+ 2*bp/d ) #- 4* aa*aa*bx*bx/d 
        
        #compute ancilary values
        #The normalization of alph_f may change. This normalization uses Ryu
        #& Jones, but Balsara may be more robust.

        betay=np.zeros_like(by)+sqrt2i
        betaz=np.zeros_like(by)+sqrt2i
        alph_f=np.zeros_like(by)+1
        alph_s=np.zeros_like(by)+1
        #non-degenerate points
        non_deg = np.logical_or( by != 0.0 , bz != 0.0 )
        bt =np.zeros_like(by)
        if is_iterable(by):
            bt[non_deg] = np.sqrt( by[non_deg]*by[non_deg] + bz[non_deg]*bz[non_deg] );
            betay[non_deg] = by[non_deg]/bt[non_deg];
            betaz[non_deg] = bz[non_deg]/bt[non_deg];
            alph_f[non_deg] = np.sqrt(np.abs((cs[non_deg]*cs[non_deg]-aa[non_deg]*aa[non_deg])/\
                                             (cf[non_deg]*cf[non_deg]-cs[non_deg]*cs[non_deg])) );
            alph_s[non_deg] = np.sqrt(np.abs((cf[non_deg]*cf[non_deg]-aa[non_deg]*aa[non_deg])/
                                             (cf[non_deg]*cf[non_deg]-cs[non_deg]*cs[non_deg])));
        else:
            if non_deg:
                bt = np.sqrt( by*by + bz*bz );
                betay = by/bt;
                betaz = bz/bt;
                alph_f = np.sqrt( (aa*aa-cs*cs)/(cf*cf-cs*cs) );
                alph_s = np.sqrt( (cf*cf-aa*aa)/(cf*cf-cs*cs) );
      
        self.alph_f = alph_f
        self.alph_s = alph_s
        self.betay = betay
        self.betaz = betaz
        self.cf = cf
        self.ca = ca
        self.cs = cs
        self.aa = aa
        #the vectors
        right={}
        left={}
        for w in ['f-','a-','s-','c','f+','a+','s+']:
            right[w] = {}
            left[w] = {}
            #for f in ['d','vx','vy','vz','bx','by','bz','e', 'p']
            for f in ['d','px','py','pz','hx','hy','hz', 'p']:
                right[w][f]=np.zeros_like(hy)
                left[w][f]=np.zeros_like(hy)
      
        #fast, left
        over_two_a2 =  1./(2*aa**2)
        self.over_two_a2 = over_two_a2
        right['f-']['d'] = alph_f*d
        right['f-']['px'] = -1*alph_f*cf         # + alph_f*vx; #these second terms
        right['f-']['py'] = +1*alph_s*betay*cs*sbx#+ alph_f*vy  #may be
        right['f-']['pz'] = +1*alph_s*betaz*cs*sbx#+ alph_f*vz  #incorrect. 
        right['f-']['hy'] = alph_s*sqrtD*aa*betay
        right['f-']['hz'] = alph_s*sqrtD*aa*betaz
        right['f-']['p']  = alph_f*d*aa*aa

        #left['f-']['d'] = 0. #already zero
        left['f-']['px'] = -alph_f*cf*over_two_a2
        left['f-']['py'] = +alph_s*cs*betay*sbx*over_two_a2
        left['f-']['pz'] = +alph_s*cs*betaz*sbx*over_two_a2
        left['f-']['hy']  = alph_s*aa*betay*sqrtDi*over_two_a2
        left['f-']['hz']  = alph_s*aa*betaz*sqrtDi*over_two_a2
        left['f-']['p']   = alph_f/d*over_two_a2


        #fast, right
        right['f+']['d'] = alph_f*d
        right['f+']['px'] = +1*alph_f*cf           + alph_f*vx; #these second terms
        right['f+']['py'] = -1*alph_s*betay*cs*sbx + alph_f*vy  #may be
        right['f+']['pz'] = -1*alph_s*betaz*cs*sbx + alph_f*vz  #incorrect. 
        right['f+']['hy'] = alph_s*sqrtDi*aa*betay
        right['f+']['hz'] = alph_s*sqrtDi*aa*betaz
        right['f+']['p'] =  alph_f*d*aa*aa

        over_two_a2 = 1./(2*aa**2)
        #left['f+']['d'] = 0. #already zero.
        left['f+']['px'] = +alph_f*cf*over_two_a2
        left['f+']['py'] = -alph_s*cs*betay*sbx*over_two_a2
        left['f+']['pz'] = -alph_s*cs*betaz*sbx*over_two_a2
        left['f+']['hy'] =  alph_s*aa*betay*sqrtDi*over_two_a2
        left['f+']['hz'] =  alph_s*aa*betaz*sqrtDi*over_two_a2
        left['f+']['p']  =  alph_f/d*over_two_a2

        #slow, right
        right['s+']['d']  = alph_s*d
        right['s+']['px'] = +1*alph_s*cs          # + alph_f*vx; #these second terms
        right['s+']['py'] = +1*alph_f*betay*cf*sbx# + alph_f*vy  #may be
        right['s+']['pz'] = +1*alph_f*betaz*cf*sbx# + alph_f*vz  #incorrect. 
        right['s+']['hy'] = -alph_f*sqrtD*aa*betay
        right['s+']['hz'] = -alph_f*sqrtD*aa*betaz
        right['s+']['p']  =  alph_s*d*aa*aa

        over_two_a2 = 1./(2*aa**2)
        #left['s+']['d'] = 0 #already zero.
        left['s+']['px'] = +alph_s*cs*over_two_a2
        left['s+']['py'] = +alph_f*cf*betay*sbx*over_two_a2
        left['s+']['pz'] = +alph_f*cf*betaz*sbx*over_two_a2
        left['s+']['hy'] = -alph_f*aa*betay*sqrtDi*over_two_a2
        left['s+']['hz'] = -alph_f*aa*betaz*sqrtDi*over_two_a2
        left['s+']['p']  =  alph_s/d*over_two_a2

        #aaa
        #slow, left
        right['s-']['d']  = alph_s*d
        right['s-']['px'] = -1.*alph_s*cs          # + alph_f*vx; #these second terms
        right['s-']['py'] = -1.*alph_f*betay*cf*sbx# + alph_f*vy  #may be
        right['s-']['pz'] = -1.*alph_f*betaz*cf*sbx# + alph_f*vz  #incorrect. 
        right['s-']['hy'] = -alph_f*sqrtD*aa*betay
        right['s-']['hz'] = -alph_f*sqrtD*aa*betaz
        right['s-']['p']  =  alph_s*d*aa*aa

        over_two_a2 = 1./(2*aa**2)
        #left['s-']['d'] = 0. #already zero
        left['s-']['px'] = -1.*alph_s*cs*over_two_a2
        left['s-']['py'] = -1.*alph_f*cf*betay*sbx*over_two_a2
        left['s-']['pz'] = -1.*alph_f*cf*betaz*sbx*over_two_a2
        left['s-']['hy'] = -alph_f*aa*betay*sqrtDi*over_two_a2
        left['s-']['hz'] = -alph_f*aa*betaz*sqrtDi*over_two_a2
        left['s-']['p']  =  alph_s/d*over_two_a2
        #bbb




        #alfven][left
        #right['a-']['d'] = 0; already zero
        #right['a-']['px'] = 0; already zero.
        right['a-']['py'] = -1.*betaz*sbx;
        right['a-']['pz'] = +1.*betay*sbx;
        right['a-']['hy'] = -betaz*sqrtD;
        right['a-']['hz'] =  betay*sqrtD;
        #right['a-']['p'] = 0; already zero

        #left['a-']['d'] = 0; already zero
        #left['a-']['px'] = 0; already zero.
        left['a-']['py'] = -1.*betaz*sbx*0.5;
        left['a-']['pz'] = +1.*betay*sbx*0.5;
        left['a-']['hy'] = -1.*betaz*sqrtDi*sbx*0.5;
        left['a-']['hz'] = +1.*betay*sqrtDi*sbx*0.5;
        #left['a-']['p'] = 0; already zero

        #alfven][right
        #right['a+']['d'] = 0; already zero
        #right['a+']['px'] = 0; already zero.
        right['a+']['py'] = +1.*betaz*sbx;
        right['a+']['pz'] = -1.*betay*sbx;
        right['a+']['hy'] = -betaz*sqrtD;
        right['a+']['hz'] =  betay*sqrtD;
        #right['a+']['p'] = 0; already zero

        #left['a+']['d'] = 0; already zero
        #left['a+']['px'] = 0; already zero.
        left['a+']['py'] = +1.*betaz*sbx*0.5;
        left['a+']['pz'] = -1.*betay*sbx*0.5;
        left['a+']['hy'] = -1.*betaz*sqrtDi*sbx*0.5;
        left['a+']['hz'] = +1.*betay*sqrtDi*sbx*0.5;
        #left['a+']['p'] = 0; already zero

      
        right['c']['d'] = 1;
        #right['c']['px'] = 0. #vx; Already zero.
        #right['c']['py'] = 0. #vy;
        #right['c']['pz'] = 0. #vz;
        #right['c']['hy'] = 0;
        #right['c']['hz'] = 0;
        #right['c']['p'] =  0;

        left['c']['d'] = 1;
        #left['c']['px'] = 0.; Already zero
        #left['c']['py'] = 0.;
        #left['c']['pz'] = 0.;
        #left['c']['hy'] = 0;
        #left['c']['hz'] = 0;
        left['c']['p'] =  -1./aa**2; #0.5*(vx*vx+vy*vy+vz*vz);

        return left,right
      
    def eigen(self,d,vx,vy,vz,bx,by,bz,eng, EquationOfState = 0, Gamma=5./3, IsothermalSoundSpeed=1, obj=None, form='rj95'):
        #B = H in this function
        if form == 'rj95':
            #as formulated by Ryu and Jones 1995
            left, right= self.eigen_ryu_jones(d,vx,vy,vz,bx,by,bz,eng, \
                    EquationOfState, Gamma, IsothermalSoundSpeed)
        elif form == 'rb96':
            #as formulated by Roe andBalsara, 1996
            left, right = self.right_roe_balsara(d,vx,vy,vz,bx,by,bz,eng, \
                    EquationOfState, Gamma, IsothermalSoundSpeed)
        return left, right


    def eigen_ryu_jones(self,d,vx,vy,vz,bx,by,bz,eng, EquationOfState = 0, Gamma=5./3, IsothermalSoundSpeed=1):
        """B=h in this function"""

        #float d, float vx, float vy, float vz, 
        #float bx, float by, float bz, float eng, 
        #float right[][7])

        #Eigen vectors taken from Ryu & Jones, ApJ 442:228-258, 1995
        #Normalization starting on p. 231.

        p = -1
        bp = 0.5*(bx*bx + by*by + bz*bz );
        sqrt2 = np.sqrt(2.0);
        sqrt2i= 1.0/sqrt2;
        sqrtD = np.sqrt(d);
        sqrtDi = 1.0/sqrtD;
        sbx  = np.sign(bx);
        og1  = 1.0/(Gamma - 1);

        if EquationOfState == 0 :
            p = (Gamma -1.0 ) * (eng - 0.5* d * (vx*vx + vy*vy + vz * vz ) - 0.5*(bx*bx + by*by + bz*bz));
        else:
            p = IsothermalSoundSpeed*IsothermalSoundSpeed * d;

        #compute wave speeds
        if EquationOfState == 0:
            aa = np.sqrt( 1.0*Gamma* p/d )
        else:
            aa = IsothermalSoundSpeed;

        cs =np.sqrt( 0.5*( aa*aa + 2*bp/d -np.sqrt( pow( (aa*aa + 2*bp/d ),2) - 4* aa*aa*bx*bx/d ) ) );
        ca =np.sqrt( bx*bx/d ); 
        cf =np.sqrt( 0.5*( aa*aa + 2*bp/d +np.sqrt( pow( (aa*aa + 2*bp/d ),2) - 4* aa*aa*bx*bx/d ) ) );
        #compute ancilary values
        #The normalization of alph_f may change. This normalization uses Ryu
        #& Jones, but Balsara may be more robust.

        betay=np.zeros_like(by)+sqrt2i
        betaz=np.zeros_like(by)+sqrt2i
        alph_f=np.zeros_like(by)+1
        alph_s=np.zeros_like(by)+1
        #non-degenerate points
        non_deg = np.logical_or( by != 0.0 , bz != 0.0 )
        bt =np.zeros_like(by)
        if is_iterable(by): #hasattr(by,'__getitem__'):
            bt[non_deg] = np.sqrt( by[non_deg]*by[non_deg] + bz[non_deg]*bz[non_deg] );
            betay[non_deg] = by[non_deg]/bt[non_deg];
            betaz[non_deg] = bz[non_deg]/bt[non_deg];
            alph_f[non_deg] = np.sqrt( (cf[non_deg]*cf[non_deg]-ca[non_deg]*ca[non_deg])/(cf[non_deg]*cf[non_deg]-cs[non_deg]*cs[non_deg]) );
            alph_s[non_deg] = np.sqrt( (cf[non_deg]*cf[non_deg]-aa[non_deg]*aa[non_deg])/(cf[non_deg]*cf[non_deg]-cs[non_deg]*cs[non_deg]) );
        else:
            if non_deg:
                bt = np.sqrt( by*by + bz*bz );
                betay = by/bt;
                betaz = bz/bt;
                alph_f = np.sqrt( (cf*cf-ca*ca)/(cf*cf-cs*cs) );
                alph_s = np.sqrt( (cf*cf-aa*aa)/(cf*cf-cs*cs) );
      
#    if by == 0.0 and bz == 0.0:
#        betay = sqrt2i
#        betaz = sqrt2i
#        alph_f = 1
#        alph_s = 1
#    else:
#        bt = np.sqrt( by*by + bz*bz );
#        betay = by/bt;
#        betaz = bz/bt;
#        alph_f = np.sqrt( (cf*cf-ca*ca)/(cf*cf-cs*cs) );
#        alph_s = np.sqrt( (cf*cf-aa*aa)/(cf*cf-cs*cs) );
      
        self.alph_f = alph_f
        self.alph_s = alph_s
        self.betay = betay
        self.betaz = betaz
        #the vectors
        right={}
        for w in ['f-','a-','s-','c','f+','a+','s+']:
            right[w] = {}
            for f in ['d','px','py','pz','hx','hy','hz','e']:
                right[w][f]=np.zeros_like(by)
      
        #fast, left
        right['f-']['d'] = alph_f
        if EquationOfState == 0 :
            right['f-']['e'] =  alph_f*0.5*(vx*vx+vy*vy+vz*vz) 
            right['f-']['e'] += alph_f*cf*cf*og1 - alph_f*cf*vx + alph_s*ca*sbx*(betay*vy + betaz*vz) 
            right['f-']['e'] += (Gamma-2)*og1*alph_f*(cf*cf-aa*aa)
        right['f-']['px'] = alph_f*(vx - cf);
        right['f-']['py'] = alph_f*vy + alph_s*betay*ca*sbx;
        right['f-']['pz'] = alph_f*vz + alph_s*betaz*ca*sbx;
        right['f-']['hy'] = alph_s*betay*cf*sqrtDi;
        right['f-']['hz'] = alph_s*betaz*cf*sqrtDi;


        #alfven][left
        #right['a-']['d'] = 0; already zero
        if EquationOfState == 0:
            right['a-']['e'] = 1*(betaz*vy - betay*vz)*sbx;
        #right['a-']['px'] = 0; already zero.
        right['a-']['py'] =  1*betaz*sbx;
        right['a-']['pz'] = -1*betay*sbx;
        right['a-']['hy'] = betaz*sqrtDi;
        right['a-']['hz'] = -betay*sqrtDi;

        #alfven,right
        #right['a+']['d'] = 0.; already zero
        if EquationOfState == 0.:
            right['a+']['e'] = -1*(betaz*vy - betay*vz)*sbx;
        #right['a+']['px'] = 0.; already zero.
        right['a+']['py'] = -1*betaz*sbx;
        right['a+']['pz'] = +1*betay*sbx;
        right['a+']['hy'] = betaz*sqrtDi;
        right['a+']['hz'] = -betay*sqrtDi;

      
        #slow,left
        right['s-']['d'] = alph_s;
        if EquationOfState == 0:
            right['s-']['e'] = alph_s*0.5*(vx*vx+vy*vy+vz*vz) + \
                alph_s*cs*cs*og1 - alph_s*cs*vx - alph_f*aa*sbx*(betay*vy + betaz*vz) +\
                (Gamma-2)*og1*alph_s*(cs*cs - aa*aa );
        right['s-']['px'] = alph_s*(vx-cs);
        right['s-']['py'] = alph_s*vy - alph_f*betay*aa*sbx;
        right['s-']['pz'] = alph_s*vz - alph_f*betaz*aa*sbx;
        right['s-']['hy'] = -alph_f*betay*aa*aa*sqrtDi/cf;
        right['s-']['hz'] = -alph_f*betaz*aa*aa*sqrtDi/cf;
      
        #entropy (no entropy wave in isothermal MHD.)(Or hydro,for that matter)
        if EquationOfState == 0:
            right['c']['d'] = 1;
            right['c']['e'] = 0.5*(vx*vx+vy*vy+vz*vz);
            right['c']['px'] = vx;
            right['c']['py'] = vy;
            right['c']['pz'] = vz;
            right['c']['hy'] = 0;
            right['c']['hz'] = 0;
      
        #slow,right
        right['s+']['d'] = alph_s;
        if EquationOfState  == 0:
            right['s+']['e'] = alph_s*0.5*(vx*vx+vy*vy+vz*vz) + alph_s*cs*cs*og1 + \
                alph_s*cs*vx + alph_f*aa*sbx*(betay*vy + betaz*vz) + (Gamma-2)*og1*alph_s*(cs*cs - aa*aa );
        right['s+']['px'] = alph_s*(vx+cs);
        right['s+']['py'] = alph_s*vy + alph_f*betay*aa*sbx;
        right['s+']['pz'] = alph_s*vz + alph_f*betaz*aa*sbx;
        right['s+']['hy'] = -alph_f*betay*aa*aa*sqrtDi/cf;
        right['s+']['hz'] = -alph_f*betaz*aa*aa*sqrtDi/cf;
      


        #fast, right
        right['f+']['d'] = alph_f;
        if EquationOfState  == 0:
            right['f+']['e'] = alph_f*0.5*(vx*vx+vy*vy+vz*vz) + \
            alph_f*cf*cf*og1 + alph_f*cf*vx - alph_s*ca*sbx*(betay*vy + betaz*vz) + \
            (Gamma-2)*og1*alph_f*(cf*cf-aa*aa);
        right['f+']['px'] = alph_f*(vx + cf);
        right['f+']['py'] = alph_f*vy - alph_s*betay*ca*sbx;
        right['f+']['pz'] = alph_f*vz - alph_s*betaz*ca*sbx;
        right['f+']['hy'] = alph_s*betay*cf*sqrtDi;
        right['f+']['hz'] = alph_s*betaz*cf*sqrtDi;
        return None, right
    def speeds(self,d,vx,vy,vz,bx,by,bz,eng=None,p=None, EquationOfState = 0, Gamma=5./3, IsothermalSoundSpeed=1):

        #float d, float vx, float vy, float vz, 
        #float bx, float by, float bz, float eng, 
        #float right[][7])

        #Eigen vectors taken from Ryu & Jones, ApJ 442:228-258, 1995
        #Normalization starting on p. 231.

        bp = 0.5*(bx*bx + by*by + bz*bz );
        sqrt2 = np.sqrt(2.0);
        sqrt2i= 1.0/sqrt2;
        sqrtD = np.sqrt(d);
        sqrtDi = 1.0/sqrtD;
        sbx  = np.sign(bx);
        og1  = 1.0/(Gamma - 1);
        if p is not None:
            Egas=p/((Gamma-1)*d)
            eng = 0.5*(vx*vx+vy*vy+vz*vz)+0.5*(bx*bx+by*by+bz*bz)/d + Egas


        if EquationOfState == 0 :
            p = (Gamma -1.0 ) * (eng - 0.5* d * (vx*vx + vy*vy + vz * vz ) - 0.5*(bx*bx + by*by + bz*bz));
        else:
            p = IsothermalSoundSpeed*IsothermalSoundSpeed * d;

        #compute wave speeds
        if EquationOfState == 0:
            aa = np.sqrt( 1.0*Gamma* p/d )
        else:
            aa = IsothermalSoundSpeed;

        cs =np.sqrt( 0.5*( aa*aa + 2*bp/d -np.sqrt( pow( (aa*aa + 2*bp/d ),2) - 4* aa*aa*bx*bx/d ) ) );
        ca =np.sqrt( bx*bx/d ); 
        cf =np.sqrt( 0.5*( aa*aa + 2*bp/d +np.sqrt( pow( (aa*aa + 2*bp/d ),2) - 4* aa*aa*bx*bx/d ) ) );
        return {'cf':cf,'cs':cs,'ca':ca,'aa':aa}


