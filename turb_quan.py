import taxi
from GL import *
import p49_fields
import p49_labels
import p49_QU2EB
import xtra_energy_fields
import astropy.io.fits as pyfits
import glob
import os
import numpy as np
import time
import glob
import h5py
import re
nar = np.array
import xtra_energy_fields
reload(xtra_energy_fields)
import cmbtools_handler as cmbtools
#if '/home/dcollins/repos/p49c/p49_eigenmodes' not in sys.path:
#    sys.path.append('/home/dcollins/repos/p49c/p49_eigenmodes')
#import p49_print_tools
import matplotlib.colors as colors
def plotter2(arrays,fname,norm=None,labs=['Q','U','E','B'],axis_labels=None,npx=2,
             share=True,zmin=None,**args):
    np = len(arrays)
    npy = np//npx
    if np%npx: npy+=1
    fig, axes = plt.subplots(npx,npy,sharex=share,sharey=share)
    print("NPX %d NPY %d"%(npx,npy))
    max_val = max([arr.max() for arr in arrays])
    min_all = min([arr.min() for arr in arrays])
    if zmin is not None:
        min_all = zmin
    if norm is 'positive':
        min_val = min_all
        args['norm']=colors.LogNorm(vmin=min_val,vmax=max_val)
        print("WTF",min_val)
    elif norm is 'symlog' :
        min_val = min_all
        args['norm']=colors.SymLogNorm(linthresh=1e-10,vmin=min_val,vmax=max_val)
    elif norm is 'ind':
        args['norm']=None
    else:
        min_val = min_all
        args['norm']=colors.Normalize(vmin=min_val,vmax=max_val)

    args['interpolation']='nearest'
    args['origin']='lower'
    oot=[]
    for n,arr in enumerate(arrays):
        ipy = n//npx
        ipx = n%npx
        if npx == 1:
            this_ax = axes
        elif npy == 1:
            this_ax = axes[ipx]
        else:
            this_ax = axes[ipx][ipy]
        oot.append(this_ax.imshow(arr,**args))
        if norm is 'ind': cb=fig.colorbar(oot[-1],ax=this_ax)
        this_ax.set_title(labs[n])
#   oot.append(axes[1][0].imshow(U,**args))
#   if norm is 'ind': cb=fig.colorbar(oot[-1],ax=axes[1][0])
#   axes[1][0].set_title(labs[1])
#   oot.append(axes[0][1].imshow(E,**args))
#   if norm is 'ind': cb=fig.colorbar(oot[-1],ax=axes[0][1])

#   axes[0][1].set_title(labs[2])
#   oot.append(axes[1][1].imshow(B,**args))
#   if norm is 'ind': cb=fig.colorbar(oot[-1],ax=axes[1][1])
#   axes[1][1].set_title(labs[3])
    if norm is not 'ind':
        cb=fig.colorbar(oot[-1],ax=axes)
        if norm is 'positive':
            cb.cmap.set_under('w')
    if axis_labels:
        for n in range(np):
            axes[n//2][n%2].set_xlabel(axis_labels[0])
            axes[n//2][n%2].set_ylabel(axis_labels[1])
    fig.savefig(fname)
    print(fname)
    plt.close(fig)

def plotter(Q,U,E,B,fname,norm=None,labs=['Q','U','E','B'],axis_labels=None,**args):
    fig, axes = plt.subplots(2,2,sharex=True,sharey=True)
    max_val = max([Q.max(),U.max(),E.max(),B.max()])
    if norm is 'positive':
        min_val = 1e-15
        args['norm']=colors.LogNorm(vmin=min_val,vmax=max_val)
    elif norm is 'symlog' :
        min_val = min([Q.min(),U.min(),E.min(),B.min()])
        args['norm']=colors.SymLogNorm(linthresh=1e-10,vmin=min_val,vmax=max_val)
    elif norm is 'ind':
        args['norm']=None
    else:
        min_val = min([Q.min(),U.min(),E.min(),B.min()])
        args['norm']=colors.Normalize(vmin=min_val,vmax=max_val)

    args['interpolation']='nearest'
    args['origin']='lower'
    oot=[]
    oot.append(axes[0][0].imshow(Q,**args))
    if norm is 'ind': cb=fig.colorbar(oot[-1],ax=axes[0][0])
    axes[0][0].set_title(labs[0])
    oot.append(axes[1][0].imshow(U,**args))
    if norm is 'ind': cb=fig.colorbar(oot[-1],ax=axes[1][0])
    axes[1][0].set_title(labs[1])
    oot.append(axes[0][1].imshow(E,**args))
    if norm is 'ind': cb=fig.colorbar(oot[-1],ax=axes[0][1])

    axes[0][1].set_title(labs[2])
    oot.append(axes[1][1].imshow(B,**args))
    if norm is 'ind': cb=fig.colorbar(oot[-1],ax=axes[1][1])
    axes[1][1].set_title(labs[3])
    if norm is not 'ind':
        cb=fig.colorbar(oot[-1],ax=axes)
        if norm is 'positive':
            cb.cmap.set_under('w')
    if axis_labels:
        for n in range(4):
            axes[n//2][n%2].set_xlabel(axis_labels[0])
            axes[n//2][n%2].set_ylabel(axis_labels[1])
    fig.savefig(fname)
    plt.close(fig)

frbname = p49_QU2EB.frbname
#fleet = [aw11]
class quan_box():
    def __init__(self,car=None, plot_format='png'):
        
        self.car = car
        if car is not None:
            self.name = car.name
        else:
            self.name='NAME'
        self.keys=['ex','ey','ez','vx','vy','vz','px','py','pz','mach']
        self.magkeys=['Bx','By','Bz','bx','by','bz','bx2','by2','bz2','beta','AlfMach','AlfvenSpeed']
        self.all_fields = ['vx','vy','vz','mach','px','py','pz','ex','ey','ez','t','bx','by','bz','bx2','by2','bz2']
        self.all_fields +=['Bx','By','Bz','Bfield_strength','AlfMach','beta','AlfvenSpeed','frames']
        self.all_fields +=['ke_tot','ke_rel','grav_pot','grav_pot_2','gas_work', 'density','density2']
        self.potential_written = False
        self.clobber=False
        self.stuff={}
        self.plot_format=plot_format
        #self.EBSlopePower={}

    def merge(self,dict2):
        dict1 = self.stuff
        #dict2 = other_quan.stuff
        if 'frames' in dict1 and 'frames' in dict2:
            frames1 = dict1['frames']
            frames2 = dict2['frames']
            for i, frame in enumerate(frames2):
                if frame not in frames1:
                    for key in dict1:
                        if key in [ 'EB', 'grav_pot_2', 'tdyn']:
                            continue
                        if key in ['grav_pot'] and not self.potential_written:
                            continue
                        dict1[key].append(dict2[key][i])

        if 'EB' in dict2:
            if 'EB' not in dict1:
                dict1['EB'] = {}
            for frame in dict2['EB']:
                if frame not in dict1['EB']:
                    dict1['EB'][frame] = dict2['EB'][frame]
    def dump(self,h5name=None):
        if h5name is None:
            h5_name = 'quan_box_%s.h5'%self.car.outname
        #if os.path.exists(pickle_name):
        #    #other_pickle = fPickle.load(pickle_name)
        #    self.merge(other_pickle)
        fptr = h5py.File(h5_name,'a')
        frames_write = np.ones_like(self.stuff.get('frames',[]),dtype=bool)
        ebframes_write = np.ones_like(self.stuff.get('EBcycles',[]),dtype=bool)
        if 'frames' in fptr and 'frames' in self.stuff:
            frames_write = np.zeros_like(self.stuff['frames'],dtype=bool)
            for i, n in enumerate(self.stuff['frames']):
                if n not in fptr['frames']:
                    frames_write[i] = True
                else:
                    frames_write[i] = False

        if 'EBcycles' in fptr and 'EBcycles' in self.stuff:
            for i, n in enumerate(self.stuff['EBcycles']):
                if n not in fptr['EBcycles']:
                    ebframes_write[i] = True
                else:
                    ebframes_write[i] = False

        for key in self.stuff:
            if key in ['tdyn']:
                continue
            if len(self.stuff[key]) == 0:
                continue
            if key[0:4] in ["Eamp", "Bamp", "Eslo", "Bslo", "EBcy"]:
                to_write = ebframes_write
            else:
                to_write = frames_write
            newsize = (to_write.sum(),)

            if key not in fptr:
                dset = fptr.create_dataset(key,newsize,maxshape=(None,))
                dset[:]=np.array(self.stuff[key])[to_write]
            else:
                dset = fptr[key]
                size1 = dset.size
                dset.resize( (size1 + newsize[0],))
                dset[size1:] = np.array(self.stuff[key])[to_write]

        fptr.close()
            

    def load(self, h5_name=None):
        if h5_name is None:
            h5_name = 'quan_box_%s.h5'%self.car.outname
        if len(glob.glob(h5_name)):
            fptr = h5py.File(h5_name,'r')
            try:
                for k in fptr:
                    self.stuff[k]=copy.copy(fptr[k].value)
            except:
                pdb.set_trace()
                raise
            finally:
                fptr.close()
        else:
            print("NO SUCH FILE", h5_name)

    def __getitem__(self,key):
        return self.__dict__[key]
    def show(self,format="%10.2e"):
        N = len(self.frames)
        keys_to_print = self.keys
        if hasattr(self,'Bx'):
            keys_to_print += self.magkeys
        for key in keys_to_print:
            print("%5s"%key,)
            print(format*N%tuple(self[key]))



    def plot_eb_vs_stuff(self,eb_quantity, axis, other_quantity):
        quad_frames = self.stuff['frames']
        eb_frames = list(self.stuff['EB'].keys())
        all_frames = np.unique(nar(quad_frames + eb_frames))


    def EBall(self,frames=None):
        if 'EBcycles' not in self.stuff:
            self.stuff['EBcycles']=[]
            for ax in 'xyz':
                self.stuff['Eamp_%s'%ax]=[]
                self.stuff['Bamp_%s'%ax]=[]
                self.stuff['Eslope_%s'%ax]=[]
                self.stuff['Bslope_%s'%ax]=[]
        if frames is None:
            frames = self.car.return_frames()
        for frame in frames:
            if frame not in self.stuff['EBcycles']:
                self.make_frbs(frame)
                self.QUEB(frame)
                self.EBslopes(frame)
        self.dump()
    def EBslopes(self,frame):
        EBSlopePower=p49_QU2EB.slopes_powers(self.car,frame, plot_format=self.plot_format)
        if 'EBcycles' not in self.stuff:
            self.stuff['EBcycles']=[]
        self.stuff['EBcycles'].append(frame)
        for ax in 'xyz':
            self.stuff['Eamp_%s'%ax  ].append(EBSlopePower['Eamp'][ax])
            self.stuff['Bamp_%s'%ax  ].append(EBSlopePower['Bamp'][ax])
            self.stuff['Eslope_%s'%ax].append(EBSlopePower['Eslope'][ax])
            self.stuff['Bslope_%s'%ax].append(EBSlopePower['Bslope'][ax])
    def GetQUEB(self,frame):

        frb_dir = "%s/%s/"%(self.car.directory,frbname)
        Qlist = glob.glob(frb_dir+'/DD%04d_Q[xyz]*.fits'%frame)
        Ulist = []
        # write the output near the input
        self.QUEBarr = {'Q':{}, 'U':{}, 'E':{}, 'B':{}}
        for Qfile in Qlist:
            mo = re.match('(.*/DD[0-9]{4}_)Q([xyz].*)(.fits)',Qfile)
            Ufile = mo.group(1)+'U'+mo.group(2)+'.fits'
            Ulist.append(Ufile)
            mo = re.match('(.*/DD[0-9]{4}_)Q([xyz].*)(.fits)',Qfile)
            outroot = mo.group(1)
            outsuf = mo.group(2)
            Efile = outroot+'E'+outsuf+'.fits'
            Bfile = outroot+'B'+outsuf+'.fits'
            Clfile = outroot+'Cl'+outsuf+'.dat'
            #self.QUEBarr['Q'][Qfile] = np.array(pyfits.open(Qfile)[0].data,dtype=np.double)
            #self.QUEBarr['U'][Ufile] = np.array(pyfits.open(Ufile)[0].data,dtype=np.double)
            #self.QUEBarr['E'][Efile] = np.array(pyfits.open(Efile)[0].data,dtype=np.double)
            #self.QUEBarr['B'][Bfile] = np.array(pyfits.open(Bfile)[0].data,dtype=np.double)
            self.QUEBarr['Q'][outsuf] = np.array(pyfits.open(Qfile)[0].data,dtype=np.double)
            self.QUEBarr['U'][outsuf] = np.array(pyfits.open(Ufile)[0].data,dtype=np.double)
            self.QUEBarr['E'][outsuf] = np.array(pyfits.open(Efile)[0].data,dtype=np.double)
            self.QUEBarr['B'][outsuf] = np.array(pyfits.open(Bfile)[0].data,dtype=np.double)

    def PlotQUEBharm(self,frame, style=1):
        self.GetQUEB(frame)

        if style == 0:
            for field in 'QUEB':
                for f in qb.QUEBarr[field]:
                    oot = "%s_%s"%(selfcar.outname,f.split('/')[-1])
                    plt.clf()
                    #Qharm = cmbtools.map2harm(Q,Delta)
                    #Uharm = cmbtools.map2harm(U,Delta)
                    pdb.set_trace()
                    f = plt.imshow(field_hat,interpolation='nearest',origin='lower')
                    plt.colorbar(f)
                    print("did a color")
                    plt.title(oot)
                    plt.savefig(oot+".png")
                    print(oot)
        else:
            for ax in 'xyz':
                outname = 'P49_QUEB_harm_4p_%s_%s.%s'%(self.car.outname,ax,self.plot_format)
                hats = {}
                for field in 'QUEB':
                    arr =  self.QUEBarr[field][ax]
                    Delta = np.array([5*np.pi/180.]*2)/np.array(arr.shape)
                    hats[field] = np.abs(cmbtools.map2harm(arr,Delta))
                    #print(p49_print_tools.nz(hats[field]))
                #print(p49_print_tools.nonzero(hats['Q'][:10,:10]))
                #print(p49_print_tools.nonzero(hats['U'][:10,:10]))
                #print(p49_print_tools.nonzero(hats['E'][:10,:10]))
                #print(p49_print_tools.nonzero(hats['B'][:10,:10]))


                plotter(hats['Q'],
                        hats['U'],
                        hats['E'],
                        hats['B'],
                        outname,norm='positive')
                print("Save %s"%outname)

    def PlotQUEB(self,frame, style=1):
        self.GetQUEB(frame)
        if style == 0:
            for field in 'QUEB':
                for f in qb.QUEBarr[field]:
                    oot = "%s_%s"%(selfcar.outname,f.split('/')[-1])
                    plt.clf()
                    f = plt.imshow(qb.QUEBarr[field][f],interpolation='nearest',origin='lower')
                    plt.colorbar(f)
                    print("did a color")
                    plt.title(oot)
                    plt.savefig(oot+".png")
                    print(oot)
        else:
            for ax in 'xyz':
                outname = 'P49_QUEB_4p_%s_%s.%s'%(self.car.outname,ax,self.plot_format)
                plotter(self.QUEBarr['Q'][ax],
                        self.QUEBarr['U'][ax],
                        self.QUEBarr['E'][ax],
                        self.QUEBarr['B'][ax],
                        outname,norm='ind')
                print("Save %s"%outname)
                fig,axobj=plt.subplots(2)
                for field in 'QUEB':
                    axobj[0].hist(self.QUEBarr[field][ax].flatten(),histtype='step',label=field)
                    xxx=self.QUEBarr[field][ax].flatten()-np.mean(self.QUEBarr[field][ax])
                    xxx=np.abs(xxx)
                    axobj[1].hist(xxx,histtype='step',label=field)
                axobj[0].legend(loc=0)
                axobj[1].legend(loc=0)
                axobj[1].set_xscale('log')
                outname = 'P49_QUEB_hist_%s_%s.%s'%(self.car.outname,ax,self.plot_format)
                fig.savefig(outname)
                plt.close(fig)



    def QUEB(self, frame):
        #ds = self.car.load(frame)
        frb_dir = "%s/%s/"%(self.car.directory,frbname)
        p49_QU2EB.QU2EB(frb_dir,frame)
#        if frames is None:
#            frames = self.car.return_frames()
#        for frame in frames:
#            self.make_frbs(frame)
#            self.fit_slopes()
#            self.plot_eebb()
    def make_frbs(self,frame, axes=['x','y','z'], ds=None):
        fields=[]
        for axis in axes:
          n0=1; p=1 #n0 in [19,39,1945] and p=0
          #fields.append( (axis,'Q%s_n0-%04d_p-%d'%(axis,n0,p))   )
          #fields.append( (axis,'U%s_n0-%04d_p-%d'%(axis,n0,p))   )
          fields.append( (axis,'Q%s'%(axis))   )
          fields.append( (axis,'U%s'%(axis))   )
          fields.append( (axis,'density') )

        for axis, field in fields :
            outputdir = "%s/%s/"%(self.car.directory,frbname)
            if not os.access(outputdir, os.F_OK):
                os.mkdir(outputdir)
            #Hm.  Q and U have the name in the field, but others don't.
            if field[0] in 'QU' and field[1] in 'xyz':
                field_name = field
            else:
                field_name = field + "_"+axis
            outfile = outputdir+"/DD%.4d_%s.fits" %(frame,field_name)
            #move this in the 'make' conditional?
            #if ds is None:
            #    ds = self.car.load(frame)
            #    res = ds.parameters['TopGridDimensions'][2 + ord('x') - ord(axis)] # zyx order
            if os.access(outfile, os.F_OK) and not self.clobber:
                print("FRB exists: %s"%outfile)
            else:
                if ds is None:
                    ds = self.car.load(frame)
                    res = ds.parameters['TopGridDimensions'][2 + ord('x') - ord(axis)] # zyx order
                print("FRB being produced: %s"%outfile)
                res = ds.parameters['TopGridDimensions'][0] #2 + ord('x') - ord(axis)]
                proj = ds.proj(field,axis)
                frb = proj.to_frb(1,res)
                hdu = pyfits.PrimaryHDU(frb[field])
                hdulist = pyfits.HDUList([hdu])
                hdulist.writeto(outfile,clobber=True)
                print("wrote", outfile)



        #car.plot()
    def __call__(self, ds=None, tdyn=1, frames=None):
        if frames is None:
            frames = self.car.return_frames()
        for frame in frames:
            self.quadratic(self.car,frame,tdyn,ds=ds)
        #car.plot()
    def quadratic(self,car, frame, tdyn=1, ds=None):

        for field in self.all_fields:
            if field not in self.stuff:
                self.stuff[field] = []

        if frame in self.stuff['frames'] and self.clobber == False:
            return
        if 'tdyn' in self.stuff:
            self.tdyn = self.stuff['tdyn']
        else:
            self.tdyn=tdyn
        self.stuff['tdyn']=self.tdyn


        if ds is None:
            print(car.name)
            car.region_type='all'
            ds=car.load(frame)
            UseMHD = car.ds['HydroMethod'] in [4,6,9000]
            this_time=car.ds['InitialTime']/(tdyn) 
        else:
            UseMHD = True
            this_time=frame

        self.potential_written = False
        try:
            if car.ds['SelfGravity']:
                car.ds.create_field_info()
                self.potential_written = 'PotentialField' in [k[1] for k in list(car.ds.field_info.keys())]
        except:
            self.potential_written =False

        if frame in self.stuff['frames'] and self.clobber == False:
            return
        self.stuff['frames'].append(frame)
        self.stuff['t'].append(this_time)
        if ds is not None:
            reg=ds.all_data()
        else:
            reg = car.get_region(frame)
            ds = car.ds
        xtra_energy_fields.dave_add_field(ds) #adds many fields.
        total_volume = reg['cell_volume'].sum()
        volume = reg['cell_volume']
        print("volume", volume)
        #self.stuff['mass'].append( (reg.quantities['WeightedAverageQuantity']('density', 'cell_volume')*total_volume).in_units('code_mass'))
        self.stuff['density'].append( reg.quantities['WeightedAverageQuantity']('density','cell_volume').in_units('code_density').v  )
        self.stuff['density2'].append( ((( reg['density'].in_units('code_density').v-self.stuff['density'][-1])**2*reg['cell_volume']).sum()/total_volume)  )
        self.stuff['mach'].append( np.sqrt(reg.quantities['WeightedAverageQuantity']('mean_square_velocity','cell_volume').in_units('code_velocity**2').v ) )
        self.stuff['vx'].append(reg.quantities['WeightedAverageQuantity']('velocity_x','cell_volume').in_units('code_velocity').v)
        self.stuff['vy'].append(reg.quantities['WeightedAverageQuantity']('velocity_y','cell_volume').in_units('code_velocity').v)
        self.stuff['vz'].append(reg.quantities['WeightedAverageQuantity']('velocity_z','cell_volume').in_units('code_velocity').v)
        self.stuff['px'].append(reg.quantities['WeightedAverageQuantity']('momentum_x','cell_volume').in_units('code_density*code_velocity').v)
        self.stuff['py'].append(reg.quantities['WeightedAverageQuantity']('momentum_y','cell_volume').in_units('code_density*code_velocity').v)
        self.stuff['pz'].append(reg.quantities['WeightedAverageQuantity']('momentum_z','cell_volume').in_units('code_density*code_velocity').v)
        self.stuff['ex'].append(reg.quantities['WeightedAverageQuantity']('eng_x','cell_volume').in_units('code_density*code_velocity**2').v)
        self.stuff['ey'].append(reg.quantities['WeightedAverageQuantity']('eng_y','cell_volume').in_units('code_density*code_velocity**2').v)
        self.stuff['ez'].append(reg.quantities['WeightedAverageQuantity']('eng_z','cell_volume').in_units('code_density*code_velocity**2').v)
        self.stuff['ke_tot'].append(reg.quantities['WeightedAverageQuantity']('kinetic_energy','cell_volume').in_units('code_density*code_velocity**2').v)
        reg.set_field_parameter('bulk_velocity',ds.arr([self.stuff['vx'][-1],self.stuff['vy'][-1], self.stuff['vz'][-1]],'code_velocity'))
        self.stuff['ke_rel'].append(reg.quantities['WeightedAverageQuantity']('rel_kinetic_energy','cell_volume').in_units('code_density*code_velocity**2').v)
        if self.potential_written:
            self.stuff['grav_pot'].append(reg.quantities['WeightedAverageQuantity']('grav_pot','cell_volume').in_units('code_density*code_velocity**2').v)
        self.stuff['gas_work'].append(reg.quantities['WeightedAverageQuantity']('gas_work','cell_volume').in_units('code_density*code_velocity**2').v)

        if UseMHD:
            self.stuff['Bx'].append( (reg['Bx']*volume).sum()/total_volume)
            self.stuff['By'].append( (reg['By']*volume).sum()/total_volume)
            self.stuff['Bz'].append( (reg['Bz']*volume).sum()/total_volume)

            self.stuff['bx'].append( ( (reg['Bx']-self.stuff['Bx'][-1])*volume).sum()/total_volume)
            self.stuff['by'].append( ( (reg['By']-self.stuff['By'][-1])*volume).sum()/total_volume)
            self.stuff['bz'].append( ( (reg['Bz']-self.stuff['Bz'][-1])*volume).sum()/total_volume)

            self.stuff['bx2'].append( np.sqrt(( (reg['Bx']-self.stuff['Bx'][-1])**2*volume).sum()/total_volume) )
            self.stuff['by2'].append( np.sqrt(( (reg['By']-self.stuff['By'][-1])**2*volume).sum()/total_volume) )
            self.stuff['bz2'].append( np.sqrt(( (reg['Bz']-self.stuff['Bz'][-1])**2*volume).sum()/total_volume) )

            self.stuff['Bfield_strength'].append( (reg['magnetic_field_strength']*volume).sum()/total_volume)
            self.stuff['AlfvenSpeed'].append( (volume*reg['magnetic_field_strength']/np.sqrt(np.pi*4*reg['density']) ).sum()/total_volume)

            self.stuff['AlfMach'].append( self.stuff['mach'][-1]/self.stuff['AlfvenSpeed'][-1])
            self.stuff['beta'].append( 2*(self.stuff['AlfMach'][-1]/self.stuff['mach'][-1])**2 )
        #picklename = "turb_quan_temp_%s.pickle"%car.name
        self.dump()
        #self.plot()





    def plot(self, HydroMethod = None):
        print( "Hydro Method", HydroMethod)
        def nom_string(val):
            if val > 10 or val < 0.1:
                return "%0.1e"%val
            else:
                return "%0.1f"%val
        print("here's the thing")

        car = self.car
        if HydroMethod is None:
            ds=car.load()
            HydroMethod = car.ds['HydroMethod']
        tn = p49_labels.nominal.get(car.name,None)
        times = nar(self.stuff['t'])
        ar = np.argsort(times)
        times=times[ar]
        for thing in self.stuff:
            if thing in ['tdyn', 'EB']:
                continue
            if len(self.stuff[thing]) == 0:
                continue
            self.stuff[thing] = nar(self.stuff[thing])[ar]
        time_label  =r'$t [\rm{code}]$'
        if tn:
            time_label += r' $(t_{\rm{cross}}= %s)$'% nom_string(0.5/tn['mach'])

        sqrtfourpi=np.sqrt(4*np.pi)

        if HydroMethod in [4,6]:
            plt.clf()
            n_points = len(times)
            my_ones = np.ones(n_points)
            
            if tn is not None:
                nom_val = tn['field_cgs']/sqrtfourpi
                st = "nominal "+ (nom_string(nom_val))
                plt.plot(times, my_ones*nom_val,label=st,c=[0.5]*4)
            plt.plot(times,self.stuff['Bx'],label='Bx')
            plt.plot(times,self.stuff['By'],label='By')
            plt.plot(times,self.stuff['Bz'],label='Bz')
            plt.plot(times,self.stuff['bx'],label='bx')
            plt.plot(times,self.stuff['by'],label='by')
            plt.plot(times,self.stuff['bz'],label='bz')
            plt.plot(times,self.stuff['bx2'],label='bx rms')
            plt.plot(times,self.stuff['by2'],label='by rms')
            plt.plot(times,self.stuff['bz2'],label='bz rms')
            plt.legend(loc=0)
            plt.ylabel('MagneticField'); plt.xlabel(time_label)
            outname = '%s_quan_field_strength.%s'%(car.name,self.plot_format)
            plt.savefig(outname)
            print(outname)

            plt.clf()
            c='r'
            if tn is not None:
                plt.plot(times, my_ones*tn['AlfMach'],label='nominal',c=c,linestyle='--')
            plt.plot(times, self.stuff['AlfMach'], label='MA',c=c )
            c='g'
            if tn is not None:
                plt.plot(times, my_ones*tn['mach'],c=c,linestyle='--')
            plt.plot(times, self.stuff['mach'], label="M", c='g')
            c='c'
            if tn is not None:
                plt.plot(times, my_ones*(tn['mach']/tn['AlfMach']),c=c,linestyle='--')
            plt.plot(times, self.stuff['AlfvenSpeed'],label='Va',c=c)
            c='b'
            if tn is not None:
                plt.plot(times, my_ones*(10**tn['logbeta']),c=c,linestyle='--')
            plt.plot(times, self.stuff['beta'], label='beta',c=c)
            plt.ylabel('Dimensionless'); plt.xlabel(time_label)
            plt.legend(loc=0)
            outname = '%s_quan_MaM.%s'%(car.name,self.plot_format)
            plt.savefig(outname)
            print(outname)
            plt.yscale('log')
            outname = '%s_quan_MaM_log.%s'%(car.name,self.plot_format)
            plt.savefig(outname)
            print(outname)

        plt.clf()
        plt.plot(times,self.stuff['ex'],label='ex')
        plt.plot(times,self.stuff['ey'],label='ey')
        plt.plot(times,self.stuff['ez'],label='ez')
        plt.ylabel('Partial Energies'); plt.xlabel(time_label)
        plt.legend(loc=0)
        outname = '%s_quan_eng.%s'%(car.name,self.plot_format)
        plt.savefig(outname)
        print(outname)

        plt.clf()
        plt.plot(times,self.stuff['px'],label='px')
        plt.plot(times,self.stuff['py'],label='py')
        plt.plot(times,self.stuff['pz'],label='pz')
        plt.legend(loc=0)
        outname = '%s_quan_mom.%s'%(car.name,self.plot_format)
        plt.ylabel('Momentum'); plt.xlabel(time_label)
        plt.savefig(outname)
        print(outname)

        plt.clf()
        plt.plot(times,self.stuff['vx'],label='vx',c='r')
        plt.plot(times,self.stuff['vy'],label='vy',c='g')
        plt.plot(times,self.stuff['vz'],label='vz',c='b')
        c='k'
        if tn is not None:
            plt.plot(times, my_ones*tn['mach'],c=c,linestyle='--')
        plt.plot(times,self.stuff['mach'],label='mach',c=c)
        plt.ylabel('velocities'); plt.xlabel(time_label)
        plt.legend(loc=0)
        outname = '%s_quan_vel.%s'%(car.name,self.plot_format)
        plt.savefig(outname)
        print(outname)

        plt.clf()
        plt.plot(times,self.stuff['ke_tot'],label='ke_tot')
        plt.plot(times,self.stuff['ke_rel'],label='ke_rel')
        if self.potential_written: 
            plt.plot(times,self.stuff['grav_pot'],label='grav_pot')
        plt.plot(times,self.stuff['gas_work'],label='gas_work')
        plt.ylabel('Energies'); plt.xlabel(time_label)
        plt.legend(loc=0)
        outname = '%s_quan_energy.%s'%(car.name,self.plot_format)
        plt.savefig(outname)
        print(outname)
