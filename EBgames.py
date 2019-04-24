from go import *
if '/home/dcollins/repos/p49c/p49_eigenmodes' not in sys.path:
    sys.path.append('/home/dcollins/repos/p49c/p49_eigenmodes')

import p49_eigen
reload(p49_eigen)
import p49_plot_tools
reload(p49_plot_tools)
import p49_QU2EB
import turb_quan
reload(turb_quan)
from p49_print_tools import *

def EB_predictor(n0,H0,Hx,Hy,dhx,dhy,dhz,dn,theta,psi,gamma=-2):
    #Psi is from the positive x.
    #H0 is in the x-z plane.
    #z is the line of sight.
    #ep = eq + i eu
    #ep = AnH^gamm( Hx+iHy)^2
    #eq = AnH^gamma ( Hx^2-Hy^2)
    #eu = AnH^gamma 2 Hx Hy
    #E+iB = (Q+iU)exp(-2 i Psi)
    #
    #MODE WISE
    A=1
    cos2psi = np.cos(2*psi)
    sin2psi = np.sin(2*psi)
    sintheta= np.sin(theta)
    costheta= np.cos(theta)
    Q0 = A*n0*H0**gamma*(Hx**2-Hy**2)
    U0 = A*n0*H0**gamma*(2*Hx*Hy)
    E0 = Q0*cos2psi+U0*sin2psi
    B0 = U0*cos2psi-Q0*sin2psi

    dH = sintheta*dhx+costheta*dhz
    dQ = A*n0*H0**(2+gamma)*(2*sintheta+dhx/H0+
                             gamma*sintheta**2*dH/H0+
                             sintheta**2*dn/n0)
    dU = A*n0*H0**(2+gamma)*(2*sintheta*dhy/H0)
    dE = dQ*cos2psi+dU*sin2psi
    dB = dU*cos2psi-dQ*sin2psi

    output={}
    output['Q0']=Q0
    output['U0']=U0
    output['E0']=E0
    output['B0']=B0
    output['dQ']=dQ
    output['dU']=dU
    output['dE']=dE
    output['dB']=dB
    return output

def norm(vec):
    return np.sqrt( (vec*vec).sum() )
def state_to_caldwell(stuff, line_of_sight):
    #xyz and 012 refer to these components.
    #XX YY ZZ refer to the Caldwell coordinates,
    #   theta is the angle of H0 relative to the line of sight.
    #   Psi is the angle of k relative to the plane-of-sky field.
    #   where ZZ is the line of sight
    #   XX is the plane-of-sky field
    #   YY is the other one.
    #   H_caldwell = H0 ( sin theta, 0, costheta)
    #   K_caldwell = k (cos Psi, sin Psi, 0)
    b0 = stuff['s'].b0+0
    Hmag = np.sqrt((b0*b0).sum())
    b_plane = b0+0
    b_plane[line_of_sight]=0

    mask = stuff['s'].kstuff['mask']
    n_k = mask.sum()
    kx = stuff['s'].kstuff['k'][0,...][mask]
    ky = stuff['s'].kstuff['k'][1,...][mask]
    kz = stuff['s'].kstuff['k'][2,...][mask]
    k_norm = (kx**2+ky**2+kz**2)**0.5

    cospsi = np.zeros(n_k)
    for n,kxyz in enumerate(zip(kx,ky,kz)):
        cospsi[n] = (kxyz*b_plane).sum()/(norm(nar(kxyz))*norm(b_plane))
        print("cosPSI",cospsi, np.arccos(cospsi))
    
def another_dumb_thing(psi,k0=1,h0=1,theta=0,wave='f-'):
    r2=np.sqrt(2)
    #state = p49_eigen.waves(hx=h0*np.cos(theta),hy=h0*np.sin(theta)/r2,hz=h0*np.sin(theta)/r2,p=0.6,this_wave=wave, HydroMethod=4)
    state = p49_eigen.waves(hx=h0*np.sin(theta),hy=0,hz=h0*np.cos(theta),p=0.6,this_wave=wave, HydroMethod=4)


    
    k = np.zeros([3,psi.size])
    k[0,:]=k0*np.cos(psi)
    k[1,:]=k0*np.sin(psi)
    k[2,:]=0
    state.wave_to_fields(k,wave)



    eps=1e-3
    EB=EB_predictor(state.quan['d'],
                 norm(state.b0),
                 state.quan['hx'],
                 state.quan['hy'],
                 eps*state.rot['hx'],
                 eps*state.rot['hy'],
                 eps*state.rot['hz'],
                 eps*state.rot['d'],
                 theta,psi,gamma=-2)


    plt.plot(psi,EB['E0'],label='E0',c='r')
    plt.plot(psi,EB['B0'],label='B0',c='g')
    plt.plot(psi,EB['dE'],label='dE',c='r',linestyle="--")
    plt.plot(psi,EB['dB'],label='dB',c='g',linestyle="--")
    plt.legend(loc=0)
    #plt.ylim([2*EB['B0'].min(),2*EB['B0'].max()])
    output={}
    output['s']=state
    return state
psi = np.arange(0,2*np.pi,0.1)
plt.clf()
for n in [0.2]: #[0.1,0.2,0.25]:
    state=another_dumb_thing(psi,theta=n*np.pi)
plt.savefig('p49c_plots/psi_test_E0.png')


#state_to_caldwell(stuff,0)
