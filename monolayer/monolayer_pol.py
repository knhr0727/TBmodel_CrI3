import numpy as np
import matplotlib.pyplot as plt

Bohr = 0.52917721092 #Angstrom
Hartree = 27.2113860217 #eV

def Rsocp(theta):
    R = np.array([[ np.cos(theta),-np.sin(theta), 0.],\
                  [ np.sin(theta), np.cos(theta), 0.],\
                  [ 0.           , 0.           , 1.]])
    return R

def Rsocd(theta):
    R1 = np.array([[ np.cos(2*theta),-np.sin(2*theta)],\
                   [ np.sin(2*theta), np.cos(2*theta)]])
    R2 = np.array([[ np.cos(theta),-np.sin(theta)],\
                   [ np.sin(theta), np.cos(theta)]])
    R = np.zeros((5,5),dtype=float)
    R[0,0] = 1.
    R[1:3,1:3] = R1
    R[3:5,3:5] = R2
    return R

def kpath(k1,k2,n):
    path = []
    for x in np.linspace(0.,1.,n):
        path.append(k1*(1.-x)+k2*x)
    return path

def PlotBand(Band_list, kticks_label=None, yrange=None, shift=False,
             eV=False, EF=None, highlight=None, save=False,
             fname=None, c1='b', c2='r', figsize=None):
    #if (save):
    #    #del(sys.modules['matplotlib'])
    #    import matplotlib
    #    matplotlib.use('Agg')
    #import matplotlib.pyplot as plt
    try:
        len(Band_list[0][0])
    except:
        Band_list = [Band_list]
    fig = plt.figure(figsize=figsize)
    ax = plt.subplot()
    kticks = []
    kstart = 0.
    i1,i2 = 1,-1
    if (highlight != None):
        i1,i2 = highlight[0],highlight[1]
    e0 = 0.0
    if ((EF!=None)and(shift==True)):
        e0 = EF
    for Bands in Band_list:
        klist = Bands[0]+kstart
        kticks.append(klist[0])
        for i in range(len(Bands)-1):
            E = Bands[i+1].copy() - e0
            if (eV == True): E *= Hartree
            col = c1
            if ((i1<=i)and(i<=i2)): col = c2
            plt.plot(klist,E,color=col)
        kstart = klist[-1]
    kticks.append(kstart)
    ax.set_xticks(kticks, minor=False)
    ax.xaxis.grid(True, which='major')
    if (kticks_label != None):
        ax.set_xticklabels(kticks_label)
    if (EF != None):
        if (eV == True):
            EF *= Hartree
            e0 *= Hartree
        plt.plot([0.,kstart],[EF-e0,EF-e0],lw=0.25,color='gray',ls='--')
    plt.xlim(0,kstart)
    if (yrange != None):
        plt.ylim([yrange[0],yrange[1]])
    if(save):
        if(fname != None):
            plt.savefig(fname)
        else:
            plt.savefig('./pymx_band.png')
    else: plt.show()


#Atoms.UnitVectors.Unit             Ang
#<Atoms.UnitVectors
a = np.array([    6.9852500000,       0.0000000000,       0.0000000000])/Bohr
b = np.array([    3.4926250000,       6.0494039518,       0.0000000000])/Bohr
c = np.array([    0.0000000000,       0.0000000000,      24.0000000000])/Bohr
#Atoms.UnitVectors>
#Atoms.Number                       8
#Atoms.SpeciesAndCoordinates.Unit   frac
#<Atoms.SpeciesAndCoordinates
Cr1 =  0.6666666667*a +   0.6666666667*b +   0.500000000*c    
Cr2 =  0.3333333333*a +   0.3333333333*b +   0.500000000*c   
I3  =  0.3572400000*a +   0.6427600000*b +   0.565970000*c   
I4  =  0.0000000000*a +   0.3572400000*b +   0.565970000*c   
I5  =  0.6427600000*a +   0.0000000000*b +   0.565970000*c   
I6  =  0.0000000000*a +   0.6427600000*b +   0.434030000*c   
I7  =  0.3572400000*a +   0.0000000000*b +   0.434030000*c   
I8  =  0.6427600000*a +   0.3572400000*b +   0.434030000*c   
#Atoms.SpeciesAndCoordinates>
#
V = np.dot(a,np.cross(b,c))
ga = 2.*np.pi*np.cross(b,c)/V
gb = 2.*np.pi*np.cross(c,a)/V
gc = 2.*np.pi*np.cross(a,b)/V

theta3 = np.pi*7./6.
theta4 = -np.pi/6.
theta5 = np.pi/2.
theta8 = np.pi*7./6.
theta6 = -np.pi/6.
theta7 = np.pi/2.

Efermi = -3.0455/Hartree 

lambda_d = 0.0507/Hartree 
lambda_p = 0.6242/Hartree 
socd = lambda_d/2.
socp = lambda_p/2.

socp00 = np.zeros((3,3),dtype=complex)
socp00[0,1],socp00[1,0] = -1.j,1.j
socp01 = np.zeros((3,3),dtype=complex)
socp01[0,2],socp01[1,2],socp01[2,0],socp01[2,1] = 1.,-1.j,-1.,1.j
socp11 = -1.*socp00.copy()
socp10 = socp01.transpose().conjugate()
r3 = np.sqrt(3.)
socd00 = np.zeros((5,5),dtype=complex)
socd00[1,2],socd00[2,1],socd00[3,4],socd00[4,3] = -2.j,2.j,-1.j,1.j
socd01 = np.zeros((5,5),dtype=complex)
socd01[0,3],socd01[0,4],socd01[1,3],socd01[1,4],socd01[2,3],socd01[2,4] = \
        -r3,1.j*r3,1.,1.j,-1.j,1.
socd01[3,0],socd01[4,0],socd01[3,1],socd01[4,1],socd01[3,2],socd01[4,2] = \
        r3,-1.j*r3,-1.,-1.j,1.j,-1.
socd11 = -1.*socd00.copy()
socd10 = socd01.transpose().conjugate()


#TB

atomnum = [1,2,3,4,5,6,7,8]
atomr = {1:Cr1, 2:Cr2, 3:I3, 4:I4, 5:I5, 6:I6, 7:I7, 8:I8}
atomlen = [0,5,5,3,3,3,3,3,3] #first is dummy
Nh = 0
atomind = [0] #first is dummy
for i in atomlen:
    Nh += i
    atomind.append(Nh)
cellind = [[0,0],[-1,0],[1,0],[0,-1],[0,1],[-1,-1],[-1,1],[1,-1],[1,1]]
#           0     1      2     3      4     5       6      7      8
Ncell = 9
cell = []
for ab in cellind:
    R = ab[0]*a + ab[1]*b
    cell.append(R)
cell = np.array(cell).transpose()
#print(cell)
def cellr(i):
    return cell[:,i].flatten()

def pos(ac):
    ra = atomr[ac[0]]
    rc = cellr(ac[1])
    return ra+rc

def phase_tau(k):
    ekt = np.zeros((Nh,Nh),dtype=complex)
    ii = 0
    for i in atomnum:
        t = atomr[i]
        l = atomlen[i]
        kt = np.exp(1.j*np.dot(k,t))
        for j in range(l):
            ekt[ii,ii] = kt
            ii += 1
    return np.kron(np.identity(2),ekt)

f = np.load('../HR.R.0.npz')
HR0 = f['HR']
f.close()
f = np.load('../HR.R.1.npz')
HR1 = f['HR']
f.close()

thetas = {3:theta3,4:theta4,5:theta5,6:theta6,7:theta7,8:theta8}

HR = np.zeros((Ncell,2*Nh,2*Nh),dtype=complex)
for i in range(Ncell):
    HR[i,:Nh,:Nh] = HR0[i,:,:]
    HR[i,Nh:,Nh:] = HR1[i,:,:]
for i in atomnum:
    iind = atomind[i]
    iend = iind+atomlen[i]
    if (i==1) or (i==2):
        HR[0,iind:iend,iind:iend] += socd*socd00
        HR[0,iind:iend,iind+Nh:iend+Nh] += socd*socd01
        HR[0,iind+Nh:iend+Nh,iind:iend] += socd*socd10
        HR[0,iind+Nh:iend+Nh,iind+Nh:iend+Nh] += socd*socd11
    else:
        t = thetas[i]
        Rp = Rsocp(t)
        RpT = Rsocp(-t)
        HR[0,iind:iend,iind:iend] += socp*np.linalg.multi_dot([RpT,socp00,Rp])
        HR[0,iind:iend,iind+Nh:iend+Nh] += socp*np.linalg.multi_dot([RpT,socp01,Rp])
        HR[0,iind+Nh:iend+Nh,iind:iend] += socp*np.linalg.multi_dot([RpT,socp10,Rp])
        HR[0,iind+Nh:iend+Nh,iind+Nh:iend+Nh] += socp*np.linalg.multi_dot([RpT,socp11,Rp])
        

def Hk(k):
    exp_vec = np.exp(1.j*np.dot(k,cell))
    ekt = phase_tau(k)
    H0 = np.tensordot(exp_vec,HR,axes=((0),(0)))
    return np.linalg.multi_dot([ekt.transpose().conjugate(),H0,ekt])


def Band(k1, k2, n):
    path = kpath(k1,k2,n)
    k = 0.0
    klist = []
    Elists = []
    kbefore = path[0]
    for kvec in path:
        k += np.linalg.norm(kvec-kbefore)
        klist.append(k)
        h = Hk(kvec)
        w = np.linalg.eigvalsh(h)
        Elists.append(w)
        kbefore = kvec
    E = np.transpose(np.array(Elists))
    out = [np.array(klist)]
    for i in range(E.shape[0]):
        out.append(E[i,:])
    #out[0]:k-axis value from 0, out[1:]:energy eigenvalues
    return out


Nocc = 42
Nz = 8
Nx = 9
frac = np.linspace(0.,1.,Nx+1)
fracz = np.linspace(0.,1.,Nz+1)
pols = []
for x in range(Nx):
    for y in range(Nx):
        det = 1.+0.j
        kx = frac[x]*ga + frac[y]*gb
        h = Hk(kx)
        w,v = np.linalg.eigh(h)
        v0 = v[:,:Nocc]
        vi = v0.copy()
        for z in range(1,Nz):
            k = kx + fracz[z]*gc
            h = Hk(k)
            w,v = np.linalg.eigh(h)
            v = v[:,:Nocc]
            Smat = np.dot(v0.conjugate().transpose(),v)
            det *= np.linalg.det(Smat)
            v0 = v.copy()
        v = np.dot(phase_tau(-gc),vi)
        Smat = np.dot(v0.conjugate().transpose(),v)
        det *= np.linalg.det(Smat)
        pols.append(-1.*np.log(det).imag)

#print(pols)
pol_sum = sum(pols)
#print(pol_sum)        

AU2Debye =  2.54174776
AU2Mucm  =  5721.4765758 
Debye2Mucm = 2251.00209

#d = 19.807*2./3./Bohr  #bilayer
d = 19.807*1./3./Bohr  #monolayer
Area = np.linalg.norm(np.cross(ga,gb))
rArea = np.linalg.norm(np.cross(a,b))
cellV = rArea*d
#print(cellV,Area)

zak_integral = Area*pol_sum/float(Nx*Nx)
Pe = -1./np.power(2*np.pi,3)*zak_integral  # e/Bohr^2 (AU)
Pe_Debye = AU2Debye*Pe*cellV
#P_tot_Debye = Pe_Debye + Cdp_proj + Bdp_proj

#core part
z_ref =  c[2]/2.
Cr_core = 6.
I_core = 5.
Cr1core = Cr_core*(Cr1[2]-z_ref)
Cr2core = Cr_core*(Cr2[2]-z_ref)
I3core = I_core*(I3[2]-z_ref)
I4core = I_core*(I4[2]-z_ref)
I5core = I_core*(I5[2]-z_ref)
I6core = I_core*(I6[2]-z_ref)
I7core = I_core*(I7[2]-z_ref)
I8core = I_core*(I8[2]-z_ref)
Pc_eBohr = Cr1core+Cr2core+I3core+I4core+I5core+I6core+I7core+I8core
Pc_Debye = Pc_eBohr*Bohr

# muC/cm^2
Pe_Mucm = Pe_Debye*AU2Mucm/AU2Debye/cellV
Pc_Mucm = Pc_Debye*AU2Mucm/AU2Debye/cellV
Ptot_Mucm = Pe_Mucm+Pc_Mucm
print("electron contribution: ",Pe_Mucm)
print("core contribution: ",Pc_Mucm)
print("total: ",Ptot_Mucm)


