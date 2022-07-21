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

# interlayer 
ilsg, ilpi =  0.4112/Hartree, -0.0603/Hartree
il_ref =  4.1284/Bohr #Bohr
il_cut = 13.200209411618099 #Bohr
il_scale_sg = 0.6464/Bohr #Bohr
il_scale_pi = 0.4446/Bohr #Bohr

def ppab(l,m,n):
    a = np.zeros((3,3),dtype=float)
    a[0,0],a[0,1],a[0,2] = l*l, l*m, l*n
    a[1,0],a[1,1],a[1,2] = m*l, m*m, m*n
    a[2,0],a[2,1],a[2,2] = n*l, n*m, n*n
    b = np.zeros((3,3),dtype=float)
    b[0,0],b[0,1],b[0,2] = 1.-l*l, -l*m, -l*n
    b[1,0],b[1,1],b[1,2] = -m*l, 1.-m*m, -m*n
    b[2,0],b[2,1],b[2,2] = -n*l, -n*m, 1.-n*n
    return [a,b]

def hpp(theta0,theta1,r):
    r_norm = np.linalg.norm(r)
    if (r_norm > il_cut):
        return np.zeros((3,3),dtype=float)
    else:
        x0 = np.array([np.cos(theta0),np.sin(theta0),0.])
        y0 = np.array([-1.*np.sin(theta0),np.cos(theta0),0.])
        x1 = np.array([np.cos(theta1),np.sin(theta1),0.])
        y1 = np.array([-1.*np.sin(theta1),np.cos(theta1),0.])
        z0  = np.array([0.,0.,1.])
        z1  = np.array([0.,0.,1.])
        r0 = r/np.linalg.norm(r)
        l,m,n = np.dot(r0,x0), np.dot(r0,y0), np.dot(r0,z0)
        a,b = ppab(l,m,n)
        hpp_sg = a*ilsg*np.exp(-(r_norm-il_ref)/il_scale_sg) 
        hpp_pi = b*ilpi*np.exp(-(r_norm-il_ref)/il_scale_pi)
        hpp = hpp_sg + hpp_pi                                
        R = np.zeros((3,3),dtype=float)
        R[0,0],R[0,1] = np.dot(x0,x1),np.dot(x0,y1)
        R[1,0],R[1,1] = np.dot(y0,x1),np.dot(y0,y1)
        R[2,2] = 1.
        return np.dot(hpp,R)


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
#Cr1 =  0.6666666667*a +   0.6666666667*b +   0.500000000*c   
#Cr2 =  0.3333333333*a +   0.3333333333*b +   0.500000000*c   
#I3  =  0.3572400000*a +   0.6427600000*b +   0.565970000*c   
#I4  =  0.0000000000*a +   0.3572400000*b +   0.565970000*c   
#I5  =  0.6427600000*a +   0.0000000000*b +   0.565970000*c   
#I6  =  0.0000000000*a +   0.6427600000*b +   0.434030000*c   
#I7  =  0.3572400000*a +   0.0000000000*b +   0.434030000*c   
#I8  =  0.6427600000*a +   0.3572400000*b +   0.434030000*c   
#Atoms.SpeciesAndCoordinates>
#
d = np.array([0.,0.,1.])*(0.651978995000 - 0.348021005000)*24.0/Bohr
Cr1b =  0.6666666667*a +   0.6666666667*b + 0.500000000*c -d/2.  
Cr2b =  0.3333333333*a +   0.3333333333*b + 0.500000000*c -d/2.  
I3b  =  0.3572400000*a +   0.6427600000*b + 0.565970000*c -d/2.  
I4b  =  0.0000000000*a +   0.3572400000*b + 0.565970000*c -d/2.  
I5b  =  0.6427600000*a +   0.0000000000*b + 0.565970000*c -d/2.  
I6b  =  0.0000000000*a +   0.6427600000*b + 0.434030000*c -d/2.  
I7b  =  0.3572400000*a +   0.0000000000*b + 0.434030000*c -d/2.  
I8b  =  0.6427600000*a +   0.3572400000*b + 0.434030000*c -d/2.  
Cr1t =  0.6666666667*a +   0.6666666667*b + 0.500000000*c +d/2. -(1./3.)*(a+b) 
Cr2t =  0.3333333333*a +   0.3333333333*b + 0.500000000*c +d/2. -(1./3.)*(a+b) 
I3t  =  0.3572400000*a +   0.6427600000*b + 0.565970000*c +d/2. -(1./3.)*(a+b) 
I4t  =  0.0000000000*a +   0.3572400000*b + 0.565970000*c +d/2. -(1./3.)*(a+b) 
I5t  =  0.6427600000*a +   0.0000000000*b + 0.565970000*c +d/2. -(1./3.)*(a+b) 
I6t  =  0.0000000000*a +   0.6427600000*b + 0.434030000*c +d/2. -(1./3.)*(a+b) 
I7t  =  0.3572400000*a +   0.0000000000*b + 0.434030000*c +d/2. -(1./3.)*(a+b) 
I8t  =  0.6427600000*a +   0.3572400000*b + 0.434030000*c +d/2. -(1./3.)*(a+b) 

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
#atomr = {1:Cr1, 2:Cr2, 3:I3, 4:I4, 5:I5, 6:I6, 7:I7, 8:I8}
atomrt = {1:Cr1t, 2:Cr2t, 3:I3t, 4:I4t, 5:I5t, 6:I6t, 7:I7t, 8:I8t}
atomrb = {1:Cr1b, 2:Cr2b, 3:I3b, 4:I4b, 5:I5b, 6:I6b, 7:I7b, 8:I8b}
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

#def pos(ac):
#    ra = atomr[ac[0]]
#    rc = cellr(ac[1])
#    return ra+rc
def pos(ac,tb):
    if (tb=='t'):
        ra = atomrt[ac[0]]
    elif (tb=='b'):
        ra = atomrb[ac[0]]
    else:
        raise Exception
    rc = cellr(ac[1])
    return ra+rc

def phase_tau(k):
    ekt = np.zeros((Nh,Nh),dtype=complex)
    ii = 0
    for i in atomnum:
        t = atomrt[i]
        l = atomlen[i]
        kt = np.exp(1.j*np.dot(k,t))
        for j in range(l):
            ekt[ii,ii] = kt
            ii += 1
    top = np.kron(np.identity(2),ekt)
    ekt = np.zeros((Nh,Nh),dtype=complex)
    ii = 0
    for i in atomnum:
        t = atomrb[i]
        l = atomlen[i]
        kt = np.exp(1.j*np.dot(k,t))
        for j in range(l):
            ekt[ii,ii] = kt
            ii += 1
    bot = np.kron(np.identity(2),ekt)
    zero = np.zeros((2*Nh,2*Nh),dtype=complex)
    return np.block([[bot,zero],[zero,top]])

#method cellperiodic
taux11 = np.zeros((Nh,Nh),dtype=float)
tauy11 = np.zeros((Nh,Nh),dtype=float)
tauz11 = np.zeros((Nh,Nh),dtype=float)
taux22 = np.zeros((Nh,Nh),dtype=float)
tauy22 = np.zeros((Nh,Nh),dtype=float)
tauz22 = np.zeros((Nh,Nh),dtype=float)
taux12 = np.zeros((Nh,Nh),dtype=float)
tauy12 = np.zeros((Nh,Nh),dtype=float)
tauz12 = np.zeros((Nh,Nh),dtype=float)
taux21 = np.zeros((Nh,Nh),dtype=float)
tauy21 = np.zeros((Nh,Nh),dtype=float)
tauz21 = np.zeros((Nh,Nh),dtype=float)
ii = 0
for i in atomnum:
    tb = atomrb[i]
    tt = atomrt[i]
    l = atomlen[i]
    for j in range(l):
        taux11[ii,:] -= tb[0]
        taux11[:,ii] += tb[0]
        tauy11[ii,:] -= tb[1]
        tauy11[:,ii] += tb[1]
        tauz11[ii,:] -= tb[2]
        tauz11[:,ii] += tb[2]
        taux22[ii,:] -= tt[0]
        taux22[:,ii] += tt[0]
        tauy22[ii,:] -= tt[1]
        tauy22[:,ii] += tt[1]
        tauz22[ii,:] -= tt[2]
        tauz22[:,ii] += tt[2]
        taux12[ii,:] -= tb[0]
        taux12[:,ii] += tt[0]
        tauy12[ii,:] -= tb[1]
        tauy12[:,ii] += tt[1]
        tauz12[ii,:] -= tb[2]
        tauz12[:,ii] += tt[2]
        taux21[ii,:] -= tt[0]
        taux21[:,ii] += tb[0]
        tauy21[ii,:] -= tt[1]
        tauy21[:,ii] += tb[1]
        tauz21[ii,:] -= tt[2]
        tauz21[:,ii] += tb[2]
        ii += 1
taux11 = np.kron(np.ones((2,2)),taux11)
tauy11 = np.kron(np.ones((2,2)),tauy11)
tauz11 = np.kron(np.ones((2,2)),tauz11)
taux22 = np.kron(np.ones((2,2)),taux22)
tauy22 = np.kron(np.ones((2,2)),tauy22)
tauz22 = np.kron(np.ones((2,2)),tauz22)
taux12 = np.kron(np.ones((2,2)),taux12)
tauy12 = np.kron(np.ones((2,2)),tauy12)
tauz12 = np.kron(np.ones((2,2)),tauz12)
taux21 = np.kron(np.ones((2,2)),taux21)
tauy21 = np.kron(np.ones((2,2)),tauy21)
tauz21 = np.kron(np.ones((2,2)),tauz21)
taux = np.block([[taux11,taux12],[taux21,taux22]])
tauy = np.block([[tauy11,tauy12],[tauy21,tauy22]])
tauz = np.block([[tauz11,tauz12],[tauz21,tauz22]])


f = np.load('../HR.R.0.npz')
HR0 = f['HR']
f.close()
f = np.load('../HR.R.1.npz')
HR1 = f['HR']
f.close()

il_theta = {3:theta3,4:theta4,5:theta5,6:theta6,7:theta7,8:theta8}

HRbot = np.zeros((Ncell,2*Nh,2*Nh),dtype=complex)
for i in range(Ncell):
    HRbot[i,:Nh,:Nh] = HR0[i,:,:]
    HRbot[i,Nh:,Nh:] = HR1[i,:,:]
for i in atomnum:
    iind = atomind[i]
    iend = iind+atomlen[i]
    if (i==1) or (i==2):
        t = 0.
        Rd = Rsocd(t)
        RdT = Rsocd(-t)
        HRbot[0,iind:iend,iind:iend] += socd*np.linalg.multi_dot([RdT,socd00,Rd])
        HRbot[0,iind:iend,iind+Nh:iend+Nh] += socd*np.linalg.multi_dot([RdT,socd01,Rd])
        HRbot[0,iind+Nh:iend+Nh,iind:iend] += socd*np.linalg.multi_dot([RdT,socd10,Rd])
        HRbot[0,iind+Nh:iend+Nh,iind+Nh:iend+Nh] += socd*np.linalg.multi_dot([RdT,socd11,Rd])
    else:
        t = il_theta[i]
        Rp = Rsocp(t)
        RpT = Rsocp(-t)
        HRbot[0,iind:iend,iind:iend] += socp*np.linalg.multi_dot([RpT,socp00,Rp])
        HRbot[0,iind:iend,iind+Nh:iend+Nh] += socp*np.linalg.multi_dot([RpT,socp01,Rp])
        HRbot[0,iind+Nh:iend+Nh,iind:iend] += socp*np.linalg.multi_dot([RpT,socp10,Rp])
        HRbot[0,iind+Nh:iend+Nh,iind+Nh:iend+Nh] += socp*np.linalg.multi_dot([RpT,socp11,Rp])
        

HRtop = np.zeros((Ncell,2*Nh,2*Nh),dtype=complex)
for i in range(Ncell):
    HRtop[i,:Nh,:Nh] = HR1[i,:,:]
    HRtop[i,Nh:,Nh:] = HR0[i,:,:]
for i in atomnum:
    iind = atomind[i]
    iend = iind+atomlen[i]
    if (i==1) or (i==2):
        t = 0.
        Rd = Rsocd(t)
        RdT = Rsocd(-t)
        HRtop[0,iind:iend,iind:iend] += socd*np.linalg.multi_dot([RdT,socd00,Rd])
        HRtop[0,iind:iend,iind+Nh:iend+Nh] += socd*np.linalg.multi_dot([RdT,socd01,Rd])
        HRtop[0,iind+Nh:iend+Nh,iind:iend] += socd*np.linalg.multi_dot([RdT,socd10,Rd])
        HRtop[0,iind+Nh:iend+Nh,iind+Nh:iend+Nh] += socd*np.linalg.multi_dot([RdT,socd11,Rd])
    else:
        t = il_theta[i]
        Rp = Rsocp(t)
        RpT = Rsocp(-t)
        HRtop[0,iind:iend,iind:iend] += socp*np.linalg.multi_dot([RpT,socp00,Rp])
        HRtop[0,iind:iend,iind+Nh:iend+Nh] += socp*np.linalg.multi_dot([RpT,socp01,Rp])
        HRtop[0,iind+Nh:iend+Nh,iind:iend] += socp*np.linalg.multi_dot([RpT,socp10,Rp])
        HRtop[0,iind+Nh:iend+Nh,iind+Nh:iend+Nh] += socp*np.linalg.multi_dot([RpT,socp11,Rp])
        

HR = np.zeros((Ncell,4*Nh,4*Nh),dtype=complex)  
for i in range(Ncell):
    HR[i,:2*Nh,:2*Nh] = HRbot[i,:,:]
    HR[i,2*Nh:,2*Nh:] = HRtop[i,:,:]

#interlayer
botI = [3,4,5]
topI = [6,7,8]
il_theta = {3:theta3,4:theta4,5:theta5,6:theta6,7:theta7,8:theta8}
for i in botI:
    rb = atomrb[i]
    thei = il_theta[i]
    iind = atomind[i]
    iend = iind+atomlen[i]
    for j in topI:
        rt = atomrt[j]
        thej = il_theta[j]
        jind = atomind[j]
        jend = jind+atomlen[j]
        for iR in range(Ncell):
            rR = cell[:,iR].transpose()
            r1 = rt+rR-rb
            r2 = rb+rR-rt
            h1 = hpp(thei,thej,r1)
            h2 = hpp(thej,thei,r2)
            HR[iR,iind:iend,jind+2*Nh:jend+2*Nh] += h1
            HR[iR,iind+Nh:iend+Nh,jind+3*Nh:jend+3*Nh] += h1
            HR[iR,jind+2*Nh:jend+2*Nh,iind:iend] += h2
            HR[iR,jind+3*Nh:jend+3*Nh,iind+Nh:iend+Nh] += h2


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

######## Bands plot part ######
#G = 0.000*ga + 0.000*gb + 0.000*gc
#K = (2./3.)*ga + (1./3.)*gb + 0.000*gc
#M = 0.500*ga + 0.500*gb + 0.000*gc
#
#band1 = Band(G,K,40)
#band2 = Band(K,M,30)
#band3 = Band(M,G,40)
#PlotBand([band1,band2,band3],kticks_label=['gamma','K','M','gamma'],\
#        shift=True, eV=True, EF=Efermi)
###############################


vxlist = []
vylist = []
vzlist = []
elist = []
N = 9
frac = np.linspace(0.,1.,N+1)
for x in range(N):
    for y in range(N):
        k = frac[x]*ga + frac[y]*gb
        h = Hk(k)
        w,v = np.linalg.eigh(h)
        tauxh = taux*h
        tauyh = tauy*h
        tauzh = tauz*h
        exp_vec = np.exp(1.j*np.dot(k,cell))
        ekt = phase_tau(k)
        Rx_exp_vec = exp_vec*cell[0,:]
        Ry_exp_vec = exp_vec*cell[1,:]
        Rz_exp_vec = exp_vec*cell[2,:]
        RxH0 = np.tensordot(Rx_exp_vec,HR,axes=((0),(0)))
        RyH0 = np.tensordot(Ry_exp_vec,HR,axes=((0),(0)))
        RzH0 = np.tensordot(Rz_exp_vec,HR,axes=((0),(0)))
        RxH = np.linalg.multi_dot([ekt.transpose().conjugate(),RxH0,ekt])
        RyH = np.linalg.multi_dot([ekt.transpose().conjugate(),RyH0,ekt])
        RzH = np.linalg.multi_dot([ekt.transpose().conjugate(),RzH0,ekt])
        dxH = (RxH + tauxh)*1.j
        dyH = (RyH + tauyh)*1.j
        dzH = (RzH + tauzh)*1.j
        vx = np.linalg.multi_dot([v.conjugate().transpose(),dxH,v])
        vy = np.linalg.multi_dot([v.conjugate().transpose(),dyH,v])
        vz = np.linalg.multi_dot([v.conjugate().transpose(),dzH,v])
        vxlist.append(vx)
        vylist.append(vy)
        vzlist.append(vz)
        elist.append(w)


Emax = 10. # eV
Nomg = 1001
Nocc = 42*2
Nk = N*N #81
eta = 0.1 # eV
Area = np.linalg.norm(np.cross(a,b))
d_layer = 19.807*2./3./Bohr  #bilayer

Emax = Emax/Hartree
eta = eta/Hartree
E = Emax*np.linspace(0.,1.,Nomg)
Sxx = np.zeros(Nomg,dtype=complex)
Sxy = np.zeros(Nomg,dtype=complex)
Syx = np.zeros(Nomg,dtype=complex)

for vx,vy,e in zip(vxlist,vylist,elist):
    nmax = vx.shape[0]
    sxx = np.zeros(Nomg,dtype=complex)
    sxy = np.zeros(Nomg,dtype=complex)
    syx = np.zeros(Nomg,dtype=complex)
    for m in range(nmax):
        em = e[m]
        for n in range(nmax):
            if ((m<Nocc)and(n>=Nocc)):
                en = e[n]
                vxmn = vx[m,n]
                vxnm = vx[n,m]
                vymn = vy[m,n]
                vynm = vy[n,m]
                efactorxx = 1./(em-en)/(E+em-en+1.j*eta)
                efactor = 1./((em-en)*(em-en)-(E+1.j*eta)*(E+1.j*eta))
                sxx += vxmn*vxnm*efactorxx
                sxy += (vxmn*vynm).imag*efactor
                syx += (vymn*vxnm).imag*efactor
            if ((m>=Nocc)and(n<Nocc)):
                en = e[n]
                vxmn = vx[m,n]
                vxnm = vx[n,m]
                vymn = vy[m,n]
                vynm = vy[n,m]
                efactorxx = -1./(em-en)/(E+em-en+1.j*eta)
                efactor = -1./((em-en)*(em-en)-(E+1.j*eta)*(E+1.j*eta))
                sxx += vxmn*vxnm*efactorxx
                sxy += (vxmn*vynm).imag*efactor
                syx += (vymn*vxnm).imag*efactor
    Sxx += sxx
    Sxy += sxy
    Syx += syx

Sxx = Sxx*(-1.j)/Nk/Area/d_layer
Sxy = Sxy/Nk/Area/d_layer
Syx = Syx/Nk/Area/d_layer
Sxxout = Sxx.copy()
Sxyout = Sxy.copy()
Syxout = Syx.copy()
Eout = E.copy()

### preventing division by 0 in Kerr formula
Sxy[0] = 0.+0.j
Sxx[0] = 1.+0.j
Syx[0] = 0.+0.j
E[0] = 1.

Kerr = -Sxy/(Sxx*np.sqrt(4.j*np.pi*Sxx/E+1.))

np.savez('./DATA.biAFM.npz',E=Eout,Sxx=Sxxout,Sxy=Sxyout,Syx=Syxout,Kerr=Kerr)

