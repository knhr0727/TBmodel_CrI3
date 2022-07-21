import numpy as np
import matplotlib.pyplot as plt
import copy

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

############## TBL part ##############
def rotation(D,theta):
    R = np.array([[np.cos(theta),-1.*np.sin(theta),0],[np.sin(theta),np.cos(theta),0],[0.,0.,1.]])
    return np.dot(R,D)

# Define TBL by (n.m)
n = 2
m = 1

cos = float(m*m+n*n+4*m*n)/2./float(m*m+n*n+m*n)
theta = np.arccos(cos)
if (m>n):
    theta = -theta

a1 = np.array([1.,0.,0.])
a2 = np.array([0.5,0.5*np.sqrt(3.),0.])

L11 = n*a1+m*a2
L21 = -m*a1+(n+m)*a2

def x111(y):
    return L21[0]*(y/L21[1])
def x121(y):
    return L11[0]*(y/L11[1])
def x211(y):
    return L21[0] + L11[0]*((y-L21[1])/L11[1])
def x221(y):
    return L11[0] + L21[0]*((y-L11[1])/L21[1])

L12 = m*a1+n*a2
L22 = -n*a1+(m+n)*a2

def x112(y):
    return L22[0]*(y/L22[1])
def x122(y):
    return L12[0]*(y/L12[1])
def x212(y):
    return L22[0] + L12[0]*((y-L22[1])/L12[1])
def x222(y):
    return L12[0] + L22[0]*((y-L12[1])/L22[1])

cell_points1 = [] #top layer
pair1 = []
for j in range(n+2*m):
    for i in range(-m,n):
        r = i*a1+j*a2
        x = r[0]
        y = r[1]
        if ((j==0)and(i==0)):
            cell_points1.append(r)
            pair1.append([i,j])
        elif((j==n+m)and(i==-m)):
            pass
        elif(0<j<=m):
            if (x111(y)<=x<x121(y)):
                cell_points1.append(r)
                pair1.append([i,j])
        elif(m<j<=n+m):
            if (x111(y)<=x<x221(y)):
                cell_points1.append(r)
                pair1.append([i,j])
        elif(n+m<j):
            if (x211(y)<=x<x221(y)):
                cell_points1.append(r)
                pair1.append([i,j])

cell_points2 = [] #bottom layer
pair2 = []
for j in range(m+2*n):
    for i in range(-n,m):
        r = i*a1+j*a2
        x = r[0]
        y = r[1]
        if ((j==0)and(i==0)):
            cell_points2.append(r)
            pair2.append([i,j])
        elif((j==m+n)and(i==-n)):
            pass
        elif(0<j<=n):
            if (x112(y)<=x<x122(y)):
                cell_points2.append(r)
                pair2.append([i,j])
        elif(n<j<=m+n):
            if (x112(y)<=x<x222(y)):
                cell_points2.append(r)
                pair2.append([i,j])
        elif(m+n<j):
            if (x212(y)<=x<x222(y)):
                cell_points2.append(r)
                pair2.append([i,j])

#print(cell_points1)
#print(cell_points2)
#print(len(cell_points1))
#print(len(cell_points2))
#print(pair1)
#print(pair2)
inTBL = m*m+n*n+m*n
#print(inTBL)

neighbor1 = []
neighbor2 = []
cellind = [[0,0],[-1,0],[1,0],[0,-1],[0,1],[-1,-1],[-1,1],[1,-1],[1,1]]
#           0     1      2     3      4     5       6      7      8
#L11 = n*a1+m*a2
#L21 = -m*a1+(n+m)*a2
for P in pair1: #top
    i,j = P
    nlist = []
    for NP in cellind:
        ni,nj = NP
        I,J = i+ni,j+nj  #neighbor of [i,j]
        IJcell = None
        IJind = None
        for ii in range(-1,2): #L1
            for jj in range(-1,2): #L2
                II = I - n*ii + m*jj 
                JJ = J - m*ii - (n+m)*jj  
                if ([II,JJ] in pair1):
                    IJcell = cellind.index([ii,jj])
                    IJind = pair1.index([II,JJ])
        nlist.append([IJind,IJcell])
    neighbor1.append(copy.deepcopy(nlist))
#print(len(neighbor1))
#for l in neighbor1:
#    print(l)
#fig1 = plt.figure()
#ax1 = fig1.add_subplot(111,aspect='equal')
#colors = ['r','g','b','y','c','m','k']
#for nn,co in zip(neighbor1,colors):
#    i0,j0 = pair1[nn[0][0]]
#    x0,y0 = (i0*a1+j0*a2)[:2]
#    for ac in nn:
#        i,j = pair1[ac[0]]
#        I,J = cellind[ac[1]]
#        x,y = (i*a1+j*a2+I*L11+J*L21)[:2]
#        plt.plot([x0,x],[y0,y],c=co)
#plt.show()

#L12 = m*a1+n*a2
#L22 = -n*a1+(m+n)*a2
for P in pair2: #bot
    i,j = P
    nlist = []
    for NP in cellind:
        ni,nj = NP
        I,J = i+ni,j+nj  #neighbor of [i,j]
        IJcell = None
        IJind = None
        for ii in range(-1,2): #L1
            for jj in range(-1,2): #L2
                II = I - m*ii + n*jj 
                JJ = J - n*ii - (m+n)*jj  
                if ([II,JJ] in pair2):
                    IJcell = cellind.index([ii,jj])
                    IJind = pair2.index([II,JJ])
        nlist.append([IJind,IJcell])
    neighbor2.append(copy.deepcopy(nlist))
#print(len(neighbor2))
#for l in neighbor2:
#    print(l)
#fig1 = plt.figure()
#ax1 = fig1.add_subplot(111,aspect='equal')
#colors = ['r','g','b','y','c','m','k']
#for nn,co in zip(neighbor2,colors):
#    i0,j0 = pair2[nn[0][0]]
#    x0,y0 = (i0*a1+j0*a2)[:2]
#    for ac in nn:
#        i,j = pair2[ac[0]]
#        I,J = cellind[ac[1]]
#        x,y = (i*a1+j*a2+I*L12+J*L22)[:2]
#        plt.plot([x0,x],[y0,y],c=co)
#plt.show()




######################################

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
        hpp_sg = a*ilsg*np.exp(-(r_norm-il_ref)/il_scale_sg) #modified from here
        hpp_pi = b*ilpi*np.exp(-(r_norm-il_ref)/il_scale_pi)
        hpp = hpp_sg + hpp_pi                                #to here
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
A = 6.9852500000/Bohr
#Atoms.UnitVectors>
#Atoms.Number                       8
#Atoms.SpeciesAndCoordinates.Unit   frac
#<Atoms.SpeciesAndCoordinates
d = np.array([0.,0.,1.])*(0.651978995000 - 0.348021005000)*24.0/Bohr
Cr1b =  0.6666666667*a +   0.6666666667*b +   0.500000000*c -d/2. -(1./3.)*(a+b)                   
Cr2b =  0.3333333333*a +   0.3333333333*b +   0.500000000*c -d/2. -(1./3.)*(a+b)                   
I3b  =  0.3572400000*a +   0.6427600000*b +   0.434030000*c -d/2. -(1./3.)*(a+b)                   
I4b  =  0.0000000000*a +   0.3572400000*b +   0.434030000*c -d/2. -(1./3.)*(a+b)                   
I5b  =  0.6427600000*a +   0.0000000000*b +   0.434030000*c -d/2. -(1./3.)*(a+b)                   
I6b  =  0.0000000000*a +   0.6427600000*b +   0.565970000*c -d/2. -(1./3.)*(a+b)                   
I7b  =  0.3572400000*a +   0.0000000000*b +   0.565970000*c -d/2. -(1./3.)*(a+b)                   
I8b  =  0.6427600000*a +   0.3572400000*b +   0.565970000*c -d/2. -(1./3.)*(a+b)                   
Cr1t =  0.6666666667*a +   0.6666666667*b +   0.500000000*c +d/2. -(1./3.)*(a+b)    
Cr2t =  0.3333333333*a +   0.3333333333*b +   0.500000000*c +d/2. -(1./3.)*(a+b)    
I3t  =  0.3572400000*a +   0.6427600000*b +   0.434030000*c +d/2. -(1./3.)*(a+b)    
I4t  =  0.0000000000*a +   0.3572400000*b +   0.434030000*c +d/2. -(1./3.)*(a+b)    
I5t  =  0.6427600000*a +   0.0000000000*b +   0.434030000*c +d/2. -(1./3.)*(a+b)    
I6t  =  0.0000000000*a +   0.6427600000*b +   0.565970000*c +d/2. -(1./3.)*(a+b)    
I7t  =  0.3572400000*a +   0.0000000000*b +   0.565970000*c +d/2. -(1./3.)*(a+b)    
I8t  =  0.6427600000*a +   0.3572400000*b +   0.565970000*c +d/2. -(1./3.)*(a+b)    

############## TBL part ##############

L = A*np.sqrt(float(n*n+m*m+m*n))
L1 = L*np.array([1.,0.,0.])
L2 = L*np.array([0.5,0.5*np.sqrt(3.),0.])
L3 = 24./Bohr*np.array([0.,0.,1.])

L1t = n*a+m*b 
thetacell = np.arccos(np.dot(L1,L1t)/(L*np.linalg.norm(L1t)))
a1t = rotation(a ,-thetacell)
a2t = rotation(b ,-thetacell)
a1b = rotation(a ,-thetacell-theta)
a2b = rotation(b ,-thetacell-theta)
Cr1b  = rotation(Cr1b ,-thetacell-theta) 
Cr2b  = rotation(Cr2b ,-thetacell-theta) 
I3b   = rotation(I3b  ,-thetacell-theta) 
I4b   = rotation(I4b  ,-thetacell-theta) 
I5b   = rotation(I5b  ,-thetacell-theta) 
I6b   = rotation(I6b  ,-thetacell-theta) 
I7b   = rotation(I7b  ,-thetacell-theta) 
I8b   = rotation(I8b  ,-thetacell-theta) 
Cr1t  = rotation(Cr1t ,-thetacell) 
Cr2t  = rotation(Cr2t ,-thetacell) 
I3t   = rotation(I3t  ,-thetacell) 
I4t   = rotation(I4t  ,-thetacell) 
I5t   = rotation(I5t  ,-thetacell) 
I6t   = rotation(I6t  ,-thetacell) 
I7t   = rotation(I7t  ,-thetacell) 
I8t   = rotation(I8t  ,-thetacell) 

V = np.dot(L1,np.cross(L2,L3))
ga = 2.*np.pi*np.cross(L2,L3)/V
gb = 2.*np.pi*np.cross(L3,L1)/V
gc = 2.*np.pi*np.cross(L1,L2)/V

theta3b = np.pi/6.     -thetacell-theta 
theta4b = np.pi*5./6.  -thetacell-theta    
theta5b = -np.pi/2.    -thetacell-theta  
theta8b = np.pi/6.     -thetacell-theta 
theta6b = np.pi*5./6.  -thetacell-theta    
theta7b = -np.pi/2.    -thetacell-theta  
theta3t = np.pi/6.     -thetacell 
theta4t = np.pi*5./6.  -thetacell    
theta5t = -np.pi/2.    -thetacell  
theta8t = np.pi/6.     -thetacell 
theta6t = np.pi*5./6.  -thetacell    
theta7t = -np.pi/2.    -thetacell  

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
# TBL hamiltonian size:
Nhtbl = Nh*2*2*inTBL  # spin * top/bot * in each layer
dtop = int(Nhtbl/2)
#cellind = [[0,0],[-1,0],[1,0],[0,-1],[0,1],[-1,-1],[-1,1],[1,-1],[1,1]]
##           0     1      2     3      4     5       6      7      8
Ncell = 9
cell = []
for ab in cellind:
    R = ab[0]*L1 + ab[1]*L2   
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

###################### check positions ######################
#fig1 = plt.figure()
#ax1 = fig1.add_subplot(111,aspect='equal')
#xc0,yc0 = 0,0
#xc1,yc1 = L1[:2]
#xc2,yc2 = (L1+L2)[:2]
#xc3,yc3 = L2[:2]
#plt.plot([xc0,xc1,xc2,xc3,xc0],[yc0,yc1,yc2,yc3,yc0],c='k')
#for i,j in pair1:
#    R = i*a1t + j*a2t
#    for ia in atomnum:
#        r = R + atomrt[ia]
#        if (ia==1)or(ia==2):
#            co = 'b'
#        else:
#            co = 'r'
#        plt.scatter(r[0],r[1],c=co)
#for i,j in pair2:
#    R = i*a1b + j*a2b
#    for ia in atomnum:
#        r = R + atomrb[ia]
#        if (ia==1)or(ia==2):
#            co = 'b'
#        else:
#            co = 'r'
#        plt.scatter(r[0],r[1],c=co)
#plt.show()
#############################################################

taub = []
for i in atomnum:
    t = atomrb[i]
    l = atomlen[i]
    for j in range(l):
        taub.append(t)
taut = []
for i in atomnum:
    t = atomrt[i]
    l = atomlen[i]
    for j in range(l):
        taut.append(t)
tau_list = []
for i,j in pair2: #bottom
    R = i*a1b + j*a2b
    for t in taub: #for spin0 component
        tau_list.append(t+R)
    for t in taub: #for spin1 component
        tau_list.append(t+R)
for i,j in pair1: #top
    R = i*a1t + j*a2t
    for t in taut: #for spin0 component
        tau_list.append(t+R)
    for t in taut: #for spin1 component
        tau_list.append(t+R)
tau_array = np.array(tau_list).transpose()

def phase_tau(k):
    ekt = np.zeros((Nhtbl,Nhtbl),dtype=complex)
    exp_vec = np.exp(1.j*np.dot(k,tau_array))
    for i in range(Nhtbl):
        ekt[i,i] = exp_vec[i]
    return ekt

taux = np.zeros((Nhtbl,Nhtbl),dtype=float)
tauy = np.zeros((Nhtbl,Nhtbl),dtype=float)
tauz = np.zeros((Nhtbl,Nhtbl),dtype=float)
for i in range(Nhtbl):
    tx,ty,tz = tau_list[i]
    taux[i,:] -= tx
    taux[:,i] += tx
    tauy[i,:] -= ty
    tauy[:,i] += ty
    tauz[i,:] -= tz
    tauz[:,i] += tz


f = np.load('../HR.L.0.npz')
HR0b = f['HR']
f.close()
f = np.load('../HR.L.1.npz')
HR1b = f['HR']
f.close()

f = np.load('../HR.L.1.npz')
HR0t = f['HR']
f.close()
f = np.load('../HR.L.0.npz')
HR1t = f['HR']
f.close()

il_thetab = {3:theta3b,4:theta4b,5:theta5b,6:theta6b,7:theta7b,8:theta8b}
il_thetat = {3:theta3t,4:theta4t,5:theta5t,6:theta6t,7:theta7t,8:theta8t}

HRbot = np.zeros((Ncell,2*Nh,2*Nh),dtype=complex)
for i in range(Ncell):
    HRbot[i,:Nh,:Nh] = HR0b[i,:,:]
    HRbot[i,Nh:,Nh:] = HR1b[i,:,:]
for i in atomnum:
    iind = atomind[i]
    iend = iind+atomlen[i]
    if (i==1) or (i==2):
        t = -thetacell-theta
        Rd = Rsocd(t)
        RdT = Rsocd(-t)
        HRbot[0,iind:iend,iind:iend] += socd*np.linalg.multi_dot([RdT,socd00,Rd])
        HRbot[0,iind:iend,iind+Nh:iend+Nh] += socd*np.linalg.multi_dot([RdT,socd01,Rd])
        HRbot[0,iind+Nh:iend+Nh,iind:iend] += socd*np.linalg.multi_dot([RdT,socd10,Rd])
        HRbot[0,iind+Nh:iend+Nh,iind+Nh:iend+Nh] += socd*np.linalg.multi_dot([RdT,socd11,Rd])
    else:
        t = il_thetab[i]
        Rp = Rsocp(t)
        RpT = Rsocp(-t)
        HRbot[0,iind:iend,iind:iend] += socp*np.linalg.multi_dot([RpT,socp00,Rp])
        HRbot[0,iind:iend,iind+Nh:iend+Nh] += socp*np.linalg.multi_dot([RpT,socp01,Rp])
        HRbot[0,iind+Nh:iend+Nh,iind:iend] += socp*np.linalg.multi_dot([RpT,socp10,Rp])
        HRbot[0,iind+Nh:iend+Nh,iind+Nh:iend+Nh] += socp*np.linalg.multi_dot([RpT,socp11,Rp])
        

HRtop = np.zeros((Ncell,2*Nh,2*Nh),dtype=complex)
for i in range(Ncell):
    HRtop[i,:Nh,:Nh] = HR0t[i,:,:]
    HRtop[i,Nh:,Nh:] = HR1t[i,:,:]
for i in atomnum:
    iind = atomind[i]
    iend = iind+atomlen[i]
    if (i==1) or (i==2):
        t = -thetacell
        Rd = Rsocd(t)
        RdT = Rsocd(-t)
        HRtop[0,iind:iend,iind:iend] += socd*np.linalg.multi_dot([RdT,socd00,Rd])
        HRtop[0,iind:iend,iind+Nh:iend+Nh] += socd*np.linalg.multi_dot([RdT,socd01,Rd])
        HRtop[0,iind+Nh:iend+Nh,iind:iend] += socd*np.linalg.multi_dot([RdT,socd10,Rd])
        HRtop[0,iind+Nh:iend+Nh,iind+Nh:iend+Nh] += socd*np.linalg.multi_dot([RdT,socd11,Rd])
    else:
        t = il_thetat[i]
        Rp = Rsocp(t)
        RpT = Rsocp(-t)
        HRtop[0,iind:iend,iind:iend] += socp*np.linalg.multi_dot([RpT,socp00,Rp])
        HRtop[0,iind:iend,iind+Nh:iend+Nh] += socp*np.linalg.multi_dot([RpT,socp01,Rp])
        HRtop[0,iind+Nh:iend+Nh,iind:iend] += socp*np.linalg.multi_dot([RpT,socp10,Rp])
        HRtop[0,iind+Nh:iend+Nh,iind+Nh:iend+Nh] += socp*np.linalg.multi_dot([RpT,socp11,Rp])
        

HR = np.zeros((Ncell,Nhtbl,Nhtbl),dtype=complex) 
for i in range(inTBL):
    nn = neighbor2[i] # bot 
    for j in range(Ncell):
        ii,jj = nn[j] #cell ind small and large
        HR[jj,2*Nh*i:2*Nh*(i+1),2*Nh*ii:2*Nh*(ii+1)] += HRbot[j,:,:]
    nn = neighbor1[i] # top 
    for j in range(Ncell):
        ii,jj = nn[j] #cell ind small and large
        HR[jj,2*Nh*i+dtop:2*Nh*(i+1)+dtop,2*Nh*ii+dtop:2*Nh*(ii+1)+dtop] += HRtop[j,:,:]


###interlayer
botI = [6,7,8]
topI = [3,4,5]
#il_theta = {3:theta3,4:theta4,5:theta5,6:theta6,7:theta7,8:theta8}
il_thetab = {3:theta3b,4:theta4b,5:theta5b,6:theta6b,7:theta7b,8:theta8b}
il_thetat = {3:theta3t,4:theta4t,5:theta5t,6:theta6t,7:theta7t,8:theta8t}
# bot to top
for i in range(inTBL):
    i1,i2 = pair2[i] #bot
    for ii in botI:
        rb = atomrb[ii] + i1*a1b + i2*a2b
        thei = il_thetab[ii]
        iind = atomind[ii]
        iend = iind+atomlen[ii]
        for j in range(inTBL):
            j1,j2 = pair1[j] #top
            for jj in topI:
                rt = atomrt[jj] + j1*a1t + j2*a2t
                thej = il_thetat[jj]
                jind = atomind[jj]
                jend = jind+atomlen[jj]
                for iR in range(Ncell):
                    rR = cell[:,iR].transpose()
                    r1 = rt+rR-rb
                    r2 = rb+rR-rt
                    h1 = hpp(thei,thej,r1)
                    h2 = hpp(thej,thei,r2)
                    ib1 = 2*Nh*i + iind
                    ib2 = 2*Nh*i + iend
                    jt1 = 2*Nh*j + dtop + jind
                    jt2 = 2*Nh*j + dtop + jend
                    HR[iR,ib1:ib2,jt1:jt2] += h1
                    HR[iR,ib1+Nh:ib2+Nh,jt1+Nh:jt2+Nh] += h1
                    HR[iR,jt1:jt2,ib1:ib2] += h2
                    HR[iR,jt1+Nh:jt2+Nh,ib1+Nh:ib2+Nh] += h2

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


sx = 0.5*np.array([[0.,1.],[1.,0.]],dtype=complex)
sy = 0.5*np.array([[0.,-1.j],[1.j,0.]],dtype=complex)
sz = 0.5*np.array([[1.,0.],[0.,-1.]],dtype=complex)
sx = np.kron(sx,np.identity(Nh))
sy = np.kron(sy,np.identity(Nh))
sz = np.kron(sz,np.identity(Nh))
Isize = 2*inTBL
sx = np.kron(np.identity(Isize),sx)
sy = np.kron(np.identity(Isize),sy)
sz = np.kron(np.identity(Isize),sz)


klist = []
ST = [[],[],[],[]]
E = [[],[],[],[]]
Nocc = 42*2*inTBL
N = 9
Nk = N*N 
frac = np.linspace(0.,1.,N+1)
for x in range(N):
    for y in range(N):
        k = frac[x]*ga + frac[y]*gb
        klist.append(k)
        h = Hk(k)
        w,v = np.linalg.eigh(h)
        vh = v.T.conjugate()
        for ib,n in zip(range(4),range(Nocc-2,Nocc+2)):
            v1 = vh[n,:]
            v2 = v[:,n]
            Sx = np.linalg.multi_dot([v1,sx,v2]).real 
            Sy = np.linalg.multi_dot([v1,sy,v2]).real
            Sz = np.linalg.multi_dot([v1,sz,v2]).real
            ST[ib].append(np.array([Sx,Sy,Sz]))
            E[ib].append(w[n])


np.savez('./spintexture_TBL_p321.npz',Nocc=Nocc,klist=klist,\
         val2=ST[0],val1=ST[1],con1=ST[2],con2=ST[3],\
         Eval2=E[0],Eval1=E[1],Econ1=E[2],Econ2=E[3],ga=ga,gb=gb)


