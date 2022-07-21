import numpy as np
import matplotlib.pyplot as plt

Bohr = 0.52917721092 #Angstrom
Hartree = 27.2113860217 #eV

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

def dpab(l,m,n): #x,y,z / z2, x2y2, xy, xz, yz
    l,m,n = -l,-m,-n
    r3 = np.sqrt(3.)
    a = np.zeros((3,5),dtype=float)  #sigma
    a[0,0] = l*(n*n-(l*l+m*m)/2.)  #x,z2
    a[1,0] = m*(n*n-(l*l+m*m)/2.)  #y,z2
    a[2,0] = n*(n*n-(l*l+m*m)/2.)  #z,z2
    a[0,1] = r3/2.*l*(l*l-m*m)  #x,x2y2
    a[1,1] = r3/2.*m*(l*l-m*m)  #y,x2y2
    a[2,1] = r3/2.*n*(l*l-m*m)  #z,x2y2
    a[0,2] = r3*l*l*m  #x,xy
    a[1,4] = r3*m*m*n  #y,yz
    a[2,3] = r3*n*n*l  #z,zx
    a[0,4] = r3*l*m*n  #x,yz
    a[1,3] = r3*l*m*n  #y,xz
    a[2,2] = r3*l*m*n  #z,xy
    a[0,3] = r3*l*l*n  #x,zx
    a[1,2] = r3*m*m*l  #y,xy
    a[2,4] = r3*n*n*m  #z,yz
    b = np.zeros((3,5),dtype=float)  #pi
    b[0,0] = -r3*l*n*n  #x,z2
    b[1,0] = -r3*m*n*n  #y,z2
    b[2,0] = r3*n*(l*l+m*m)  #z,z2
    b[0,1] =  l*(1.-l*l+m*m)  #x,x2y2
    b[1,1] = -m*(1.+l*l-m*m)  #y,x2y2
    b[2,1] = -n*(l*l-m*m)  #z,x2y2
    b[0,2] = m*(1.-2.*l*l)  #x,xy
    b[1,4] = n*(1.-2.*m*m)  #y,yz
    b[2,3] = l*(1.-2.*n*n)  #z,zx
    b[0,4] = -2.*l*m*n  #x,yz
    b[1,3] = -2.*l*m*n  #y,xz
    b[2,2] = -2.*l*m*n  #z,xy
    b[0,3] = n*(1.-2.*l*l)  #x,zx
    b[1,2] = l*(1.-2.*m*m)  #y,xy
    b[2,4] = m*(1.-2.*n*n)  #z,yz
    return [np.transpose(a),np.transpose(b)]

def hdp(theta,r,Vlist):
    x0 = np.array([1.,0.,0.])
    y0 = np.array([0.,1.,0.])
    z0 = np.array([0.,0.,1.])
    x1 = np.array([np.cos(theta),np.sin(theta),0.])
    y1 = np.array([-1.*np.sin(theta),np.cos(theta),0.])
    z1  = np.array([0.,0.,1.])
    R = np.zeros((3,3),dtype=float)
    R[0,0],R[0,1] = np.dot(x0,x1),np.dot(x0,y1)
    R[1,0],R[1,1] = np.dot(y0,x1),np.dot(y0,y1)
    R[2,2] = 1.
    r0 = r/np.linalg.norm(r)
    l,m,n = r0[0],r0[1],r0[2]
    a,b = dpab(l,m,n)
    hdp1 = np.dot(a*Vlist[0 ] + b*Vlist[1 ],R) 
    hdp2 = np.dot(a*Vlist[2 ] + b*Vlist[3 ],R)
    hdp3 = np.dot(a*Vlist[4 ] + b*Vlist[5 ],R)
    hdp = np.zeros((5,3),dtype=float)
    hdp[0,:] = hdp1[0,:]
    hdp[1:3,:] = hdp2[1:3,:]
    hdp[3:5,:] = hdp3[3:5,:]
    return hdp

def ddabc(l,m,n): #z2, x2y2, xy, xz, yz
    r3 = np.sqrt(3.)
    a = np.zeros((5,5),dtype=float)  #sigma
    a[0,0] = (n*n-(l*l+m*m)/2.)**2  #z2,z2
    a[1,1] = 3./4.*(l*l-m*m)**2  #x2y2,x2y2
    a[1,0] = r3/2.*(l*l-m*m)*(n*n-(l*l+m*m)/2.)  #x2y2,z2
    a[0,1] = a[1,0]
    a[0,2] = r3*l*m*(n*n-(l*l+m*m)/2.)  #z2,xy
    a[0,3] = r3*n*l*(n*n-(l*l+m*m)/2.)  #z2,xz
    a[0,4] = r3*m*n*(n*n-(l*l+m*m)/2.)  #z2,yz
    a[2,0] = a[0,2]
    a[3,0] = a[0,3]
    a[4,0] = a[0,4]
    a[1,2] = 3./2.*l*m*(l*l-m*m)  #x2y2,xy
    a[1,3] = 3./2.*n*l*(l*l-m*m)  #x2y2,xz
    a[1,4] = 3./2.*m*n*(l*l-m*m)  #x2y2,yz
    a[2,1] = a[1,2]
    a[3,1] = a[1,3]
    a[4,1] = a[1,4]
    a[2,2],a[3,3],a[4,4] = 3.*l*l*m*m, 3.*n*n*l*l, 3.*m*m*n*n #xy,xz,yz
    a[2,4] = 3.*l*m*m*n  #xy,yz
    a[4,3] = 3.*m*n*n*l  #yz,xz
    a[3,2] = 3.*n*l*l*m  #xz,xy
    a[4,2] = a[2,4]
    a[3,4] = a[4,3]
    a[2,3] = a[3,2]
    b = np.zeros((5,5),dtype=float)  #pi
    b[0,0] = 3.*n*n*(l*l+m*m)  #z2,z2
    b[1,1] = (l*l+m*m-(l*l-m*m)**2)  #x2y2,x2y2
    b[1,0] = r3*n*n*(m*m-l*l)  #x2y2,z2
    b[0,1] = b[1,0]
    b[0,2] = -2.*r3*l*m*n*n  #z2,xy
    b[0,3] = r3*n*l*(l*l+m*m-n*n)  #z2,xz
    b[0,4] = r3*m*n*(l*l+m*m-n*n)  #z2,yz
    b[2,0] = b[0,2]
    b[3,0] = b[0,3]
    b[4,0] = b[0,4]
    b[1,2] = 2.*l*m*(m*m-l*l)  #x2y2,xy
    b[1,3] = n*l*(1.-2.*(l*l-m*m))  #x2y2,xz
    b[1,4] = -m*n*(1.+2.*(l*l-m*m))   #x2y2,yz
    b[2,1] = b[1,2]
    b[3,1] = b[1,3]
    b[4,1] = b[1,4]
    b[2,2] = (l*l+m*m-4.*l*l*m*m) #xy
    b[3,3] = (n*n+l*l-4.*n*n*l*l) #xz
    b[4,4] = (m*m+n*n-4.*m*m*n*n) #yz
    b[2,4] = l*n*(1.-4.*m*m)  #xy,yz
    b[4,3] = m*l*(1.-4.*n*n)  #yz,xz
    b[3,2] = n*m*(1.-4.*l*l)  #xz,xy
    b[4,2] = b[2,4]
    b[3,4] = b[4,3]
    b[2,3] = b[3,2]
    c = np.zeros((5,5),dtype=float)  #delta
    c[0,0] = 3./4.*(l*l+m*m)**2  #z2,z2
    c[1,1] = n*n+((l*l-m*m)**2)/4.  #x2y2,x2y2
    c[1,0] = r3/4.*(1.+n*n)*(l*l-m*m)  #x2y2,z2
    c[0,1] = c[1,0]
    c[0,2] = r3*l*m*(1.+n*n)/2.  #z2,xy
    c[0,3] = -r3*n*l*(l*l+m*m)/2.  #z2,xz
    c[0,4] = -r3*m*n*(l*l+m*m)/2.  #z2,yz
    c[2,0] = c[0,2]
    c[3,0] = c[0,3]
    c[4,0] = c[0,4]
    c[1,2] = l*m*(l*l-m*m)/2.  #x2y2,xy
    c[1,3] = -n*l*(1.-(l*l-m*m)/2.)  #x2y2,xz
    c[1,4] = m*n*(1.+(l*l-m*m)/2.)  #x2y2,yz
    c[2,1] = c[1,2]
    c[3,1] = c[1,3]
    c[4,1] = c[1,4]
    c[2,2] = n*n+l*l*m*m   #xy
    c[3,3] = m*m+n*n*l*l   #xz
    c[4,4] = l*l+m*m*n*n   #yz
    c[2,4] = l*n*(m*m-1.)  #xy,yz
    c[4,3] = m*l*(n*n-1.)  #yz,xz
    c[3,2] = n*m*(l*l-1.)  #xz,xy
    c[4,2] = c[2,4]
    c[3,4] = c[4,3]
    c[2,3] = c[3,2]
    return [a,b,c]

def hdd(r,Vlist,Vlist2,sign):
    r0 = r/np.linalg.norm(r)
    l,m,n = r0[0],r0[1],r0[2]
    a,b,c = ddabc(l,m,n)
    hdd = np.zeros((5,5),dtype=float)
    hdd1 = a*Vlist[0 ] + b*Vlist[1 ] + c*Vlist[2 ]
    hdd2 = a*Vlist[3 ] + b*Vlist[4 ] + c*Vlist[5 ]
    hdd3 = a*Vlist[6 ] + b*Vlist[7 ] + c*Vlist[8 ]
    hdd4 = a*Vlist[9 ] + b*Vlist[10] + c*Vlist[11]
    hdd[0,0] = hdd1[0,0]
    hdd[0,1:3] = hdd2[0,1:3]
    hdd[1:3,0] = hdd2[1:3,0]
    hdd[1:3,1:3] = hdd3[1:3,1:3]
    hdd[3:5,3:5] = hdd4[3:5,3:5]
    #chirality dependent terms
    hdd[0,3] = sign*(-1.)*m*Vlist2[0]
    hdd[3,0] = sign*(-1.)*m*Vlist2[0]
    hdd[0,4] = sign*l*Vlist2[0]
    hdd[4,0] = sign*l*Vlist2[0]
    hdd[1,3] = sign*((-2.)*m*l*l*Vlist2[1] +(l*l-m*m)*m      *Vlist2[2])  
    hdd[1,4] = sign*((-2.)*m*m*l*Vlist2[1] +(-1.)*(l*l-m*m)*l*Vlist2[2]) 
    hdd[2,3] = sign*((l*l-m*m)*l*Vlist2[1] +(2.)*m*m*l       *Vlist2[2]) 
    hdd[2,4] = sign*((l*l-m*m)*m*Vlist2[1] +(-2.)*m*l*l      *Vlist2[2]) 
    hdd[3,1] = sign*((-2.)*m*l*l*Vlist2[1] +(l*l-m*m)*m      *Vlist2[2])  
    hdd[4,1] = sign*((-2.)*m*m*l*Vlist2[1] +(-1.)*(l*l-m*m)*l*Vlist2[2]) 
    hdd[3,2] = sign*((l*l-m*m)*l*Vlist2[1] +(2.)*m*m*l       *Vlist2[2]) 
    hdd[4,2] = sign*((l*l-m*m)*m*Vlist2[1] +(-2.)*m*l*l      *Vlist2[2]) 
    return hdd

def Psgn(l1,l2,H):
    I1 = np.identity(3,dtype=float)
    I2 = np.identity(3,dtype=float)
    for i in range(3):
        I1[i,i] *= l1[i]
        I2[i,i] *= l2[i]
    return np.linalg.multi_dot([I1,H,I2]) 

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

Ez2   = -5.762742/Hartree 
Exy   = -5.3424765/Hartree 
Ezxyz = -3.841907/Hartree 
Ex    = -5.36134983/Hartree 
Ey    = -4.89230816/Hartree 
Ez    = -5.03214283/Hartree 
dxyz  = -0.856814/Hartree 
pyz   = 0.25314066/Hartree 
hpp1 = np.array(\
[[-2.46529333e-01, -9.96666667e-05,  2.53666667e-04],
 [-9.93333333e-05,  2.15872000e-01, -4.74854000e-01],
 [ 2.53666667e-04, -4.74853667e-01,  5.87302000e-01]])
hpp2 = np.array(\
[[-0.42788367, -0.06623917,  0.05873317],
 [-0.51984167,  0.04266233,  0.04036717],
 [ 0.01561217,  0.02680167, -0.10570967]])
hpp3 = np.array(\
[[ 0.02363933, -0.024254  , -0.000591  ],
 [ 0.024309  ,  0.019826  , -0.00129267],
 [ 0.00061167, -0.00129233, -0.03433567]])
hpp5 = np.array(\
[[-0.23146733, -0.03387   , -0.28284733],
 [ 0.03372133,  0.07594933,  0.103291  ],
 [ 0.28293567,  0.10311   ,  0.32818933]])
hpp6 = np.array(\
[[ 0.04747933,  0.10176733, -0.01067233],
 [-0.10169   , -0.23631467, -0.00213567],
 [ 0.01067667, -0.00213167, -0.02906533]])
hpp4 = np.transpose(hpp2)
## Simplifying / symmetrize
hpp1 = np.array(\
[[-0.2465,  0.0000,  0.0000],
 [ 0.0000,  0.2159, -0.4749],
 [ 0.0000, -0.4749,  0.5873]])/Hartree
hpp2 = np.array(\
[[-0.4279, -0.0662,  0.0587],
 [-0.5198,  0.0427,  0.0404],
 [ 0.0156,  0.0268, -0.1057]])/Hartree
hpp3 = np.array(\
[[ 0.0236, -0.0243, -0.0006],
 [ 0.0243,  0.0198, -0.0013],
 [ 0.0006, -0.0013, -0.0343]])/Hartree #rotation
hpp5 = np.array(\
[[-0.2315, -0.0338, -0.2829],
 [ 0.0338,  0.0759,  0.1032],
 [ 0.2829,  0.1032,  0.3282]])/Hartree #rotation
hpp6 = np.array(\
[[ 0.0475,  0.1017, -0.0107],
 [-0.1017, -0.2363, -0.0021],
 [ 0.0107, -0.0021, -0.0291]])/Hartree #mirror
hpp4 = np.transpose(hpp2)

# m means mirror
# no hpp1m
hpp2m = Psgn([-1,1,1],[-1,1,1],hpp2) #=3-(10)4
hpp3m = Psgn([-1,1,1],[-1,1,1],hpp3) #=3-(10)6
hpp4m = Psgn([-1,1,1],[-1,1,1],hpp4) #=3-(01)5
hpp5m = Psgn([-1,1,1],[-1,1,1],hpp5) #=3-(01)7
hpp6m = Psgn([-1,1,1],[-1,1,1],hpp6) #=3-(01)4
# r means rotation (2-fold)
hpp1r = Psgn([1,-1,-1],[1,-1,-1],hpp1) #=8-(00)3
hpp2r = Psgn([1,-1,-1],[1,-1,-1],hpp2) #=8-(00)6
hpp3r = Psgn([1,-1,-1],[1,-1,-1],hpp3) #=8-(00)4
hpp4r = Psgn([1,-1,-1],[1,-1,-1],hpp4) #=8-(00)7
hpp5r = Psgn([1,-1,-1],[1,-1,-1],hpp5) #=8-(00)5
hpp6r = Psgn([1,-1,-1],[1,-1,-1],hpp6) #=8-(1-1)6
hpp2rm = hpp2 #=8-(01)7
hpp3rm = hpp3 #=8-(01)5
hpp4rm = hpp4 #=8-(10)6
hpp5rm = hpp5 #=8-(10)4
hpp6rm = hpp6 #=8-(10)7

dpsg1_1, dppi1_1 = -4.78689331/Hartree, 0.42815845/Hartree
dpsg1_2, dppi1_2 = -0.77434264/Hartree, 0.46654074/Hartree
dpsg1_3, dppi1_3 = -1.15277907/Hartree, 0.34849574/Hartree
dpsg2_1, dppi2_1 =  0.02260047/Hartree, -0.02950182/Hartree
dpsg2_2, dppi2_2 =  0.00965214/Hartree, -0.02024346/Hartree
dpsg2_3, dppi2_3 = -0.00360908/Hartree, -0.01343560/Hartree

ddsg1,  ddpi1,  dddt1  =   -0.0641856/Hartree, 0.0, 0.0
ddsg2,  ddpi2,  dddt2  =   -0.07790281/Hartree, 0.0, 0.0
ddsg3,  ddpi3,  dddt3  =   -0.01851115/Hartree, 0.00411336/Hartree, 0.0
ddsg4,  ddpi4,  dddt4  =   0.0, -0.07902342/Hartree, -0.00081791/Hartree
ddVlist = [ddsg1,  ddpi1,  dddt1,\
           ddsg2,  ddpi2,  dddt2,\
           ddsg3,  ddpi3,  dddt3,\
           ddsg4,  ddpi4,  dddt4]  
Vdd1 = -0.02593499/Hartree
Vdd2, Vdd3 = -0.03800515/Hartree, -0.02618482/Hartree
ddVlist2 = [Vdd1,Vdd2,Vdd3]

hd = np.zeros((5,5),dtype=float)
hd[0,0],hd[1,1],hd[2,2],hd[3,3],hd[4,4] = Ez2,Exy,Exy,Ezxyz,Ezxyz
hd[1,3],hd[3,1],hd[2,4],hd[4,2] = dxyz,dxyz,-dxyz,-dxyz
hp = np.zeros((3,3),dtype=float)
hp[0,0],hp[1,1],hp[2,2],hp[1,2],hp[2,1] = Ex,Ey,Ez,pyz,pyz

def hpp(i,jR):
    if ([i,0]==jR):
        return hp
    if  (i==3):
        if  (jR==[8,0]):
            return hpp1
        elif(jR==[5,0]):
            return hpp2
        elif(jR==[7,0]):
            return hpp3
        elif(jR==[4,0]):
            return hpp4
        elif(jR==[6,0]):
            return hpp5
        elif(jR==[5,6]):
            return hpp6
        elif(jR==[4,2]):
            return hpp2m
        elif(jR==[6,2]):
            return hpp3m
        elif(jR==[5,4]):
            return hpp4m
        elif(jR==[7,4]):
            return hpp5m
        elif(jR==[4,4]):
            return hpp6m
    elif(i==4):
        if  (jR==[6,0]):
            return hpp1
        elif(jR==[3,0]):
            return hpp2
        elif(jR==[8,0]):
            return hpp3
        elif(jR==[5,0]):
            return hpp4
        elif(jR==[7,0]):
            return hpp5
        elif(jR==[3,3]):
            return hpp6
        elif(jR==[5,6]):
            return hpp2m
        elif(jR==[7,6]):
            return hpp3m
        elif(jR==[3,1]):
            return hpp4m
        elif(jR==[8,1]):
            return hpp5m
        elif(jR==[5,1]):
            return hpp6m
    elif(i==5):
        if  (jR==[7,0]):
            return hpp1
        elif(jR==[4,0]):
            return hpp2
        elif(jR==[6,0]):
            return hpp3
        elif(jR==[3,0]):
            return hpp4
        elif(jR==[8,0]):
            return hpp5
        elif(jR==[4,2]):
            return hpp6
        elif(jR==[3,3]):
            return hpp2m
        elif(jR==[8,3]):
            return hpp3m
        elif(jR==[4,7]):
            return hpp4m
        elif(jR==[6,7]):
            return hpp5m
        elif(jR==[3,7]):
            return hpp6m
    elif(i==6):
        if  (jR==[4,0]):
            return hpp1r
        elif(jR==[7,0]):
            return hpp2r
        elif(jR==[5,0]):
            return hpp3r
        elif(jR==[8,0]):
            return hpp4r
        elif(jR==[3,0]):
            return hpp5r
        elif(jR==[7,4]):
            return hpp6r
        elif(jR==[8,1]):
            return hpp2rm
        elif(jR==[3,1]):
            return hpp3rm
        elif(jR==[7,6]):
            return hpp4rm
        elif(jR==[5,6]):
            return hpp5rm
        elif(jR==[8,6]):
            return hpp6rm
    elif(i==7):
        if  (jR==[5,0]):
            return hpp1r
        elif(jR==[8,0]):
            return hpp2r
        elif(jR==[3,0]):
            return hpp3r
        elif(jR==[6,0]):
            return hpp4r
        elif(jR==[4,0]):
            return hpp5r
        elif(jR==[8,1]):
            return hpp6r
        elif(jR==[6,7]):
            return hpp2rm
        elif(jR==[4,7]):
            return hpp3rm
        elif(jR==[8,3]):
            return hpp4rm
        elif(jR==[3,3]):
            return hpp5rm
        elif(jR==[6,3]):
            return hpp6rm
    elif(i==8):
        if  (jR==[3,0]):
            return hpp1r
        elif(jR==[6,0]):
            return hpp2r
        elif(jR==[4,0]):
            return hpp3r
        elif(jR==[7,0]):
            return hpp4r
        elif(jR==[5,0]):
            return hpp5r
        elif(jR==[6,7]):
            return hpp6r
        elif(jR==[7,4]):
            return hpp2rm
        elif(jR==[5,4]):
            return hpp3rm
        elif(jR==[6,2]):
            return hpp4rm
        elif(jR==[4,2]):
            return hpp5rm
        elif(jR==[7,2]):
            return hpp6rm
# 3 [8,0],[5,0],[7,0],[4,0],[6,0],[5,6],[4,2],[6,2],[5,4],[7,4],[4,4]
# 4 [6,0],[3,0],[8,0],[5,0],[7,0],[3,3],[5,6],[7,6],[3,1],[8,1],[5,1]
# 5 [7,0],[4,0],[6,0],[3,0],[8,0],[4,2],[3,3],[8,3],[4,7],[6,7],[3,7]
# 6 [4,0],[7,0],[5,0],[8,0],[3,0],[7,4],[8,1],[3,1],[7,6],[5,6],[8,6]
# 7 [5,0],[8,0],[3,0],[6,0],[4,0],[8,1],[6,7],[4,7],[8,3],[3,3],[6,3]
# 8 [3,0],[6,0],[4,0],[7,0],[5,0],[6,7],[7,4],[5,4],[6,2],[4,2],[7,2]

cri = np.linalg.norm(Cr2-I3)+0.3
def dp_param(r):
    if (np.linalg.norm(r)<cri):
        return [dpsg1_1,dppi1_1,dpsg1_2,dppi1_2,dpsg1_3,dppi1_3]
    else:
        return [dpsg2_1,dppi2_1,dpsg2_2,dppi2_2,dpsg2_3,dppi2_3]

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

def hijR(i,jR):
    j,Ri = jR
    R = cellr(Ri)
    ri,rj = atomr[i],atomr[j]+R
    r = rj-ri
    if ((i==1)or(i==2))and((j==1)or(j==2)):
        if (i==j)and(Ri==0):
            return hd
        elif((i==1)and(j==2)):
            return hdd(r,ddVlist,ddVlist2,-1.) #LtoR
        elif((i==2)and(j==1)):
            return hdd(r,ddVlist,ddVlist2,1.) #RtoL
    if ((i!=1)and(i!=2)and(j!=1)and(j!=2)):
        return hpp(i,jR)
    if ((i==1)or(i==2))and((j!=1)and(j!=2)):
        if  (j==3)or(j==8):
            theta = theta3
        elif(j==4)or(j==6):
            theta = theta4
        elif(j==5)or(j==7):
            theta = theta5
        dpVlist = dp_param(r)
        return hdp(theta,r,dpVlist)
    if ((j==1)or(j==2))and((i!=1)and(i!=2)):
        if  (i==3)or(i==8):
            theta = theta3
        elif(i==4)or(i==6):
            theta = theta4
        elif(i==5)or(i==7):
            theta = theta5
        dpVlist = dp_param(r)
        return hdp(theta,-r,dpVlist).transpose()

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
    return ekt



#nearest neighbor/ [atom,cell]
NN = [\
[[1,0],[5,4],[7,4],[3,0],[8,0],[4,2],[6,2],[2,4],[2,0],[2,2],[8,4],[4,4],[6,0],[5,0],[7,2],[3,2]],
[[2,0],[8,0],[3,0],[6,0],[4,0],[7,0],[5,0],[1,0],[1,1],[1,3],[3,3],[6,7],[4,2],[7,4],[5,6],[8,1]],
[[3,0],[1,0],[2,0],[2,4],[1,1],[8,0],[5,0],[7,0],[4,0],[6,0],[5,6],[4,2],[6,2],[5,4],[7,4],[4,4]],
[[4,0],[1,1],[2,0],[2,1],[1,3],[6,0],[3,0],[8,0],[5,0],[7,0],[3,3],[5,6],[7,6],[3,1],[8,1],[5,1]],
[[5,0],[1,3],[2,0],[2,7],[1,0],[7,0],[4,0],[6,0],[3,0],[8,0],[4,2],[3,3],[8,3],[4,7],[6,7],[3,7]],
[[6,0],[2,0],[1,1],[1,0],[2,6],[4,0],[8,1],[3,1],[7,6],[5,6],[8,6],[7,0],[5,0],[8,0],[3,0],[7,4]],
[[7,0],[2,0],[1,3],[1,1],[2,3],[5,0],[6,7],[4,7],[8,3],[3,3],[6,3],[8,0],[3,0],[6,0],[4,0],[8,1]],
[[8,0],[2,0],[1,0],[1,3],[2,2],[3,0],[7,4],[5,4],[6,2],[4,2],[7,2],[6,0],[4,0],[7,0],[5,0],[6,7]]]

HR = np.zeros((Ncell,Nh,Nh),dtype=float)
for i,nn in zip(atomnum,NN):
    iind = atomind[i]
    iend = iind+atomlen[i]
    #print(i)
    for jR in nn:
        j,R = jR
        jind = atomind[j]
        jend = jind+atomlen[j]
        h = hijR(i,jR)
        HR[R,iind:iend,jind:jend] = h
        #print(jR,iind,iend,jind,jend)
#print(HR[0,:,:])

def Hk(k):
    exp_vec = np.exp(1.j*np.dot(k,cell))
    ekt = phase_tau(k)
    H0 = np.tensordot(exp_vec,HR,axes=((0),(0)))
    return np.linalg.multi_dot([ekt.transpose().conjugate(),H0,ekt])
    #return H0

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

G = 0.000*ga + 0.000*gb + 0.000*gc
K = (2./3.)*ga + (1./3.)*gb + 0.000*gc
M = 0.500*ga + 0.500*gb + 0.000*gc

band1 = Band(G,K,40)
band2 = Band(K,M,30)
band3 = Band(M,G,40)
PlotBand([band1,band2,band3],kticks_label=['gamma','K','M','gamma'],\
        shift=True, eV=True, EF=Efermi)
#np.savez('./banddat.1.npz',band1=band1,band2=band2,band3=band3)


np.savez('../HR.R.1.npz',HR = HR)
