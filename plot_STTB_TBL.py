import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
import copy
import sys
import numpy as np
#import scipy.linalg as scipylinalg

npzfile = sys.argv[1]
f = np.load(npzfile)
Nocc = f['Nocc']
klist = f['klist']
val2 = f['val2'] 
val1 = f['val1'] 
con1 = f['con1'] 
con2 = f['con2'] 
Eval2 = f['Eval2'] 
Eval1 = f['Eval1'] 
Econ1 = f['Econ1'] 
Econ2 = f['Econ2'] 
ga = f['ga']
gb = f['gb']
f.close()

Nocc = int(Nocc)
nk = klist.shape[0]
kx,ky = [],[]
for i in range(nk):
    x,y = klist[i,:2]
    kx.append(x)
    ky.append(y)
kx = np.array(kx)
ky = np.array(ky)

val2x,val2y,val2z = [],[],[]
for i in range(nk):
    x,y,z = val2[i,:]
    val2x.append(x)
    val2y.append(y)
    val2z.append(z)
val2x = np.array(val2x)
val2y = np.array(val2y)
val2z = np.array(val2z)

val1x,val1y,val1z = [],[],[]
for i in range(nk):
    x,y,z = val1[i,:]
    val1x.append(x)
    val1y.append(y)
    val1z.append(z)
val1x = np.array(val1x)
val1y = np.array(val1y)
val1z = np.array(val1z)

con1x,con1y,con1z = [],[],[]
for i in range(nk):
    x,y,z = con1[i,:]
    con1x.append(x)
    con1y.append(y)
    con1z.append(z)
con1x = np.array(con1x)
con1y = np.array(con1y)
con1z = np.array(con1z)

con2x,con2y,con2z = [],[],[]
for i in range(nk):
    x,y,z = con2[i,:]
    con2x.append(x)
    con2y.append(y)
    con2z.append(z)
con2x = np.array(con2x)
con2y = np.array(con2y)
con2z = np.array(con2z)


###
#fig2 = plt.figure('energy surface')
#ax2 = fig2.add_subplot(111,projection='3d',aspect='equal')
#energy_surf1 = ax2.plot_surface(kX,kY,e,alpha=0.4,color='blue')
###

### scale ###
nn = 2
mm = 1
L = np.sqrt(float(nn*nn+mm*mm+mm*nn))
sclS = 1./L
#sclS = 1.
#############

Lx = ga[0]*1.05
Ly = (gb[1]-ga[1])*1.05

X1 = ga[0]
Y1 = ga[1]
X2 = gb[0]
Y2 = gb[1]
X3 = ga[0]+gb[0]
Y3 = ga[1]+gb[1]

bz1 =  (2./3.)*ga +(1./3.)*gb
bz2 =  (1./3.)*ga +(2./3.)*gb
bz3 = -(1./3.)*ga +(1./3.)*gb
bz4 = -(2./3.)*ga -(1./3.)*gb
bz5 = -(1./3.)*ga -(2./3.)*gb
bz6 =  (1./3.)*ga -(1./3.)*gb

refx = bz1[0]
refy = 3*bz2[1]

####
#Sx,Sy,Sz = val2x,val2y,val2z
Sx,Sy,Sz = val1x,val1y,val1z
#Sx,Sy,Sz = con1x,con1y,con1z
#Sx,Sy,Sz = con2x,con2y,con2z
####

fs = (6,6) #figsize
sc = 'k'

fig1 = plt.figure('spintexture',figsize=fs)
ax1 = fig1.add_subplot(111,aspect='equal')
#spintexture1 = ax1.quiver(kx,ky,sclS*val1x,sclS*val1y,color='b',angles='xy',scale_units='xy',scale=1., zorder=3.)

sz1 = ax1.scatter(kx   ,ky   ,s=15.0,c=Sz,cmap='bwr',vmin=-0.5,vmax=0.5,zorder=2.5)
sz2 = ax1.scatter(kx-X1,ky-Y1,s=15.0,c=Sz,cmap='bwr',vmin=-0.5,vmax=0.5,zorder=2.5)
sz3 = ax1.scatter(kx-X2,ky-Y2,s=15.0,c=Sz,cmap='bwr',vmin=-0.5,vmax=0.5,zorder=2.5)
sz4 = ax1.scatter(kx-X3,ky-Y3,s=15.0,c=Sz,cmap='bwr',vmin=-0.5,vmax=0.5,zorder=2.5)
fig1.colorbar(sz1)
ax1.quiver(kx   ,ky   ,sclS*Sx,sclS*Sy,color=sc,angles='xy',scale_units='xy',scale=1., zorder=3.)
ax1.quiver(kx-X1,ky-Y1,sclS*Sx,sclS*Sy,color=sc,angles='xy',scale_units='xy',scale=1., zorder=3.)
ax1.quiver(kx-X2,ky-Y2,sclS*Sx,sclS*Sy,color=sc,angles='xy',scale_units='xy',scale=1., zorder=3.)
ax1.quiver(kx-X3,ky-Y3,sclS*Sx,sclS*Sy,color=sc,angles='xy',scale_units='xy',scale=1., zorder=3.)

ref1 = ax1.quiver(refx,refy,sclS*(-0.5),0.,color='k',angles='xy',scale_units='xy',scale=1., zorder=3.)
plt.plot([bz1[0],bz2[0],bz3[0],bz4[0],bz5[0],bz6[0],bz1[0]],\
         [bz1[1],bz2[1],bz3[1],bz4[1],bz5[1],bz6[1],bz1[1]],linewidth=1.,color='black', zorder=2.)

plt.xlim([-Lx, Lx])
plt.ylim([-Ly, Ly])
ax1.set_xticks([])
ax1.set_xticklabels([])
ax1.set_yticks([])
ax1.set_yticklabels([])


plt.tight_layout()
plt.show()


