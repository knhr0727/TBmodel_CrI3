import numpy as np
import matplotlib.pyplot as plt
import sys

def printvalue(E,X):
    for e,x in zip(E,X):
        print(" %.6f   %.10f "%(e,x))

Bohr = 0.52917721092 #Angstrom
Hartree = 27.2113860217 #eV
#unit = 5.82577568e4 # (Ohm*m)^-1 from epsilon
unit = 4.599848e6 # (Ohm*m)^-1 direct calculation of sigma

mokedata = sys.argv[1]
f = np.load(mokedata)

Sxx = f['Sxx'] 
Sxy = f['Sxy'] 
Syx = f['Syx'] 
E = f['E']
Kerr = f['Kerr'] 
f.close()

Sxxplot = Sxx.copy()
Sxyplot = Sxy
Eplot = E*Hartree

mokere = Kerr.real
mokeim = Kerr.imag

sxx = Sxxplot*unit
sxy = Sxyplot*unit
E[0] = 1. # preventing division by 0 
exx = 1.+4.*np.pi*1.j*Sxxplot/E
exy = 4.*np.pi*1.j*Sxyplot/E


def deg2mrad(x):
    return 1000.*np.pi*x/180.

def mrad2deg(x):
    return 180.*0.001*x/np.pi


fig = plt.figure('MOKE',figsize=(5,4))
ax = fig.add_subplot(111)
plt.plot(Eplot,mrad2deg(1000.*mokere),'b-')
ax.set_ylabel(r'$\theta_K$ (degrees)')
ax.set_xlabel(r'$\hbar\omega$ (eV)')
plt.xlim(0.,4.)
plt.ylim(-0.01,0.01)
#plt.legend()
plt.grid()

plt.tight_layout()




fig1 = plt.figure("epsilon_xx")
ax1 = fig1.add_subplot(111)
plt.plot(Eplot,exx.real,'b-',label=r'Re $\varepsilon_{xx}/\varepsilon_0$')
plt.plot(Eplot,exx.imag,'r-',label=r'Im $\varepsilon_{xx}/\varepsilon_0$')
plt.legend()
ax1.set_xlabel(r'$\hbar\omega$ (eV)')
plt.xlim(0,4)
plt.ylim(-2,10)
plt.grid()
plt.tight_layout()


fig2 = plt.figure("epsilon_xy")
ax2 = fig2.add_subplot(111)
plt.plot(Eplot,exy.real,'b-',label=r'Re $\varepsilon_{xy}/\varepsilon_0$')
plt.plot(Eplot,exy.imag,'r-',label=r'Im $\varepsilon_{xy}/\varepsilon_0$')
plt.legend()
ax2.set_xlabel(r'$\hbar\omega$ (eV)')
plt.xlim(0,4)
plt.ylim(-0.01,0.01)
plt.grid()
plt.tight_layout()


fig3 = plt.figure("sigma_xx")
ax3 = fig3.add_subplot(111)
plt.plot(Eplot,0.0001*sxx.real,'b-',label=r'Re $\sigma_{xx}$')
plt.plot(Eplot,0.0001*sxx.imag,'r-',label=r'Im $\sigma_{xx}$')
plt.legend()
ax3.set_ylabel(r'$10^{4} (\Omega$m$)^{-1}$')
ax3.set_xlabel(r'$\hbar\omega$ (eV)')
plt.xlim(0,4)
plt.ylim(-30,30)
plt.grid()
plt.tight_layout()


fig4 = plt.figure("sigma_xy")
ax4 = fig4.add_subplot(111)
plt.plot(Eplot,0.0001*sxy.real,'b-',label=r'Re $\sigma_{xy}$')
plt.plot(Eplot,0.0001*sxy.imag,'r-',label=r'Im $\sigma_{xy}$')
plt.legend()
ax4.set_ylabel(r'$10^{4} (\Omega$m$)^{-1}$')
ax4.set_xlabel(r'$\hbar\omega$ (eV)')
plt.xlim(0,4)
plt.ylim(-0.01,0.01)
plt.grid()
plt.tight_layout()



plt.show()


#printvalue(Eplot,mrad2deg(1000.*mokere))
#printvalue(Eplot,mrad2deg(1000.*mokeim))
#printvalue(Eplot,0.0001*sxx.real)
#printvalue(Eplot,0.0001*sxx.imag)
#printvalue(Eplot,0.0001*sxy.real)
#printvalue(Eplot,0.0001*sxy.imag)
