## Bibliothèques utilisées

from numpy import *
from pylab import *

from scipy.integrate import *
from scipy.integrate.odepack import odeint


## Paramètres du modèle 

Y=1
Scons=0.01
    # Paramètres opératoires
Vl=1000
Vr=1
    # Paramètres de simulation
tmax=200000
h=5
pas_integration=0.1
nb_pts=100
nb_pts2=500
    # Paramètres initiaux
Ss0=1 ; Se0=1 ; B0=1 ; Sr0=0 ; t0=0

## Fonctions 

def mu(S):
    mu=S/(1+S)

def Q(S): 
    Q=Vr*mu(S)


## Schéma Euler

def f(X,t):
    [Sr, B, Se]=X
    dSr=-B*mu(Sr)/Y + ((Q(Sr)*Se)/Vr)-((Q(Sr)*Sr)/Vr)
    dB=B*mu(Sr) - (Q(Sr)*B/Vr)
    dSe=-Q(Sr)*Se/Vl+ Q(Sr)*Sr/Vl
    
    return(dSe, dB, dSr)


def Euler(f,t0,x0,h,vt,arg):
    t=t0
    x=array(x0)
    vx=empty((0,len(x)))
    for i in vt:
        while t<i:
            x=x+h*array(f(x,t,*arg))
            t=t+h
        vx=vstack((vx,x))
    return(vx)


## Simulation du modèle

    # Détermination débit maximal

Qm=Q(Scons)
vQ=linspace(0, Qm, num=nb_pts2)

    # Représentation du débit et de la vitesse de réaction (mu)

f1=figure()
vs=linspace(0,Ss0,num=100)
plot_mu,=plot(vs,[mu(s) for s in vs],'k')
plot_D,=plot([0,Se0],[Q,Q],'b')
xlabel('substrat')
ylabel('vitesse')
legend([plot_mu,plot_D],['mu','D'],loc='best')
show(block=False)


