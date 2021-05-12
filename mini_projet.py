# -*- coding: utf-8 -*-
"""
Created on Mon Apr 12 10:46:30 2021

@author: quentin
"""

import numpy as np
import math as m
import matplotlib.pyplot as pp
from scipy.integrate import odeint, quad

#----------------------------------------------------------
#3.2 #Benchmark utilisé
#----------------------------------------------------------

def Mat(a,b,N,K,f):
    h = (b-a)/(N-1)
    t = np.linspace(a,b,N)
    x = np.linspace(a,b,N)
    A = np.zeros((N,N))
    F = np.zeros(N)
    for i in range (N):
        F[i] = f(t[i])
        for j in range (N):
            if (j == 0) or (j == N-1):
                A[i,j] = K(x[i],t[j])
            else:
                A[i,j] = 2*K(x[i],t[j])
    M = np.eye(N)- (h/2)*A
    F = F.T
    return t,F,M

def courbe(t, U, M, F):
    sol_app = U
    sol_exa = np.linalg.solve(M, F)
    fig = pp.figure()
    fig.patch.set_facecolor('azure')
    pp.plot(t, sol_app, label = 'Solution approchée', color = "chartreuse")
    pp.plot(t, sol_exa, label = 'Solution exacte', color = "fuchsia")
    pp.legend()
    pp.title('Courbes de la solution exacte et approchée')
    pp.xlabel('t')
    pp.ylabel('solution')
    pp.grid()
    pp.show()
    
    return sol_app, sol_exa

a = -1
b = 1
N = 10

def F_benchmark(x):
    return (m.cos((m.pi*x)/2) -(8/m.pi))

def K_benchmark(x,t):
    return 2

t = np.linspace(a,b,N)
f = np.zeros((len(t),1))
U = np.zeros((len(t),1))
for i in range(len(t)):
    U[i,0] = (np.cos(0.5*np.pi*t[i]))

t, F, M = Mat(a, b, N, K_benchmark, F_benchmark)
sol_app, sol_exa = courbe(t, U, M, F)

err = np.linalg.norm((sol_app-sol_exa), 2)




#----------------------------------------------------------
#3.3 #Love électrostatique
#----------------------------------------------------------

def K_Love(x,t):
    return ((1/m.pi)*(1/(1+(x-t)**2)))

def f_Love(x):
    return 1

# Résolution avec la fonction Mat
t_Love, F_Love, M_Love = Mat(-1,1,10,K_Love,f_Love)
U_Love = np.linalg.inv(M_Love)@F_Love

# Graphique
fig = pp.figure()
fig.patch.set_facecolor('pink')
pp.plot(t_Love, U_Love, label="Solution Approchée", color = 'deeppink')
pp.title('Représentation graphique de u au cours du temps')
pp.xlabel("temps (s)")
pp.ylabel("valeurs de u")
pp.grid()
pp.show()

#----------------------------------------------------------
#4 #Circuit RLC
#----------------------------------------------------------
C = 10**-6
R = 3
L = 0.5
e = 10
h = (2-0)/200
t = np.linspace(0,2,200)

Y0 = np.zeros(2) #où i(t) vaut Y[0] et s(t) vaut Y[1]

def rlcprim(Y,t):
    C = 10**-6
    R = 3
    L = 0.5
    e = 10
    Yp0 = (e-Y[1]- R*Y[0])/L
    Yp1 = Y[0]/C
    Yp = np.array([Yp0, Yp1])
    return Yp


sol_ode = odeint(rlcprim, Y0, t)
#print(sol_ode)

#----------------------------------------------------------
#5 #Moteur à courant continu
#----------------------------------------------------------

Y0 = np.array([0,0]) #où i(t) vaut Y[0] et w(t) vaut Y[1]

C = 10**-6
R = 5
L = 50*10**-3
Ke = 0.2
Kc = 0.1
Fm = 0.01
Jm = 0.05

t = np.linspace(0,80,1000)
h = 80/8000

def u(t):
    if t>=10 and t<=50:
        u = 5
    else:
        u = 0
    return u

def moteurCC(Y,t):
    y1 = np.array([[-R/L, -Ke/L], [Kc/Jm, -Fm/Jm]])
    y2 = np.array([u(t), 0])
    Yp = np.dot(y1, Y) + y2
    return Yp

Y = odeint(moteurCC, Y0, t) #1e colonne est i(t) et 2e colonne est w(t) 


#question 5 : couple moteur
fig = pp.figure()
pp.plot(t, Kc*Y[:,0], color = 'orange')
fig.patch.set_facecolor('bisque')

pp.title("Évolution du couple moteur Cm en fonction du temps")
pp.xlabel('temps (s)')
pp.ylabel('Cm (N.m)')
pp.grid()
pp.show()

#question 6 : vitesse angulaire
fig = pp.figure()
pp.plot(t, Y[:,1], color = 'deeppink')
fig.patch.set_facecolor('pink')
pp.title("Évolution de la vitesse angulaire w en fonction du temps")
pp.xlabel('temps (s)')
pp.ylabel('w (rad/s)')
pp.grid()
pp.show()

#----------------------------------------------------------
#6 #Mouvement d'une fusée
#----------------------------------------------------------

t = np.linspace(0,160,2000)
Y0 = np.array([0,400,0])

def fusee(Y, t) :
    D = 4
    a0 = 8*10**3
    g = 9.81
    k0 = 0.1
    u = 2*10**3
    Yp = np.zeros(3)
    if (Y[1]<80):
        Y[1] = 80
        D = 0
    Yp0= ((D*u)/Y[1])-g-k0*m.exp(-Y[2]/a0)*((Y[0]**2)/Y[1])
    Yp1= -D
    Yp2= Y[0]
    Yp = np.array([Yp0, Yp1, Yp2])
    
    return Yp

Y = odeint(fusee, Y0, t)


#question 3 : vitesse
fig = pp.figure()
pp.plot(t, Y[:,0], color = 'red')
fig.patch.set_facecolor('lightcoral')

pp.title("Évolution de la vitesse v en fonction du temps")
pp.xlabel('temps (s)')
pp.ylabel('vitesse (m/s)')
pp.grid()
pp.show()

#question 4 : trajectoire
t = np.linspace(0,80,2000)
fig = pp.figure()
pp.plot(t, Y[:,2], color = 'blue')
fig.patch.set_facecolor('skyblue')

pp.title("Trajectoire de la fusée lors de la phase de propulsion")
pp.xlabel('temps (s)')
pp.ylabel('altitude (m)')
pp.grid()
pp.show()


#----------------------------------------------------------
#7 #Modèle proie-prédateur
#----------------------------------------------------------

t = np.linspace(0, 2, 500)
alpha1, beta1 = 3, 1
alpha2, beta2 = 2, 1
y1_0, y2_0 = 5, 3

#question 1 : proies sans prédateurs
y1 = []
for i in range(len(t)):
    y1.append(y1_0*m.exp(alpha1*t[i]))

fig = pp.figure()
pp.plot(t, y1, color = 'blue')
fig.patch.set_facecolor('skyblue')

pp.title("Évolution du nombre de proies en l'absence de prédateurs")
pp.xlabel('temps (années)')
pp.ylabel('nombre de proies (unité)')
pp.grid()
pp.show()

#question 2 : prédateurs sans proies
y2 = []
for i in range(len(t)):
    y2.append(y2_0*m.exp(-alpha2*t[i]))

fig = pp.figure()
pp.plot(t, y2, color = 'blue')
fig.patch.set_facecolor('skyblue')

pp.title("Évolution du nombre de prédateurs en l'absence de proies")
pp.xlabel('temps (années)')
pp.ylabel('nombre de prédateurs (unité)')
pp.grid()
pp.show()


#question 3 : méthode d'Euler proies/prédateurs par EULER
def Euler_Exp(f,y0,h):
    Y =  np.zeros(shape=(len(t),2))
    y = y0
    for i in range (0,len(t)):
        Y[i,0] = y[0]
        Y[i,1] = y[1]
        y = y + h*f(y,t[i])
    return Y

def modele(Y,t):    
    Yp0 = alpha1*Y[0]-beta1*Y[0]*Y[1]
    Yp1 = -alpha2*Y[1]+beta2*Y[0]*Y[1]
    Yp = np.array([Yp0, Yp1])
    return Yp

Y0 = np.array([5,3]) #où y1(t) vaut Y[0] et y2(t) vaut Y[1]
t = np.linspace(0,10,500)
h = 10/500

eul = Euler_Exp(modele, Y0, h)

fig = pp.figure()
fig.patch.set_facecolor('skyblue')
pp.plot(t, eul[:,0], label = 'proies')
pp.plot(t, eul[:,1], label = 'prédateurs')
pp.title('Évolution des proies et prédateurs par Euler explicite')
pp.xlabel('temps (années)')
pp.ylabel('nombre (unité)')
pp.legend()
pp.grid()
pp.show()

#question 4 : méthode d'Euler proies/prédateurs par ODEINT
ode = odeint(modele, Y0, t)

fig = pp.figure()
fig.patch.set_facecolor('skyblue')
pp.plot(t, ode[:,0], label = 'proies')
pp.plot(t, ode[:,1], label = 'prédateurs')
pp.title('Évolution des proies et prédateurs par Odeint')
pp.xlabel('temps (années)')
pp.ylabel('nombre (unité)')
pp.legend()
pp.grid()
pp.show()

#question 5 : Portrait de phase 
fig = pp.figure()
fig.patch.set_facecolor('skyblue')
pp.plot(eul[:,0], eul[:,1], label = 'par euler explicite')
pp.plot(ode[:,0], ode[:,1], label = 'par odeint')
pp.title('Portrait de phase')
pp.xlabel('nombre de proies (unité)')
pp.ylabel('Nombre de prédateurs (unité)')
pp.legend()
pp.grid()
pp.show()

#question 6 : Variation des solutions
alpha1, beta1 = 3, 1
alpha2, beta2 = 2, 1
y1_0, y2_0 = 4, 8
Y0 = np.array([y1_0,y2_0]) #où y1(t) vaut Y[0] et y2(t) vaut Y[1]
ode = odeint(modele, Y0, t)

fig = pp.figure()
fig.patch.set_facecolor('skyblue')
pp.plot(t, ode[:,0], label = 'proies')
pp.plot(t, ode[:,1], label = 'prédateurs')
pp.title('Évolution des proies et prédateurs par Odeint')
pp.xlabel('temps (années)')
pp.ylabel('nombre (unité)')
pp.legend()
pp.grid()
pp.show()

print("On constate qu'en prenant y1(0) > y2(0),\
      les courbes d'évolution des proies et des prédateurs\
          sont légèrement déphasées par rapport à celles\
              obtenues précédemment. Ce qui est logique puisque\
                  l'on débute la simulation avec plus de prédateurs\
                      que de proies cette fois-ci.")

alpha1, beta1 = 10, 1
alpha2, beta2 = 2, 1
y1_0, y2_0 = 5, 3
Y0 = np.array([y1_0,y2_0]) #où y1(t) vaut Y[0] et y2(t) vaut Y[1]
ode = odeint(modele, Y0, t)

fig = pp.figure()
fig.patch.set_facecolor('skyblue')
pp.plot(t, ode[:,0], label = 'proies')
pp.plot(t, ode[:,1], label = 'prédateurs')
pp.title('Évolution des proies et prédateurs par Odeint')
pp.xlabel('temps (années)')
pp.ylabel('nombre (unité)')
pp.legend()
pp.grid()
pp.show()

print("En prenant alpha1 plus grand, on accroît nettement\
      la différence entre le nombre de proies et de prédateurs.")
      
alpha1, beta1 = 2, 1
alpha2, beta2 = 10, 1
y1_0, y2_0 = 5, 3
Y0 = np.array([y1_0,y2_0]) #où y1(t) vaut Y[0] et y2(t) vaut Y[1]
ode = odeint(modele, Y0, t)

fig = pp.figure()
fig.patch.set_facecolor('skyblue')
pp.plot(t, ode[:,0], label = 'proies')
pp.plot(t, ode[:,1], label = 'prédateurs')
pp.title('Évolution des proies et prédateurs par Odeint')
pp.xlabel('temps (années)')
pp.ylabel('nombre (unité)')
pp.legend()
pp.grid()
pp.show()

print("En prenant alpha2 plus grand, on accroît tellement\
      la différence entre le nombre de proies et de prédateurs\
          que le nombre de proies minimum tend vers le nombre de\
              prédateurs maximum.")
              
alpha1, beta1 = 3, 10
alpha2, beta2 = 2, 1
y1_0, y2_0 = 5, 3
Y0 = np.array([y1_0,y2_0]) #où y1(t) vaut Y[0] et y2(t) vaut Y[1]
ode = odeint(modele, Y0, t)

fig = pp.figure()
fig.patch.set_facecolor('skyblue')
pp.plot(t, ode[:,0], label = 'proies')
pp.plot(t, ode[:,1], label = 'prédateurs')
pp.title('Évolution des proies et prédateurs par Odeint')
pp.xlabel('temps (années)')
pp.ylabel('nombre (unité)')
pp.legend()
pp.grid()
pp.show()

print("En prenant beta1 plus grand, on diminue la période de temps\
      qui s'écoule entre chaque cycle. Autrement dit, on retrouve\
          moins de variations sur la période de 10 ans considérée.")
          
alpha1, beta1 = 3, 1
alpha2, beta2 = 2, 10
y1_0, y2_0 = 5, 3
Y0 = np.array([y1_0,y2_0]) #où y1(t) vaut Y[0] et y2(t) vaut Y[1]
ode = odeint(modele, Y0, t)

fig = pp.figure()
fig.patch.set_facecolor('skyblue')
pp.plot(t, ode[:,0], label = 'proies')
pp.plot(t, ode[:,1], label = 'prédateurs')
pp.title('Évolution des proies et prédateurs par Odeint')
pp.xlabel('temps (années)')
pp.ylabel('nombre (unité)')
pp.legend()
pp.grid()
pp.show()

print("En prenant beta2 plus grand, on fait en sorte d'avancer\
      la variation sur la période considérée.")










