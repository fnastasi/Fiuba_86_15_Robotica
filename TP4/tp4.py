#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 22 15:28:29 2019

@author: franco
"""

# In[]
"""
INFO IMPORTANTE
    roll: http://softmc.servotronix.com/wiki/SCARA_robot

"""
# In[]
import numpy as np
from numpy import *
import matplotlib.pyplot as plt
#from tp1_modificado import *

# In[]
a1 = 200
a2 = 200
a = 200 # Como a1 = a2 se define a = a1 = a2

tol = 1e-6

# In[]
# Defino una función que me devuelve la matriz homogenea que representa 
# la rototraslación a partir del criterio Denavit-Hartenberg 

# ángulos deben estar en radianes
def DH_hom_mat(theta,d,a,alpha):
    
    
    
    R = array(
        [[cos(theta), -sin(theta)*cos(alpha), sin(theta)*sin(alpha)],
        [sin(theta), cos(theta)*cos(alpha), -cos(theta)*sin(alpha)],
        [0,sin(alpha) , cos(alpha)]])
    P = array([[a*cos(theta),
               a*sin(theta),
               d]]).T
    A = concatenate((R,P),axis = 1) # Junta R y P -> [R | P]
    
    # Junta Ry P con el vector [0001] -> [  R    | P ]
    #                                    [ 0 0 0 | 1 ]
    
    A = concatenate((A,array([[0,0,0,1]])), axis = 0) 
    A[absolute(A)<tol] = 0
    return A

# In[]
def prob_dir_scara(t1,t2,d3,t4):
    
    [t1,t2,t4] = array([t1,t2,t4])*pi/180
    
    A_1_0  = DH_hom_mat (t1,0,a1,0)
        
    A_2_1  = DH_hom_mat (t2,0,a2, 0)
    
    A_3_2  = DH_hom_mat (0,d3,0, 0)
        
    A_4_3  = DH_hom_mat (t4,0,0, 0)
        
    A = linalg.multi_dot([A_1_0, A_2_1, A_3_2, A_4_3])

    g = 1 if sign(sin(t2)) >= 0 else -1
    
    return [A,g]


# In[]
    
def prob_inv_scara(A,g):
    
    px,py,pz = A[0:3,3]
    d3 = pz
    c2 = (px**2 + py**2 - a1**2 - a2**2 )/(2*a1*a2) 
    if abs(c2) <= 1:
        s2 = g*sqrt(1-c2**2)
        t2 = arctan2(s2,c2)
    else:
        print("Error: Posición no alcanzable")
        return [0,0,0,0]
        
    P = array([px, py])
    M = array([[a2*c2+a1, -a2*s2],[a2*s2,a2*c2+a1]])
    c1,s1 = linalg.solve(M,P)
    t1 = arctan2(s1,c1)
    
    #R_1_0  = DH_hom_mat(t1,0,a1,0)[0:3,0:3]
    #R_2_1  = DH_hom_mat(t2,0,a2, 0)[0:3,0:3]
    #R_2_0 = dot(R_1_0,R_2_1)
    #R_4_0 = A[0:3,0:3]
    #R_4_3 = dot(R_2_0.T,R_4_0)
    #t4 = arctan2(R_4_3[1,0],R_4_3[1,1])
    
    t124 = arctan2(A[1,0],A[0,0])
    t4 = t124 - t1- t2
    
    t1 = t1*180/pi
    t2 = t2*180/pi
    t4 = t4*180/pi
    
    tope_mec = t2>= 150 or t2<= -150 or d3>=-50 or d3<=-250
    if tope_mec:
        print("Atención: Para alcanzar la POSE, se sobrepasan los topes mecánicos por lo que no es alcanzable")
    return t1,t2,d3,t4

# In[]
    
def resolv_inv_scara(x,y,z,roll,g):
    R = array ([[cos(roll),-sin(roll), 0],
                [sin(roll),cos(roll),0],
                [0,0,1]])
    P = array([[x,y,z]])
    A = concatenate( (R,P.T),axis=1)
    A = concatenate((A,array([[0,0,0,1]])))
    return prob_inv_scara(A,g)


# In[]
"""    
# Prueba problema indirecto:
    x=200
    y=200
    z=-100
    roll= 0
    g=-1
    t1,t2,d3,t4 = resolv_inv_scara(x,y,z,roll,g)
"""

# In[]

# Constantes en mm y ms
    
tacc = 200
v1max = 90/1000
v2max = 180/1000
v3max = 1000/1000
v4max = 360/1000
vmax = array([v1max,v2max,v3max,v4max])
#dt = 0.1    # diferencial de tiempo en mm para hacer los gráficos.
                # Es conveniente que divida a tacc
    
# In[]

def gen_tray(POSE_vec,td_vec,dt):
    """
    Función para graficar:
        Obtiene los vectores de gen_tray_POSE que son los datos de la trayectoria y los grafica
    """
    
    global tacc
    # Si tacc no es divisible por dt, asigno a tacc el valor más cercano y
    # superior a si mismo divisible por tacc. Esto se hace para
    # que empalme bien los vectores de trayectoria de las variables articulares
    tacc = (int(tacc/dt) + 1 )*dt
    
    
    cant_pose = len(POSE_vec[:,0]) # cantidad de pose's ingresadas
    var_art_pos = zeros((4,cant_pose)) # Variable para guardar el resultado del problema inverso para cada pose
    
    # Cálculo del problema inverso para cada pose
    for id,POSE in enumerate(POSE_vec):
        var_art_pos[:,id] =  resolv_inv_scara(*POSE)
    
    
    # Cálculo de trayectoria
    var_art, var_art_der, T_vec = gen_tray_POSE(var_art_pos,td_vec,dt)
    t = arange(-tacc,np.sum(T_vec) +tacc, dt) # Variable de tiempo. Empieza de -tacc, y termina en la suma de todos los tiempos de trayectoria + tacc
    
    
    # graficos de variables articulares y velocidades
    graf_var_art(t,var_art,var_art_der)
    
    
    # graficos de trayectorias
    graf_tray(var_art)
    
    

# In[]
    
def gen_tray_POSE(var_art_pos,td,dt):
    
    
    """
    Función principal: calcula los vectores que son los valores de las variables articulares durante las trayectorias
    Debe devolver:  - La matriz con los valores de las variables articulares
                    - La matriz con los valores de las derivadas de las variables articulares
                    - vector de los valores de T (tiempo de trayectoria) en cada tramo (esto es para graficar después)
    
    
    var_art_pos es vector de los valores a donde quiere ir las variables articuladas
    var_art_pos[0] es la posición de la que sale. Debe tener minimo 2 valores   

    td tiempo deseado 
    """
    
    global tacc
    tacc = (int(tacc/dt) + 1 )*dt
    
    cant_tray = len(var_art_pos[0,:]) -1 # Cantidad de trayectorias
    
    # Inicialización de las variables articulares
    var_art = array([[],
                     [],
                     [],
                     []])
    
    var_art_der = array([[],
                         [],
                         [],
                         []])
    
    
    # Inicialización de los vectores A, B y C
    A = zeros(4)
    B = zeros(4)
    C = zeros(4)
    A = var_art_pos[:,0].reshape(4,1)
    
    # Inicialización del vector donde se guardan los valores de T de cada trayectoria
    T_vec = zeros(cant_tray)
    
    for i in arange(cant_tray):
        
        # Seteo valores de B, C y delta A y delta C en cada trayectoria
        B = var_art_pos[:,i].reshape(4,1)
        C = var_art_pos[:,i+1].reshape(4,1)
        DA = A-B
        DC = C-B
        
        # Cálculo del tiempo en que se mueven los ejes
        tmax = amax(DC.reshape(1,4)/vmax)
        T = amax([2*tacc,tmax,td[i]])
        
        
        # Cálculo de las variables articulares en la zona 1
        
        
        t = arange(-tacc,tacc,dt)
        var_art_zona1 = DC*((t+tacc)**2)/(4*T*tacc) + DA*((t-tacc)**2)/(4*tacc**2) + B
        var_art_der_zona1 = DC*(t + tacc)/ (2*T*tacc) + DA*(t-tacc)/(2*tacc**2)
        
        # Se lo agrego al vector de las variables articulares
        var_art = concatenate((var_art,var_art_zona1),axis=1)
        var_art_der = concatenate((var_art_der,var_art_der_zona1),axis=1)
        
        
        # Calculo de las variables articulares en la zona 2
        
        T = (int(T/dt) + 1 )*dt # En caso de que T no sea divisible por dt, tomo 
                                # el valor más cercano y mayor a T que sea divisible
        
        t = arange(tacc,T - tacc,dt)
        var_art_zona2 = DC*t/T + B
        var_art_der_zona2 = DC/T*ones((4,len(t)))
        
        # Se lo agrego al vector de las variables articulares
        
        var_art = concatenate((var_art,var_art_zona2),axis=1)
        var_art_der = concatenate((var_art_der,var_art_der_zona2),axis=1)
        
        # Actualizo el valor de T_vec con el valor de T usado en esta trayectoria
        T_vec[i] = T
        # Actualizo la matriz A con el valor a donde llegaron las variables articuladas
        A = var_art[:,-1].reshape(4,1)
        
        
        
    # Se hace una última zona 1 al llegar a la última POSE
    
    # Seteo valores de B, C y delta A y delta C en la última trayectoria
    B = C
    DA = A-B
    DC = C-B # Debería ser 0
    
    # Calculo de las variables articulares en la zona 1
    t = arange(-tacc,tacc,dt)
    var_art_zona1 = DC*((t+tacc)**2)/(4*T*tacc) + DA*((t-tacc)**2)/(4*tacc**2) + B
    var_art_der_zona1 = DC*(t + tacc)/ (2*T*tacc) + DA*(t-tacc)/(2*tacc**2)
    
        
    
    # Se lo agrego al vector de las variables articulares
    var_art = concatenate((var_art,var_art_zona1),axis=1)
    var_art_der = concatenate((var_art_der,var_art_der_zona1),axis=1)

    return var_art, var_art_der, T_vec
    

# In[]
    
def graf_var_art(t,var_art,var_art_der):
    
    # Gráficos de las variables articulares
    
    #### tita 1 ####
    fig, ax = plt.subplots(2,1)    
    
    ax[0].minorticks_on();ax[0].grid(which='both')
    ax[1].minorticks_on();ax[1].grid(which='both')
    
    ax[0].plot(t, var_art[0,:],label=r"$\theta_1$")
    ax[1].plot(t, var_art_der[0,:]*1e3,label=r"$\dot{\theta}_1$") # Unidades ° /s (por eso se multiplica por 1e3)
    ax[1].set_xlabel("Tiempo [ms]")
    ax[0].set_ylabel(r"$\theta_1 [°]$")
    ax[1].set_ylabel(r"$\dot{\theta}_1 [°/s]$")
    ax[0].legend(); ax[1].legend()
    #ax[0].grid() ; ax[1].grid() 
    fig.savefig("tita1.png",dpi =300)
    
    
    #### tita 2 ####
    fig, ax = plt.subplots(2,1)    
    
    ax[0].minorticks_on();ax[0].grid(which='both')
    ax[1].minorticks_on();ax[1].grid(which='both')
    
    ax[0].plot(t, var_art[1,:],label=r"$\theta_2$")
    ax[1].plot(t, var_art_der[1,:]*1e3,label=r"$\dot{\theta}_2$") # Unidades ° /s (por eso se multiplica por 1e3)
    ax[1].set_xlabel("Tiempo [ms]")
    ax[0].set_ylabel(r"$\theta_2 [°]$")
    ax[1].set_ylabel(r"$\dot{\theta}_2 [°/s]$")
    ax[0].legend(); ax[1].legend()
    #ax[0].grid() ; ax[1].grid() 
    fig.savefig("tita2.png",dpi =300)
    

    
    #### tita 2 ####
    fig, ax = plt.subplots(2,1)    
    
    ax[0].minorticks_on();ax[0].grid(which='both')
    ax[1].minorticks_on();ax[1].grid(which='both')
    
    ax[0].plot(t, var_art[2,:],label=r"$d_3$")
    ax[1].plot(t, var_art_der[2,:]*1e3,label=r"$\dot{d}_3$") # Unidades ° /s (por eso se multiplica por 1e3)
    ax[1].set_xlabel("Tiempo [ms]")
    ax[0].set_ylabel(r"$d_3 [mm]$")
    ax[1].set_ylabel(r"$\dot{d}_3 [mm/s]$")
    ax[0].legend(); ax[1].legend()
    #ax[0].grid() ; ax[1].grid() 
    fig.savefig("d3.png",dpi =300)
          
    
        #### tita 4 ####
    fig, ax = plt.subplots(2,1)    
    
    ax[0].minorticks_on();ax[0].grid(which='both')
    ax[1].minorticks_on();ax[1].grid(which='both')
    
    ax[0].plot(t, var_art[3,:],label=r"$\theta_4$")
    ax[1].plot(t, var_art_der[3,:]*1e3,label=r"$\dot{\theta}_4$") # Unidades ° /s (por eso se multiplica por 1e3)
    ax[1].set_xlabel("Tiempo [ms]")
    ax[0].set_ylabel(r"$\theta_4 [°]$")
    ax[1].set_ylabel(r"$\dot{\theta}_4 [°/s]$")
    ax[0].legend(); ax[1].legend()
    #ax[0].grid() ; ax[1].grid() 
    fig.savefig("tita4.png",dpi =300)
    
    
    
    ## Todas juntas ##
    
    fig,ax_t = plt.subplots()
    ax_t.plot(t, var_art[0,:],'r-',label=r"$\theta_1$")
    ax_t.plot(t, var_art[1,:],'b-',label=r"$\theta_2$")
    ax_t.plot(t, var_art[3,:],'g-',label=r"$\theta_4$")
    
    ax_t.set_ylim(-100,100)
    
    ax_t.minorticks_on()
    ax_t.grid(which='both')
    
    ax_d = ax_t.twinx()
    ax_d.plot(t, var_art[2,:],'k',label=r"$d_3$")
    ax_d.set_ylim(-225,-75)
    
    ax_t.set_xlabel("Tiempo [ms]")
    ax_t.set_ylabel(r"$\theta_i [°]$")
    ax_d.set_ylabel(r"$d_3 [mm]$")
    
    ax_t.legend()
    ax_d.legend(loc='upper center')
    
    fig.savefig("var_art.png",dpi =300)
    
    
    
    ##############################
    # Gráficos de las velocidades
    ##############################
    
    fig_v,ax_t_v = plt.subplots()
    ax_t_v.plot(t, var_art_der[0,:]*1e3,'r-',label=r"$\dot{\theta}_1$") # Unidades ° /s (por eso se multiplica por 1e3)
    ax_t_v.plot(t, var_art_der[1,:]*1e3,'b-',label=r"$\dot{\theta}_2$")
    ax_t_v.plot(t, var_art_der[3,:]*1e3,'g-',label=r"$\dot{\theta}_4$")
    
    ax_t_v.set_ylim(-200,200)
    
    ax_t_v.minorticks_on()
    ax_t_v.grid(which='both')
    
    ax_d_v = ax_t_v.twinx()
    ax_d_v.plot(t, var_art_der[2,:]*1e3,'k',label=r"$\dot{d}_3$") # Unidades mm/s (por eso se multiplica por 1e3)
    ax_d_v.set_ylim(-120,120)
    
    ax_t_v.set_xlabel("Tiempo [ms]")
    ax_t_v.set_ylabel(r"$\theta_i [°/s]$")
    ax_d_v.set_ylabel(r"$d_3 [mm/s]$")
    
    ax_t_v.legend()
    ax_d_v.legend()
    
    fig_v.savefig("var_art_der.png",dpi =300)


# In[]
    
def graf_tray(var_art):    
    
    #############################################################
    # Gráficos de la posición de la TCP en la terna 0 en el plano
    #############################################################
    
    TCP_x = array([])
    TCP_y = array([])
    
    # Variables para graficar la velocidad
    
    for i in arange(len(var_art.T)):
        A,g = prob_dir_scara(*var_art[:,i])
        TCP_x = append(TCP_x,A[0,3])
        TCP_y = append(TCP_y,A[1,3])        
    
    #print(var_art[:,0])    
    #print(prob_dir_scara(*var_art[:,0]))
    #print(var_art[:,-1])
    #print(prob_dir_scara(*var_art[:,-1]))
    #print(var_art[:,-1000])
    #print(prob_dir_scara(*var_art[:,-1000]))
    
    
    fig,ax = plt.subplots()
    
    ax.minorticks_on()
    ax.grid(which='both')
    
    ax.plot(TCP_x,TCP_y,label="Trayectoria TCP")
    ax.set_xlabel(r"$TCP_x$")
    ax.set_ylabel(r"$TCP_y$")    
    #ax.grid()
    ax.legend()
    fig.savefig("tray.png",dpi = 300)
    
    #fix,ax = plt.subplots()
    #print(V_x)
    #print(V_y)
    #ax.quiver(X, Y, V_x, V_y)




# In[]
"""
Cálculo del jacobiano del robot scara en la terna 0
"""    

def calc_Jac_ter0(t1,t2,d3,t4):
    J = array([[-a2*sin(t1+t2) - a1*sin(t1),    -a2*sin(t1+t2),  0 ,    0],
               [a2*cos(t1+t2) + a1*cos(t1),    a2*cos(t1+t2),  0,  0],
               [0,  0,  1,  0],
               [1,  1,  0,  1] ])
    return J    
    
