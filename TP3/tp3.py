#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 18:01:49 2019

@author: franco
"""

# In[]

# Imports 

import numpy as np
import matplotlib.pyplot as plt
from tp2 import  *  

"""
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
"""
# In[]

"""
Función para calcular el jacobiano (solo velocidad de traslación) en la terna 0
"""

def calc_Jac_tern0(t1,t2,t3):
    #t1 = t1*pi/180
    #t2 = t2*pi/180
    #t3 = t3*pi/180
    
    J = np.array([[-sin(t1)* (d4*sin(t2+t3) + a2*cos(t2) +a1), cos(t1)*(d4*cos(t2+t3)-a2*sin(t2)), d4*cos(t1)*cos(t2+t3)],
                  [cos(t1)* (d4*sin(t2+t3) + a2*cos(t2) +a1), sin(t1)*(d4*cos(t2+t3)-a2*sin(t2)), d4*sin(t1)*cos(t2+t3)],
                  [0, -(d4*sin(t2+t3) +a2*cos(t2)), -d4*sin(t2+t3)]])
    return J
# In[]

def graf_vel(xi, yi,zi,xf,yf,zf,dt,V_o_4):
    """
    xi,yi,zi -> posición inicial
    xf,yf,zf -> posición final
    dt -> paso
    V_0_4 -> Vector velocidad expresado en la terna 0 
    
    xi = 10; yi = -100; zi =500
    
    """
    
    # Inicialización de variales
    y = yi
    
    # Tiempo total y cantidad de iteraciones segun el tiempo y el valor de dt
    time = (yf-yi)/(V_o_4[1])
    cant_iter = int(time/dt)
    
    #Inicialización de thetas
    t1_vec = zeros(cant_iter)
    t2_vec = zeros(cant_iter)
    t3_vec = zeros(cant_iter)
    
    t1_punto_vec = zeros(cant_iter)
    t2_punto_vec = zeros(cant_iter)
    t3_punto_vec = zeros(cant_iter)
      
    
    for i in range(cant_iter):
        # Calculo de matriz homogénea A en la posición del robot
        A = np.array([[1, 0, 0, xi], [0,1,0,y] ,[0,0,1,zi],[0,0,0,1] ])
        
        # Calculo los valores de tita en la posición del robot con una configuración y los guardo en los vectores correspondientes
        t1,t2,t3, t4,t5,t6 = pos_prob_inv(A,1,1,1)*pi/180
        
        t1_vec[i] = t1
        t2_vec[i] = t2
        t3_vec[i] = t3
        
        # Calculo el Jacobiano 
        J =calc_Jac_tern0(t1,t2,t3)
        
        # Calculo los valores de tita punto
        t_vec_punto = linalg.solve(J,V_o_4)
        #t_vec_punto = linalg.inv(J)*V_o_4
        
        t1_punto_vec[i] = t_vec_punto[0]
        t2_punto_vec[i] = t_vec_punto[1]
        t3_punto_vec[i] = t_vec_punto[2]
        
         
        # Actualizó el valores de  y
        y = y + dt*V_o_4[1]
        
    
    # Se grafcan como varían los valores de theta 1 2 y 3 y el theta punto 1, 2 y 3
    time_vec = linspace(0,time,cant_iter)
    
    fig1,ax1 = plt.subplots()
    ax1.plot(time_vec*1e3,t1_vec,label = r"$\theta_1;\quad x = $"+str(xi))
    ax1.set_xlabel('Tiempo (ms)')
    ax1.set_ylabel(r'$\theta_1 (rad)$')
    ax1.grid(True)
    ax1.legend()
    
    fig2,ax2 = plt.subplots() 
    ax2.plot(time_vec*1e3,t2_vec,label = r"$\theta_2;\quad x = $"+str(xi))
    ax2.set_xlabel('Tiempo (ms)')
    ax2.set_ylabel(r'$\theta_2 (rad)$')
    ax2.grid(True)
    ax2.legend()
    
    
    fig3,ax3 = plt.subplots() 
    ax3.plot(time_vec*1e3,t3_vec,label = r"$\theta_3;\quad x = $"+str(xi))
    ax3.set_xlabel('Tiempo (ms)')
    ax3.set_ylabel(r'$\theta_3 (rad)$')
    ax3.grid(True)
    ax3.legend()
    
    fig1.savefig("Tp3_tita_1_xi_"+str(xi)+"_T0.png",dpi = 300)
    fig2.savefig("Tp3_tita_2_xi_"+str(xi)+"_T0.png",dpi = 300)
    fig3.savefig("Tp3_tita_3_xi_"+str(xi)+"_T0.png",dpi = 300)
    
    
    fig1,ax1 = plt.subplots()     
    ax1.plot(time_vec*1e3,t1_punto_vec,label = r"$\dot{\theta}_1;\quad x = $"+str(xi))
    ax1.set_xlabel('Tiempo (ms)')
    ax1.set_ylabel(r'$\dot{\theta_1} (rad/s)$')
    ax1.grid(True)
    ax1.legend()
    
    fig2,ax2 = plt.subplots() 
    ax2.plot(time_vec*1e3,t2_punto_vec,label = r"$\dot{\theta_2};\quad x = $"+str(xi))
    ax2.set_xlabel('Tiempo (ms)')
    ax2.set_ylabel(r'$\dot{\theta}_2 (rad/s)$')
    ax2.grid(True)
    ax2.legend()
    
    
    fig3,ax3 = plt.subplots() 
    ax3.plot(time_vec*1e3,t3_punto_vec,label = r"$\dot{\theta}_3;\quad x = $"+str(xi))
    ax3.set_xlabel('Tiempo (ms)')
    ax3.set_ylabel(r'$\dot{\theta}_3 (rad/s)$')
    ax3.grid(True)
    ax3.legend()
    
    fig1.savefig("Tp3_tita_1_punto_xi_"+str(xi)+"_T0.png",dpi = 300)
    fig2.savefig("Tp3_tita_2_punto_xi_"+str(xi)+"_T0.png",dpi = 300)
    fig3.savefig("Tp3_tita_3_punto_xi_"+str(xi)+"_T0.png",dpi = 300)
    
    
    #### Graficar las 3 velcidades
    fig4,ax4 = plt.subplots(3)
    ax4[0].plot(time_vec*1e3,t1_punto_vec,label = r"$\dot{\theta}_1;\quad x = $"+str(xi))
    ax4[1].plot(time_vec*1e3,t2_punto_vec,label = r"$\dot{\theta}_2;\quad x = $"+str(xi))
    ax4[2].plot(time_vec*1e3,t3_punto_vec,label = r"$\dot{\theta}_3;\quad x = $"+str(xi))
    ax4[2].set_xlabel('Tiempo (ms)')
    ax4[0].set_ylabel(r'$\dot{\theta}_1 (rad/s)$');ax4[1].set_ylabel(r'$\dot{\theta}_2 (rad/s)$');ax4[2].set_ylabel(r'$\dot{\theta}_3 (rad/s)$');
    ax4[0].grid(True);ax4[1].grid(True);ax4[2].grid(True);
    ax4[0].legend();ax4[1].legend();ax4[2].legend();
    fig4.savefig("Tp3_titas__punto_xi_"+str(xi)+"_T0.png",dpi = 600)
   
    
# In[]
    
# Creación de gráficos 

xi = 0.1; yi = -100; zi =500
xf = 0.1; yf = 100; zf =500   
dt = 1e-4    
V_o_4 = array([0,200,0])

graf_vel(xi, yi,zi,xf,yf,zf,dt,V_o_4)



# In[]

# Función que calcula el jacobiano en la terna 3
def calc_Jac_tern3(t1,t2,t3):
    #t1 = t1*pi/180
    #t2 = t2*pi/180
    #t3 = t3*pi/180
    
    J = np.array([[0, d4+a2*sin(t3),d4],
                   [d4*sin(t2+t3)+a2*cos(t2)+a1,0,0],
                   [0,-a2*cos(t3),0]])
    
    
    return J

# In[]

def graf_vel2(xi, yi,zi,xf,yf,zf,dt,V_o4):
    """
    xi,yi,zi -> posición inicial
    xf,yf,zf -> posición final
    dt -> paso
    V_0_4 -> Vector velocidad expresado en la terna 0 
    
    xi = 10; yi = -100; zi =500
    
    """
    
    # Inicialización de variales
    y = yi
    
    # Tiempo total y cantidad de iteraciones segun el tiempo y el valor de dt
    time = (yf-yi)/(V_o4[1])    
    cant_iter = int(time/dt)
    
    #Inicialización de thetas
    t1_vec = zeros(cant_iter)
    t2_vec = zeros(cant_iter)
    t3_vec = zeros(cant_iter)
    t1_punto_vec = zeros(cant_iter)
    t2_punto_vec = zeros(cant_iter)
    t3_punto_vec = zeros(cant_iter)
      
    
    for i in range(cant_iter):
        # Calculo de matriz homogénea A en la posición del robot
        A = np.array([[1, 0, 0, xi], [0,1,0,y] ,[0,0,1,zi],[0,0,0,1] ])
        
        # Calculo los valores de tita en la posición del robot con una configuración y los guardo en los vectores correspondientes
        t1,t2,t3, t4,t5,t6 = pos_prob_inv(A,1,1,1) * pi/180
        
        t1_vec[i] = t1
        t2_vec[i] = t2
        t3_vec[i] = t3
        
        # Calculo el Jacobiano en base de la terna 3
        J =calc_Jac_tern3(t1,t2,t3)
        
        #Cálculo de la matriz A_3_0  y con esto la matriz R_0_3 para obtener la velocidad en la terna 3
        A_1_0  = DH_hom_mat (t1,0,a1, -pi/2)
        A_2_1  = DH_hom_mat (t2,0,a2, 0)    
        A_3_2  = DH_hom_mat (t3,0,0, pi/2)
        
        A_3_0 = linalg.multi_dot([A_1_0, A_2_1, A_3_2]) 
        R_3_0 = A_3_0[0:3,0:3]
        V_o4_terna3 = dot(R_3_0.T,V_o4)

        
        # Calculo los valores de tita punto        
        t_vec_punto = linalg.solve(J,V_o4_terna3)
        #t_vec_punto = dot(linalg.inv(J),V_o4)
        
        t1_punto_vec[i] = t_vec_punto[0]
        t2_punto_vec[i] = t_vec_punto[1]
        t3_punto_vec[i] = t_vec_punto[2]
        
        # Actualizó el valores de  y
        y = y + dt*V_o4[1]
        
        
    # Vector de tiempo
    time_vec = linspace(0,time,cant_iter)
    
    
    
    # Se grafcan como varían los valores de theta 1 2 y 3 y el theta punto 1, 2 y 3
    fig1,ax1 = plt.subplots()     
    ax1.plot(time_vec*1e3,t1_vec,label = r"$\theta_1;\quad x = $"+str(xi))
    ax1.set_xlabel('Tiempo (ms)')
    ax1.set_ylabel(r'$\theta_1 (rad)$')
    ax1.grid(True)
    ax1.legend()
    
    fig2,ax2 = plt.subplots() 
    ax2.plot(time_vec*1e3,t2_vec,label = r"$\theta_2;\quad x = $"+str(xi))
    ax2.set_xlabel('Tiempo (ms)')
    ax2.set_ylabel(r'$\theta_2 (rad)$')
    ax2.grid(True)
    ax2.legend()
    
    
    fig3,ax3 = plt.subplots() 
    ax3.plot(time_vec*1e3,t3_vec,label = r"$\theta_3;\quad x = $"+str(xi))
    ax3.set_xlabel('Tiempo (ms)')
    ax3.set_ylabel(r'$\theta_3 (rad)$')
    ax3.grid(True)
    ax3.legend()
    
    fig1.savefig("Tp3_tita_1_xi_"+str(xi)+"_T3.png",dpi = 300)
    fig2.savefig("Tp3_tita_2_xi_"+str(xi)+"_T3.png",dpi = 300)
    fig3.savefig("Tp3_tita_3_xi_"+str(xi)+"_T3.png",dpi = 300)
    
    
    fig1,ax1 = plt.subplots()     
    ax1.plot(time_vec*1e3,t1_punto_vec,label = r"$\dot{\theta}_1;\quad x = $"+str(xi))
    ax1.set_xlabel('Tiempo (ms)')
    ax1.set_ylabel(r'$\dot{\theta_1} (rad/s)$')
    ax1.grid(True)
    ax1.legend()
    
    fig2,ax2 = plt.subplots() 
    ax2.plot(time_vec*1e3,t2_punto_vec,label = r"$\dot{\theta_2};\quad x = $"+str(xi))
    ax2.set_xlabel('Tiempo (ms)')
    ax2.set_ylabel(r'$\dot{\theta}_2 (rad/s)$')
    ax2.grid(True)
    ax2.legend()
    
    
    fig3,ax3 = plt.subplots() 
    ax3.plot(time_vec*1e3,t3_punto_vec,label = r"$\dot{\theta}_3;\quad x = $"+str(xi))
    ax3.set_xlabel('Tiempo (ms)')
    ax3.set_ylabel(r'$\dot{\theta}_3 (rad/s)$')
    ax3.grid(True)
    ax3.legend()
    
    fig1.savefig("Tp3_tita_1_punto_xi"+str(xi)+"_T3.png",dpi = 300)
    fig2.savefig("Tp3_tita_2_punto_xi"+str(xi)+"_T3.png",dpi = 300)
    fig3.savefig("Tp3_tita_3_punto_xi"+str(xi)+"_T3.png",dpi = 300)
    
    
# In[]    
    
# Ejecutar este bloque para  graficar a partir de los cálculos con la terna 3    
xi = 10; yi = -100; zi =500
xf = 10; yf = 100; zf =500   
dt = 1e-4    
V_o_4 = array([0,200,0])

graf_vel2(xi, yi,zi,xf,yf,zf,dt,V_o_4)


