#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 22 21:58:07 2019

@author: franco
"""

# In[]
"""
Función main del programa:
    Definir las pose
    Crear un vector con todas las pose
    Crear un vector td con N-1 componentes donde N es la cantidad de Pose's
    Elegir el tiempo dt entre las que calcula los valores de las variables articulares
    Ejecutar la función gen_tray
"""

def main():
    
    ##############
    POSE1 = array([200,200,-100,0,1])
    POSE2 = array([200,200,-200,0,-1])
    
    ##############
    POSE_vec = array([POSE1,POSE2,POSE1])
    
    ##############
    td = array([1,1]) # en ms
    
    ##############
    dt = 0.1 # en ms
    
    ##############
    gen_tray(POSE_vec,td,dt)


# In[]
    
# Imports

    
from tp4 import *  



# In[] 

if __name__ == "__main__":
    main()
    
