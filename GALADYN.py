#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 24/04/2022
@author: bego
"""


from fourier_analysis_functions import *
from accelerations import * 
from config import *
#import ConfigParser


        ##################################################
        #             GALACTIC DYNAMICS                  #
        #   This code makes a complete analysis of       #
        #     the dynamics in GARROTXA model,
        #  measures accelerations of all components      #
        ##################################################


def main():

    #TO DO: INICIALIZACION, COMPROBAR QUE EXISTEN TODOS LOS DIRECTORIOS 
    #TO DO: LEER LA CONFIGURACIÓN 
#    read_parameters()
    #COMPROBAR QUE EXISTE EL CROSSMATCHS Y ARCHIVOS DE CADA COMPONENTE
    for name in snapshots_analysis:
    #First step accelerations

        if ac_dm == 1:  
            acceleration_in_mesh_comp (name,"dm",None, plot = True, tidal = False)
        if ac_gas == 1:  
            acceleration_in_mesh_comp (name,"gas",None, plot = True, tidal = False)
        if ac_stars == 1: 
            acceleration_in_mesh_comp (name,"stars","disk", plot = True, tidal = False)
            acceleration_in_mesh_comp (name,"stars","nodisk", plot = True, tidal = False)


    #Second step: fourier on accelerations
    if fourier_accelerations ==1:
        fourier_analysis_accelerations()

    #Thirs step: fourier on z and vz
    if fourier_bending ==1:
        fourier_analysis_bending ()

    #Third step: comparison

if __name__ == "__main__":
    main()