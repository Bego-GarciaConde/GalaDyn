#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 24/04/2022
@author: B. García-Conde
"""


#from fourier_analysis_functions import *
#from accelerations import * 
from config import *
from mesh_utils import *
from satellites_utils import *
from fourier_utils import Fourier
from snapshot_definition import Snapshot
#import ConfigParser


        ##################################################
        #             GALACTIC DYNAMICS                  #
        #   This code makes a complete analysis of       #
        #     the dynamics in GARROTXA model,
        #  measures accelerations of all components      #
        ##################################################

def main ():
    #TO DO: INICIALIZACION, COMPROBAR QUE EXISTEN TODOS LOS DIRECTORIOS 
    #TO DO: LEER LA CONFIGURACIÓN 
#    read_parameters()
    #COMPROBAR QUE EXISTE EL CROSSMATCHS Y ARCHIVOS DE CADA COMPONENTE
    for name in snapshots_analysis:
    #First step accelerations
        mesh_snapshot = mesh(name)
        if ac_dm == 1:  
            dm_z, dm_r, dm_phi = mesh_snapshot.acceleration_in_mesh_comp(comp = "dm", mode_stars = None, tidal = False)
            mesh_snapshot.ac_dm = [dm_z, dm_r, dm_phi]
            mesh_snapshot.plot_acceleration_components("ac_dm")

        if ac_gas == 1:
            gas_z, gas_r, gas_phi = mesh_snapshot.acceleration_in_mesh_comp(comp = "gas", mode_stars = None, tidal = False)
            mesh_snapshot.ac_dm = [gas_z, gas_r, gas_phi]
            mesh_snapshot.plot_acceleration_components("ac_gas")




        # if ac_dm == 1:  
        #     acceleration_in_mesh_comp (name, "dm",None, plot=True, tidal=False)
        # if ac_gas == 1:  
        #     acceleration_in_mesh_comp (name, "gas", None, plot=True, tidal=False)
        # if ac_stars == 1: 
        #     acceleration_in_mesh_comp (name, "stars", "disk", plot=True, tidal=False)
        #     acceleration_in_mesh_comp (name, "stars", "nodisk", plot=True, tidal=False)
         
       # if ac_satellites_as_points == 1:
        # if ac_satellites_ids ==1 :
        #    satellites_acceleration_id (name, "all", plot = True)
    #Second step: fourier on accelerations
    # if fourier_accelerations == 1:
    #     fourier= fourier()
    #     fourier.apply_fourier_on_disk()

    #Thirs step: fourier on z and vz
    fourier= Fourier(snapshots_analysis=snapshots_analysis, lookback=lookback)
    if fourier_z == 1:
        fourier.apply_fourier_on_disk(peso = "VZ")
        fourier.apply_fourier_on_disk(peso = "Z")
    if fourier_density == 1:
        fourier.apply_fourier_on_disk()


    #Third step: comparison


if __name__ == "__main__":
    main()
