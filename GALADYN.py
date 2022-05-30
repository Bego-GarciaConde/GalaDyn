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
            mesh_snapshot.ac_gas = [gas_z, gas_r, gas_phi]
            mesh_snapshot.plot_acceleration_components("ac_gas")
        
        if ac_stars == 1:
            stars_z, stars_r, stars_phi = mesh_snapshot.acceleration_in_mesh_comp(comp = "stars", mode_stars = "disk", tidal = False)
            mesh_snapshot.ac_stars_disk= [stars_z, stars_r, stars_phi]
            mesh_snapshot.plot_acceleration_components("ac_stars", mode_stars="disk")

        if ac_all_sat == 1:
            mesh_snapshot.satellites_acceleration_id("all")
            mesh_snapshot.satellites_acceleration_id("arania")
            mesh_snapshot.satellites_acceleration_id("grillo")
            mesh_snapshot.satellites_acceleration_id("mosquito")



       # if ac_satellites_as_points == 1:
        # if ac_satellites_ids ==1 :
        #    satellites_acceleration_id (name, "all", plot = True)


    #Thirs step: fourier 
    fourier_ac= Fourier(snapshots_analysis=snapshots_analysis, lookback=lookback, maximo=40, minimo = 0, nbins = 40)
    if fourier_acceleration_dm ==1:
        fourier_ac.apply_fourier_accelerations(comp="dm")

    if fourier_acceleration_gas ==1:
        fourier_ac.apply_fourier_accelerations(comp="gas")
    
    if fourier_acceleration_stars ==1:
        fourier_ac.apply_fourier_accelerations(comp="stars_disk")

    if fourier_acceleration_satellites ==1:
        fourier_ac.apply_fourier_sat()

    fourier= Fourier(snapshots_analysis=snapshots_analysis, lookback=lookback)
    if fourier_density == 1:
        fourier.apply_fourier_on_disk()

    if fourier_z == 1:
        fourier.apply_fourier_on_disk(peso = "Z")

    if fourier_vz ==1: 
        fourier.apply_fourier_on_disk(peso = "VZ")

    if fourier_vr ==1: 
        fourier.apply_fourier_on_disk(peso = "Vr")

    if fourier_vphi ==1: 
        fourier.apply_fourier_on_disk(peso = "Vphi")





    #Third step: comparison


if __name__ == "__main__":
    main()
