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
#     TO DO: INICIALIZACION, COMPROBAR QUE EXISTEN TODOS LOS DIRECTORIOS 
#     TO DO: LEER LA CONFIGURACIÓN 
#    read_parameters()

    if ac_calculate ==1:
        for name in snapshots_analysis:
        #First step accelerations
            snapshot = Snapshot(name)
            snapshot.load_stars()
            snapshot.load_disk()
            dfA = snapshot.filter_disk_particles()
            dfA = dfA[(dfA["Z"]<1)&(dfA["Z"]>-1)]
            mesh_snapshot = mesh(name,dfA)
            if ac_dm == 1:  
                dm_z = mesh_snapshot.acceleration_in_mesh_comp(comp = "dm", mode_stars = None, tidal = False)
                mesh_snapshot.ac_dm = dm_z
            #    mesh_snapshot.plot_acceleration_components("ac_dm")

            if ac_gas == 1:
                gas_z = mesh_snapshot.acceleration_in_mesh_comp(comp = "gas", mode_stars = None, tidal = False)
                mesh_snapshot.ac_gas = gas_z
            #   mesh_snapshot.plot_acceleration_components("ac_gas")
            
            if ac_stars == 1:
                #stars_z, stars_r, stars_phi = mesh_snapshot.acceleration_in_mesh_comp(comp = "stars", mode_stars = "disk", tidal = False)
                stars_z = mesh_snapshot.acceleration_in_mesh_comp(comp = "stars", mode_stars =None, tidal = False)
                mesh_snapshot.ac_stars= stars_z


                stars_z = mesh_snapshot.acceleration_in_mesh_comp(comp = "stars", mode_stars ="disk", tidal = False)
                #mesh_snapshot.ac_stars= stars_z
            # mesh_snapshot.plot_acceleration_components("ac_stars", mode_stars=None)

            if ac_all_sat == 1:
              #  mesh_snapshot.satellites_acceleration_id("all")
                mesh_snapshot.satellites_acceleration_id("arania")
               # mesh_snapshot.satellites_acceleration_id("grillo")
               # mesh_snapshot.satellites_acceleration_id("mosquito")
    else:
        pass


      # if ac_satellites_as_points == 1:
     #   if ac_satellites_ids ==1 :
      #     satellites_acceleration_id (name, "all", plot = True)


    #Thirs step: fourier 
    fourier_ac= Fourier(snapshots_analysis=snapshots_analysis, lookback=lookback, maximo=40, minimo = 0, nbins = 40)
    if fourier_acceleration_dm ==1:
        fourier_ac.apply_fourier_accelerations(comp="dm")

    if fourier_acceleration_gas ==1:
        fourier_ac.apply_fourier_accelerations(comp="gas")

    if fourier_acceleration_total ==1:
        fourier_ac.apply_fourier_accelerations(comp="total")
    if fourier_acceleration_stars ==1:
        fourier_ac.apply_fourier_accelerations(comp="stars")
    
    if fourier_acceleration_satellites ==1:
       # fourier_ac.apply_fourier_sat(sat_name = "all")
        fourier_ac.apply_fourier_sat(sat_name = "arania")
        #fourier_ac.apply_fourier_sat(sat_name = "grillo")
        #fourier_ac.apply_fourier_sat(sat_name = "mosquito")

    #fourier_ac_stars= Fourier(snapshots_analysis=snapshots_analysis, lookback=lookback, maximo=22, minimo = 0, nbins = 22)
   # if fourier_acceleration_stars ==1:
     #   fourier_ac_stars.apply_fourier_accelerations(comp="stars_disk", nbins=22)
   #     fourier_ac_stars.apply_fourier_accelerations(comp="stars", nbins=22)


    fourier= Fourier(snapshots_analysis=snapshots_analysis, lookback=lookback)

    if fourier_density == 1:
        fourier.apply_fourier_on_disk()

    if fourier_mass == 1:
        fourier.apply_fourier_on_disk(peso = "Mass")

    if fourier_z == 1:
        fourier.apply_fourier_on_disk(peso = "Z")

    if fourier_vz ==1: 
        fourier.apply_fourier_on_disk(peso = "VZ")

    if fourier_vr ==1: 
        fourier.apply_fourier_on_disk(peso = "Vr")

    if fourier_vphi ==1: 
        fourier.apply_fourier_on_disk(peso = "Vphi")

    if fourier_bending_breathing ==1:
        fourier.apply_fourier_on_bending_breathing()
        
    if fourier_bar == 1:
        fourier= Fourier(snapshots_analysis=snapshots_analysis, lookback=lookback, maximo = 20, minimo = 0,nbins = 40, maxmode = 4)
        fourier.apply_fourier_on_bar(stars_or_dm="dm")
        print("density Fourier applied")
        fourier.apply_fourier_on_bar(peso = "Z")
        print("Z Fourier applied")
        fourier.apply_fourier_on_bar(peso = "VZ")
        print("VZ Fourier applied")
    if fourier_dm   ==1 :
        fourier= Fourier(snapshots_analysis=snapshots_analysis, lookback=lookback,maximo=20, minimo = 0, nbins = 40)
        #fourier.apply_fourier_on_bar(stars_or_dm="dm")
        #fourier.apply_fourier_on_bar(stars_or_dm="dm",peso = None)
        fourier.apply_fourier_on_bar(stars_or_dm="dm",peso = "Z")
        #fourier.apply_fourier_on_bar(stars_or_dm="dm")




    #Third step: comparison


if __name__ == "__main__":
    main()
