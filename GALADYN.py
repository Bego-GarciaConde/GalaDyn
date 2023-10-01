#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 24/04/2022
@author: B. García-Conde
"""


from config import *
from mesh import *
from satellites import *
from fourier import Fourier
from snapshot import Snapshot



        ##################################################
        #             GALACTIC DYNAMICS                  #
        #   This code makes a complete analysis of       #
        #     the dynamics in GARROTXA model,
        #  measures accelerations of all components      #
        ##################################################

def main ():


    if AC_CALCULATE ==1:
        for name in snapshots_analysis:
            mesh_snapshot = Mesh(name)

            if AC_DM == 1:  
                dm_z = mesh_snapshot.acceleration_in_mesh_comp(comp = "dm", mode_stars = None, tidal = False)
                mesh_snapshot.ac_dm = dm_z

            if AC_GAS == 1:
                gas_z = mesh_snapshot.acceleration_in_mesh_comp(comp = "gas", mode_stars = None, tidal = False)
                mesh_snapshot.ac_gas = gas_z
            
            if AC_STARS == 1:
                stars_z = mesh_snapshot.acceleration_in_mesh_comp(comp = "stars", mode_stars ="disk", tidal = False)

            if AC_STARS_ELLIPSOID == 1:
                stars_z = mesh_snapshot.acceleration_in_mesh_comp(comp = "stars", mode_stars ="ellipsoid", tidal = False)

            if AC_ALL_SAT == 1:
              #  mesh_snapshot.satellites_acceleration_id("all")
                mesh_snapshot.satellites_acceleration_id("arania")
               # mesh_snapshot.satellites_acceleration_id("grillo")
               # mesh_snapshot.satellites_acceleration_id("mosquito")
    else:
        pass

    #Third step: fourier 
    fourier_ac= Fourier(snapshots_analysis=snapshots_analysis, lookback=lookback, maximo=40, minimo = 0, nbins = 40)
    if FOURIER_ACCELERATION_DM ==1:
        fourier_ac.apply_fourier_accelerations(comp="dm")

    if FOURIER_ACCELERATION_GAS ==1:
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


    fourier= Fourier(snapshots_analysis=snapshots_analysis, lookback=lookback)

    if FOURIER_DENSITY == 1:
        fourier.apply_fourier_on_disk()

    if FOURIER_MASS == 1:
        fourier.apply_fourier_on_disk(peso = "Mass")

    if FOURIER_Z == 1:
        fourier.apply_fourier_on_disk(peso = "Z")

    if FOURIER_VZ ==1: 
        fourier.apply_fourier_on_disk(peso = "VZ")

    if FOURIER_VR ==1: 
        fourier.apply_fourier_on_disk(peso = "Vr")

    if FOURIER_VPHI ==1: 
        fourier.apply_fourier_on_disk(peso = "Vphi")

    if FOURIER_BENDING_BREATHING ==1:
        fourier.apply_fourier_on_bending_breathing()
        
    if FOURIER_BAR == 1:
        fourier= Fourier(snapshots_analysis=snapshots_analysis, lookback=lookback, maximo = 20, minimo = 0,nbins = 40, maxmode = 4)
        fourier.apply_fourier_on_bar(stars_or_dm="dm")
        print("density Fourier applied")
        fourier.apply_fourier_on_bar(peso = "Z")
        print("Z Fourier applied")
        fourier.apply_fourier_on_bar(peso = "VZ")
        print("VZ Fourier applied")

    if FOURIER_DM   ==1 :
        fourier= Fourier(snapshots_analysis=snapshots_analysis, lookback=lookback,maximo=20, minimo = 0, nbins = 40)
        #fourier.apply_fourier_on_bar(stars_or_dm="dm")
        #fourier.apply_fourier_on_bar(stars_or_dm="dm",peso = None)
        fourier.apply_fourier_on_bar(stars_or_dm="dm",peso = "Z")
        #fourier.apply_fourier_on_bar(stars_or_dm="dm")

    if CALCULATE_2ND_ALIGNMENT==1:
            for t,name in enumerate(snapshots_analysis):
                print(name)
                snapshot = Snapshot(name)
                snapshot.load_stars()
                snapshot.second_alignment()

    #Third step: comparison

if __name__ == "__main__":
    main()
