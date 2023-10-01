from snapshot_definition import Snapshot
from satellites_utils import Satellite
import numpy as np
import pandas as pd
import gc
from multiprocessing import Pool
from config import *
import matplotlib.pyplot as plt
from itertools import repeat
from functools import partial
import concurrent.futures
import multiprocessing
import sys
import uuid
from scipy import stats


def force_in_mesh_gen (coord, x_sat0, y_sat0, z_sat0, mass_sat0):
    """
    Calculation of the acceleration in the z direction by particles in the position
    x_sat0, y_sat0, z_sat0
    """"
    #a = GM/(r²)     <----

    X_resta = (x_sat0 - coord[0])
    Y_resta = (y_sat0 - coord[1])
    Z_resta = (z_sat0 - coord[2])
    
    dist_sat_part = np.sqrt(X_resta**2 + Y_resta**2 + Z_resta**2)


    #ax = (mass_sat0*X_resta/np.power(dist_sat_part,3)).astype(np.float64)
    #ay = (mass_sat0*Y_resta/np.power(dist_sat_part,3)).astype(np.float64)
    az = (mass_sat0*Z_resta/np.power(dist_sat_part,3)).astype(np.float64)

 #   ax = (mass_sat0*X_resta/(np.power(dist_sat_part + softening,2)*dist_sat_part)).astype(np.float64)
 #   ay = (mass_sat0*Y_resta/(np.power(dist_sat_part + softening,2)*dist_sat_part)).astype(np.float64)
 #   az = (mass_sat0*Z_resta/(np.power(dist_sat_part + softening,2)*dist_sat_part)).astype(np.float64)

   # return G*np.sum(ax)/(kpc_to_km**2),G*np.sum(ay)/(kpc_to_km**2),G*np.sum(az)/(kpc_to_km**2)
    return G*np.sum(az)/(kpc_to_km**2)

class Mesh:
    """
    Initialization of the mesh in which we will calcualte the accelerations
    """"
    def __init__(self, name):
        self.name = name
        self.lb = None
        self.x = None
        self.y = None
        self.z = None
        self.R = None
        self.phi = None
        self.ac_dm = None
        self.ac_gas = None
        self.ac_stars_disk = None
        self.ac_stars_nodisk = None
        self.ac_sat = None

        def snapshot_to_grid (nbins = 100):
   

        # snapshot_to_grid()

            mesh_x = np.zeros(nbins**2, dtype = np.float32)
            mesh_y =  np.zeros(nbins**2, dtype = np.float32)
            mesh_r = np.zeros(nbins**2, dtype = np.float32)
            mesh_z  = np.zeros(nbins**2, dtype = np.float32)
            x_array = y_array = np.linspace(-limit, +limit, nbins)

            m = 0
            for x in x_array:
                for y in y_array:
                    mesh_x[m],mesh_y[m] = x,y
                    mesh_z[m] = 0
                    mesh_r[m] = np.sqrt(x**2 +y**2)
                    m= m+1

            res = [idx for idx, val in enumerate(mesh_r) if val < limit]
            self.x = mesh_x[res]
            self.y = mesh_y[res]
            self.z = mesh_z[res]
            self.R = mesh_r[res]
            self.phi = np.mod(np.arctan2(self.y,self.x), 2*np.pi)

        snapshot_to_grid()


    def calculate_force(self, DM, tidal=False, multiprocess = False):
    """
    Multiprocessed calculation of accelerations
    """"

        az_halo = np.zeros(len(self.z), dtype = np.float64)

        x_sat0 = np.array(DM["X"], dtype = np.float32)
        y_sat0 = np.array(DM["Y"], dtype = np.float32)
        z_sat0 = np.array(DM["Z"], dtype = np.float32)

        mass_sat0 = np.array(DM["Mass"], dtype = np.float32)
        print("Calculating force")


        pool = Pool(6)

        coord = []
        for k in range(len(self.x)):
            coord.append([self.x[k], self.y[k], self.z[k]])
        a = zip(coord,repeat(x_sat0, len(self.x)), repeat(y_sat0, len(self.x)),
                repeat(z_sat0, len(self.x)),repeat(mass_sat0, len(self.x)))

        az_halo = pool.starmap(force_in_mesh_gen, a)

        pool.close()
        pool.join()

        
        data ={'X':self.x, 'Y':self.y,  'Z':self.z, 'R':self.R,'Phi':self.phi, 'az': az_halo}

        mesh_completa = pd.DataFrame(data)
        return mesh_completa


    def acceleration_in_mesh_comp (self, comp, mode_stars, tidal = False):
    """
    General function to call snapshot and proceed to calculation of accelerations
    """"
        print(self.name)
        snapshot = Snapshot(self.name)
        if comp == "dm":   
            snapshot.load_dm()
            data_acceleration = snapshot.dm
        
        if comp == "gas":
            snapshot.load_gas()
            data_acceleration= snapshot.gas

        if mode_stars =="disk":
            snapshot.load_stars()
            data_acceleration= snapshot.filter_disk_particles_by_age()
   
        elif mode_stars == "ellipsoid":
            snapshot.load_stars()
            data_acceleration = snapshot.filter_stellar_ellipsoid()

        mesh_completa = self.calculate_force(data_acceleration,multiprocess = True, tidal = tidal)        


        if mode_stars is not None:
            mesh_completa.to_csv(path_acceleration + f"mesh_aceleracion_{comp}_{mode_stars}_{self.name}_ytRS_{limit}.csv", sep = ",")
        else:
            mesh_completa.to_csv(path_acceleration + f"mesh_aceleracion_{comp}_{self.name}_ytRS_{limit}.csv" , sep = ",")
        
        return mesh_completa["az"]*seconds_to_Myr
            
    def acceleration_satellite_as_point (self, satelite, tidal = True):
    """
    Calculate acceleration if individial satellite as points
    """"
        R = np.sqrt(satelite.coord[0]**2, satelite.coord[1]**2,satelite.coord[2]**2)
        #a = GM/(r² Npart)     <----
        ax_re, ay_re, az_re = force_in_mesh_gen (satelite.coord[0], satelite.coord[1],satelite.coord[2],satelite.mass)

        if tidal == True:
            ax_0 = G* satelite.mass*satelite.coord[0]/((R**3)*(kpc_to_km**2))
            ay_0 =  G*satelite.mass*satelite.coord[1]/((R**3)*(kpc_to_km**2))
            az_0 = G*satelite.mass*satelite.coord[2]/((R**3)*(kpc_to_km**2))

            ax_re = ax_re - ax_0
            ay_re = ay_re - ay_0
            az_re = az_re - az_0
    
        a_r = ax_re*np.cos(self.phi) + ay_re*np.sin(self.phi)
        a_phi = +ax_re*(np.sin(self.phi)) - ay_re*np.cos(self.phi)
        
        return az_re, a_r, a_phi

    def satellites_acceleration_id (self, sat):
    """
    Calculate acceleration if individial satellite by separating them by ID
    """"
            
            if sat == "all":
                sat_arania =  Satellite(self.name, "arania")
                sat_grillo =  Satellite(self.name, "grillo")
                sat_mosquito =  Satellite(self.name, "mosquito")
                sat_arania.separate_stream_remanent()
                sat_grillo.separate_stream_remanent()
                sat_mosquito.separate_stream_remanent()
                DM_core_i = pd.concat ([sat_arania.core,sat_grillo.core, sat_mosquito.core])
                DM_stream_i = pd.concat ([sat_arania.stream,sat_grillo.stream, sat_mosquito.stream])

            else:
                satelite = Satellite(self.name, sat)
                satelite.separate_stream_remanent()
                DM_core_i = satelite.core
                DM_stream_i = satelite.stream

        #---------------------PROGENITOR------------------------------------
            print(f"calculating ac in core of {sat}  in {self.name}")
            mesh_core = self.calculate_force(DM_core_i, multiprocess = True)

        #---------------------STREAMS------------------------------

            mesh_stream= self.calculate_force(DM_stream_i, multiprocess = True, tidal = tidal)
 
            data ={'X':self.x, 'Y':self.y,  'Z':self.z,  'az_core': mesh_core["az"],'az_stream': mesh_stream["az"] }

            mesh_completa = pd.DataFrame(data)

            mesh_completa.to_csv(path_acceleration + f"mesh_aceleracion_{self.name}_{sat}_satellites_id_ytRS.csv", sep = ",")
            gc.collect()

