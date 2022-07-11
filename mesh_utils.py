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



def force_in_mesh_gen (coord, x_sat0, y_sat0, z_sat0, mass_sat0):
    #a = GM/(r²)     <----

    X_resta = (x_sat0 - coord[0])
    Y_resta = (y_sat0 - coord[1])
    Z_resta = (z_sat0 - coord[2])
    
    dist_sat_part = np.sqrt(X_resta**2 + Y_resta**2 + Z_resta**2)


    ax = (mass_sat0*X_resta/np.power(dist_sat_part,3)).astype(np.float64)
    ay = (mass_sat0*Y_resta/np.power(dist_sat_part,3)).astype(np.float64)
    az = (mass_sat0*Z_resta/np.power(dist_sat_part,3)).astype(np.float64)

 #   ax = (mass_sat0*X_resta/(np.power(dist_sat_part + softening,2)*dist_sat_part)).astype(np.float64)
 #   ay = (mass_sat0*Y_resta/(np.power(dist_sat_part + softening,2)*dist_sat_part)).astype(np.float64)
 #   az = (mass_sat0*Z_resta/(np.power(dist_sat_part + softening,2)*dist_sat_part)).astype(np.float64)

    return G*np.sum(ax)/(kpc_to_km**2),G*np.sum(ay)/(kpc_to_km**2),G*np.sum(az)/(kpc_to_km**2)
def force_in_mesh (x_sat0, y_sat0, z_sat0, mass_sat0,  x_bin, y_bin, z_bin):
    # ToDo: meter las otras variables en la def

    #In this version:
    #a = GM/(r²)   #weighed by Npart  <----
    #F = GMm/(r²)   #weighed by Npart
  #  x_bin = mesh_x[i]
  #  y_bin = mesh_y[i]
  #  z_bin  = mesh_z[i]
    X_resta = (x_sat0 - x_bin)
    Y_resta = (y_sat0 - y_bin)
    Z_resta = (z_sat0 - z_bin)
    
    dist_sat_part = np.sqrt(X_resta**2 + Y_resta**2 + Z_resta**2)

    ax = (mass_sat0*X_resta/np.power(dist_sat_part,3)).astype(np.float64)
    ay = (mass_sat0*Y_resta/np.power(dist_sat_part,3)).astype(np.float64)
    az =(mass_sat0*Z_resta/np.power(dist_sat_part,3)).astype(np.float64)

    return G*np.sum(ax)/(kpc_to_km**2),G*np.sum(ay)/(kpc_to_km**2),G*np.sum(az)/(kpc_to_km**2)


class mesh:
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
            #limit = 20
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

    def calculate_force_sat (self, DM):
        ax_halo = np.zeros(len(self.y), dtype = np.float64)
        ay_halo = np.zeros(len(self.y), dtype = np.float64)
        az_halo = np.zeros(len(self.y), dtype = np.float64)

        x_sat0 = np.array(DM["X"], dtype = np.float64)
        y_sat0 = np.array(DM["Y"], dtype = np.float64)
        z_sat0 = np.array(DM["Z"], dtype = np.float64)
        mass_sat0 = np.array(DM["Mass"], dtype = np.float64)


        print("Calculating acceleration in the center")
        R = np.sqrt(x_sat0**2 + y_sat0**2 +z_sat0**2)
        ax_0 = np.sum(G*mass_sat0*x_sat0/((R**3)*(kpc_to_km**2)))
        ay_0 = np.sum(G*mass_sat0*y_sat0/((R**3)*(kpc_to_km**2)))
        az_0 = np.sum(G*mass_sat0*z_sat0/((R**3)*(kpc_to_km**2)))

        for i in range(len(self.y)):
        #    for i in range(20):
            ax_halo[i],ay_halo[i],az_halo[i] =force_in_mesh(x_sat0, y_sat0, z_sat0, mass_sat0,  self.x[i], self.y[i],  self.z[i])

        ax_halo = ax_halo -ax_0
        ay_halo = ay_halo -ay_0
        az_halo = az_halo -az_0
        a_r_halo = ax_halo*np.cos(self.phi) + ay_halo*np.sin(self.phi)
        a_phi_halo = ax_halo*(np.sin(self.phi)) -ay_halo*np.cos(self.phi)

        data ={'X':self.x, 'Y':self.y,  'Z':self.z, 'R':self.R,'Phi':self.phi,  
                    'ax':ax_halo, 'ay':ay_halo, 'az': az_halo,
                    'ar':a_r_halo , 'aphi':a_phi_halo}
        mesh_completa = pd.DataFrame(data)
        return mesh_completa

    def calculate_force(self, DM, tidal=False, multiprocess = False):

        ax_halo = np.zeros(len(self.x), dtype = np.float64)
        ay_halo = np.zeros(len(self.y), dtype = np.float64)
        az_halo = np.zeros(len(self.z), dtype = np.float64)

        x_sat0 = np.array(DM["X"], dtype = np.float32)
        y_sat0 = np.array(DM["Y"], dtype = np.float32)
        z_sat0 = np.array(DM["Z"], dtype = np.float32)

        mass_sat0 = np.array(DM["Mass"], dtype = np.float32)
        print("Calculating force")


        pool = Pool(6)
    # res =pool.map(force_in_mesh, [i for i in range(len(mesh_y))])
        coord = []
        for k in range(len(self.x)):
            coord.append([self.x[k], self.y[k], self.z[k]])
        a = zip(coord,repeat(x_sat0, len(self.x)), repeat(y_sat0, len(self.x)),
                repeat(z_sat0, len(self.x)),repeat(mass_sat0, len(self.x)))

        res = pool.starmap(force_in_mesh_gen, a)
        ax_halo = [item[0] for item in res]
        ay_halo = [item[1] for item in res]
        az_halo = [item[2] for item in res]
        pool.close()
        pool.join()

        if tidal ==True:
            #----Acceleration in 0,0 ------
            R = np.sqrt(x_sat0**2 + y_sat0**2 +z_sat0**2)
            ax_0 = np.sum(G*mass_sat0*x_sat0/((R**3)*(kpc_to_km**2)))
            ay_0 = np.sum(G*mass_sat0*y_sat0/((R**3)*(kpc_to_km**2)))
            az_0 = np.sum(G*mass_sat0*z_sat0/((R**3)*(kpc_to_km**2)))

          #  ax_0 = np.sum(G*mass_sat0*x_sat0/(np.power((R+softening),2)*R*(kpc_to_km**2)))
          #  ay_0 = np.sum(G*mass_sat0*y_sat0/(np.power((R+softening),2)*R*(kpc_to_km**2)))
          #  az_0 = np.sum(G*mass_sat0*z_sat0/(np.power((R+softening),2)*R*(kpc_to_km**2)))

            #------------------------------
            ax_halo = ax_halo - ax_0
            ay_halo = ay_halo - ay_0
            az_halo = az_halo - az_0


        a_r_halo = ax_halo*np.cos(self.phi) + ay_halo*np.sin(self.phi)
        a_phi_halo = ax_halo*(np.sin(self.phi)) - ay_halo*np.cos(self.phi)
        
        data ={'X':self.x, 'Y':self.y,  'Z':self.z, 'R':self.R,'Phi':self.phi,  
                    'ax':ax_halo, 'ay':ay_halo, 'az': az_halo,
                    'ar':a_r_halo , 'aphi':a_phi_halo}

        mesh_completa = pd.DataFrame(data)
        return mesh_completa


    def acceleration_in_mesh_comp (self, comp, mode_stars, tidal = False):
        print(self.name)
        #TODO: read paths
        if comp == "dm" or comp=="stars":
            DM = pd.read_csv(path_csv + f"{self.name}_{comp}_Rvir.csv", sep = ",", dtype = np.float32)
            arania = pd.read_csv(path_crossmatch + f"arania_{self.name}_crossmatch_{comp}.csv", sep = ",")
            grillo = pd.read_csv(path_crossmatch + f"grillo_{self.name}_crossmatch_{comp}.csv", sep = ",")
            mosquito = pd.read_csv(path_crossmatch + f"mosquito_{self.name}_crossmatch_{comp}.csv", sep = ",")
            print("Number of particles w/o substractig satellites ",len(DM["ID"]))
            DM = DM[~ DM['ID'].isin(arania["ID"])]
            DM = DM[~ DM['ID'].isin(grillo["ID"])]
            DM = DM[~ DM['ID'].isin(mosquito["ID"])]
            print("Number of particles after substractig satellites ", len(DM["ID"]))

        
        else:
          DM = pd.read_csv(path_csv + f"Gas_{self.name}.csv", sep = ",", dtype = np.float32)

        if mode_stars =="disk":
            disc_IDs = pd.read_csv(path_disk + f"Stars_disco_{self.name}.csv", sep = ",")
            DM  = DM[DM['ID'].isin(disc_IDs["ID"])]
        elif mode_stars == "nodisk":
            disc_IDs = pd.read_csv(path_disk + f"Stars_disco_{self.name}.csv", sep = ",")
            DM  = DM[~ DM['ID'].isin(disc_IDs["ID"])]
        
        mesh_completa = self.calculate_force(DM,multiprocess = True)        


        if mode_stars is not None:
            mesh_completa.to_csv(path_acceleration + f"mesh_aceleracion_{comp}_{mode_stars}_{self.name}_ytRS_{limit}.csv", sep = ",")
        else:
            mesh_completa.to_csv(path_acceleration + f"mesh_aceleracion_{comp}_{self.name}_ytRS_{limit}.csv" , sep = ",")
        return mesh_completa["az"], mesh_completa["ar"], mesh_completa["aphi"]

            
    def plot_acceleration_components (self, comp, rango_z=1e-14, rango_r=1e-13, rango_phi=1e-14, mode_stars=None):
        if mode_stars is not None:
            acceleration_values = getattr(self, f"{comp}_{mode_stars}" )
        else:
            acceleration_values = getattr(self, f"{comp}" )
        
        plt.style.use('dark_background')
        size = 5
        ancho = limit + 5
        fig, ax = plt.subplots(1, 3, sharex=False, sharey=True,figsize = (15,4))
        az = ax[0].scatter(self.x, self.y, marker='s', c=acceleration_values[0], 
                    cmap= "seismic", s = size, vmin =-rango_z, vmax = rango_z)

        cbar_az_ax = fig.add_axes([0.36, 0.1, 0.01,0.85 ])
        fig.colorbar(az,cbar_az_ax )

        ar = ax[1].scatter(self.x, self.y, marker='s', 
                    c=acceleration_values[1], cmap= "seismic", s =size, vmin =-rango_r, vmax =rango_r)

        cbar_ar_ax = fig.add_axes([0.63, 0.1, 0.01,0.85 ])
        fig.colorbar(ar,cbar_ar_ax )

        aphi = ax[2].scatter(self.x, self.y, marker='s',
                    c=acceleration_values[2], cmap= "seismic", s = size,vmin =-rango_phi, vmax = rango_phi)

        cbar_aphi_ax = fig.add_axes([0.91,0.1, 0.01,0.85 ])
        fig.colorbar(ar,cbar_aphi_ax )

        for k in range(3):
            ax[k].set_xlabel("X [kpc]")
            ax[k].set_ylabel("Y [kpc]")
            ax[k].set_xlim(-ancho,ancho)
            ax[k].set_ylim(-ancho,ancho)

        plt.subplots_adjust(left=0.125,bottom=0.1, right=0.9, top=0.9, wspace=0.3,hspace=0.35)
        if mode_stars is not None:
            plt.savefig(path_figures + f"acceleration_{comp}_{mode_stars}_{self.name}_ytRS_{limit}.png", format='png', dpi=150, bbox_inches='tight')
        else:
            plt.savefig(path_figures + f"acceleration_{comp}_{self.name}_ytRS_{limit}.png", format='png', dpi=150, bbox_inches='tight')
        
        gc.collect()



    def acceleration_satellite_as_point (self, satelite, tidal = True):
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


    def satellites_acceleration_id (self, sat, plot = True, rango=1e-17):
            
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
            mesh_core = self.calculate_force_sat(DM_core_i)
            #az_mean_prog = np.mean(np.abs(az_core))

        #---------------------STREAMS------------------------------

            mesh_stream= self.calculate_force_sat(DM_stream_i)
        # az_mean_stream = np.mean(np.abs(az_stream))


            data ={'X':self.x, 'Y':self.y,  'Z':self.z, 
                        'ax_core':mesh_core["ax"], 'ay_core':mesh_core["ay"], 'az_core': mesh_core["az"],
                        'ar_core':mesh_core["ar"] , 'aphi_core': mesh_core["aphi"],
                        'ax_stream':mesh_stream["ax"], 'ay_stream':mesh_stream["ay"], 'az_stream': mesh_stream["az"],
                        'ar_stream':mesh_stream["ar"] , 'aphi_stream':mesh_stream["aphi"] }

            mesh_completa = pd.DataFrame(data)

            mesh_completa.to_csv(path_acceleration + f"mesh_aceleracion_{self.name}_{sat}_satellites_id_ytRS.csv", sep = ",")
            gc.collect()

