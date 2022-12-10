import numpy as np
import scipy.io
from scipy import stats
import matplotlib.pylab as plt
from scipy.stats import cosine
import pandas as pd
from multiprocessing import Pool
from mpl_toolkits.axes_grid1 import make_axes_locatable
import time
from itertools import repeat
path_csv = "/mnt/usb-TOSHIBA_EXTERNAL_USB_20220124010088F-0:0-part2/snapshots_resim_new/"
path_datos = "/home/bego/GARROTXA_copia/datos_GARROTXA_resim/"
path_results = "/home/bego/GARROTXA/GalaDyn/results/"
path_crossmatch = "/home/bego/GARROTXA/satelites_crossmatch/"
path_figures_acceleration = "/home/bego/GARROTXA/aceleration_figures/"
path_figures = "/home/bego/GARROTXA/acceleration_figures/"
path_acceleration = "/home/bego/GARROTXA/acceleration/"
path_disk = "/home/bego/GARROTXA/disco/"
seconds_to_Myr = 3.15576e+16
datos_edades =  pd.read_csv(path_datos + "edades.csv", sep = ",",index_col = 0)

def force_in_mesh_gen (coord, x_sat0, y_sat0, z_sat0, mass_sat0):
    #a = GM/(r²)     <----
    X_resta = (x_sat0 - coord[0])
    Y_resta = (y_sat0 - coord[1])
    Z_resta = (z_sat0 - coord[2])
    
    dist_sat_part = np.sqrt(X_resta*X_resta + Y_resta*Y_resta + Z_resta*Z_resta)

 #   az = (mass_sat0*Z_resta/np.power(dist_sat_part,3)).astype(np.float16)

    az = mass_sat0*Z_resta/np.power(dist_sat_part,3)
   # return G*np.sum(az)/(kpc_to_km**2)
    return np.sum(az)

limit = 40 #limit for acceleration mesh
G= 1.3273e11
kpc_to_km= 3.086e16


snapshots_analysis= [ 884, 888,
890, 892, 894, 898, 900, 902, 904, 907, 908, 910, 912, 915, 916, 918, 921, 922, 924, 927, 929, 
930, 932, 934, 937,
939, 941,942, 944, 946, 948]

for name in snapshots_analysis:
    stars = pd.read_csv(path_csv+ f"{name}_stars_Rvir.csv", sep = ",")
    disk =  pd.read_csv(path_disk+ f"Stars_disco_{name}.csv", sep = ",")
    df = stars[stars['ID'].isin(disk["ID"])]
    DM = pd.read_csv(path_csv+ f"{name}_dm_Rvir.csv", sep = ",")
    gas = pd.read_csv(path_csv+ f"Gas_{name}.csv", sep = ",")


    df_filt = df[(df["R"]>5)&(df["R"]<15)].copy()
    df_sample  = df_filt.sample(frac=0.1, replace=False, random_state=1)

#--------DM-------------------
    x_sat0 = np.array(DM["X"], dtype = np.float32)
    y_sat0 = np.array(DM["Y"], dtype = np.float32)
    z_sat0 = np.array(DM["Z"], dtype = np.float32)

    mass_sat0 = np.array(DM["Mass"], dtype = np.float32)
    #mass_sat0_r = np.around(mass_sat0, decimals=0)

    ax_halo = np.zeros(len(df_sample["X"]), dtype = np.float32)
    #az_halo = np.zeros(len(df_sample["X"]))
    pool = Pool(6)

    coord = []
    for k in range(len(df_sample["X"])):
        coord.append([df_sample["X"].iloc[k], df_sample["Y"].iloc[k], df_sample["Z"].iloc[k]])
        
    a = zip(coord,repeat(x_sat0, len(df_sample["X"])), repeat(y_sat0, len(df_sample["X"])),
            repeat(z_sat0, len(df_sample["X"])),repeat(mass_sat0, len(df_sample["X"])))

    az_halo = pool.starmap(force_in_mesh_gen, a)
    pool.close()
    pool.join()


    az_halo_i =  G*np.array(az_halo)/(kpc_to_km**2)
    df_sample["az_DM"]= az_halo_i



#----GAS------------
    x_sat0 = np.array(gas["X"], dtype = np.float32)
    y_sat0 = np.array(gas["Y"], dtype = np.float32)
    z_sat0 = np.array(gas["Z"], dtype = np.float32)

    mass_sat0 = np.array(gas["Mass"], dtype = np.float32)
    #mass_sat0_r = np.around(mass_sat0, decimals=0)

    ax_halo = np.zeros(len(df_sample["X"]), dtype = np.float32)
    #az_halo = np.zeros(len(df_sample["X"]))
    pool = Pool(6)
        
    a = zip(coord,repeat(x_sat0, len(df_sample["X"])), repeat(y_sat0, len(df_sample["X"])),
            repeat(z_sat0, len(df_sample["X"])),repeat(mass_sat0, len(df_sample["X"])))

    az_halo = pool.starmap(force_in_mesh_gen, a)
    pool.close()
    pool.join()


    az_halo_i =  G*np.array(az_halo)/(kpc_to_km**2)

    df_sample["az_gas"]= az_halo_i

    df_sample.to_csv(f"/mnt/usb-TOSHIBA_EXTERNAL_USB_20220124010088F-0:0-part2/acceleration_csv/{name}_acceleration.csv")