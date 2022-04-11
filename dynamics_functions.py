
import numpy as np
import pandas as pd
import gc
from multiprocessing import Pool
from config import *
import matplotlib.pyplot as plt

from config import*

G= 1.3273e11
kpc_to_km= 3.086e16

def snapshot_to_grid ():
    limit = 20
    nbins = 75
    mesh_x = np.zeros(nbins**2)
    mesh_y =  np.zeros(nbins**2)
    mesh_r = np.zeros(nbins**2)
    mesh_z  = np.zeros(nbins**2)
    x_array = y_array = np.linspace(-limit, +limit, nbins)

    m = 0
    for x in x_array:
        for y in y_array:
            mesh_x[m],mesh_y[m] = x,y
            mesh_z[m] = 0
            mesh_r[m] = np.sqrt(x**2 +y**2)
            m= m+1

    res = [idx for idx, val in enumerate(mesh_r) if val < limit]

    return mesh_x[res], mesh_y[res], mesh_z[res], mesh_r[res]

#----------------------FOR COMPONENTS-----------------------------

def force_in_mesh (x_bin, y_bin, z_bin, x_sat0, y_sat0, z_sat0, mass_sat0):
    #a = GM/(r²)     <----

    X_resta = (x_sat0 - x_bin)
    Y_resta = (y_sat0 - y_bin)
    Z_resta = (z_sat0 - z_bin)
    
    dist_sat_part = np.sqrt(X_resta**2 + Y_resta**2 + Z_resta**2)


    ax = (mass_sat0*X_resta/np.power(dist_sat_part,3)).astype(np.float64)
    ay = (mass_sat0*Y_resta/np.power(dist_sat_part,3)).astype(np.float64)
    az =(mass_sat0*Z_resta/np.power(dist_sat_part,3)).astype(np.float64)

    return G*np.sum(ax)/(kpc_to_km**2),G*np.sum(ay)/(kpc_to_km**2),G*np.sum(az)/(kpc_to_km**2)



def calculate_force(DM, mesh_x,mesh_y,mesh_z, mesh_phi,tidal=False, multiprocess = False):
    ax = np.zeros(len(mesh_y), dtype = np.float64)
    ay = np.zeros(len(mesh_y), dtype = np.float64)
    az = np.zeros(len(mesh_y), dtype = np.float64)

    x_sat0 = np.array(DM["X"], dtype = np.float64)
    y_sat0 = np.array(DM["Y"], dtype = np.float64)
    z_sat0 = np.array(DM["Z"], dtype = np.float64)
    mass_sat0 = np.array(DM["Mass"], dtype = np.float64)


    R = np.sqrt(x_sat0**2 + y_sat0**2 +z_sat0**2)
    ax_0 = np.sum(G*mass_sat0*x_sat0/((R**3)*(kpc_to_km**2)))
    ay_0 = np.sum(G*mass_sat0*y_sat0/((R**3)*(kpc_to_km**2)))
    az_0 = np.sum(G*mass_sat0*z_sat0/((R**3)*(kpc_to_km**2)))


    if multiprocess==True:    
        pool = Pool(4)
        res =pool.map(force_in_mesh, [i for i in range(len(mesh_y))])
        ax = [item[0] for item in res]
        ay= [item[1] for item in res]
        az = [item[2] for item in res]
        pool.close()
        pool.join()
    else:
        for i in range(len(mesh_y)):
            ax[i],ay[i],az[i] =force_in_mesh(mesh_x[i], mesh_y[i], mesh_z[i],x_sat0, y_sat0, z_sat0, mass_sat0)

            
    if tidal == True:
        ax = ax - ax_0
        ay = ay - ay_0
        az = az - az_0
    ar = ax*np.cos(mesh_phi) + ay*np.sin(mesh_phi)
    aphi = ax*(np.sin(mesh_phi)) -ay*np.cos(mesh_phi)
    return ax, ay,az,ar,aphi


def acceleration_in_mesh_comp (name,mode,mode_stars, path_disk, plot = True, tidal = False):
    print(name)
    #TODO: read paths
    mesh_x, mesh_y, mesh_z, mesh_r = snapshot_to_grid(name) 
    mesh_phi = np.mod(np.arctan2(mesh_y,mesh_x), 2*np.pi)
    if mode == "dm" or mode=="stars":
      DM = pd.read_csv(path_csv + "%s_%s_Rvir.csv" %(name, mode), sep = ",")
      arania = pd.read_csv(path_crossmatch +"arania_%s_crossmatch_%s.csv"%(name,mode) , sep = ",")
      grillo = pd.read_csv(path_crossmatch +"grillo_%s_crossmatch_%s.csv"%(name,mode) ,
                        sep = ",")
      mosquito = pd.read_csv(path_crossmatch+"mosquito_%s_crossmatch_%s.csv"%(name, mode),
                        sep = ",")
      print("Number of particles w/o substractig satellites ",len(DM["ID"]))
      DM = DM[~ DM['ID'].isin(arania["ID"])]
      DM = DM[~ DM['ID'].isin(grillo["ID"])]
      DM = DM[~ DM['ID'].isin(mosquito["ID"])]
      print("Number of particles after substractig satellites ", len(DM["ID"]))
    else:
      DM = pd.read_csv(path_csv +"Gas_%s.csv" %name, sep = ",", dtype = np.float32)

    if mode_stars =="disk":
      disc_IDs = pd.read_csv(path_disk + "Stars_disco_%s.csv" %(name), sep = ",")
      DM  = DM[DM['ID'].isin(disc_IDs["ID"])]
    elif mode_stars == "nodisk":
      disc_IDs = pd.read_csv(path_disk + "Stars_disco_%s.csv" %(name), sep = ",")
      DM  = DM[~ DM['ID'].isin(disc_IDs["ID"])]

    ax, ay, az, ar, aphi = calculate_force(DM,mesh_x, mesh_y, mesh_z, mesh_phi, tidal=True)

    data ={'X':mesh_x, 'Y':mesh_y,  'Z':mesh_z,'ax':ax, 'ay':ay, 'az': az,'ar':ar , 'aphi':aphi}
    mesh = pd.DataFrame(data)
    if mode_stars is not None:
      mesh.to_csv(path_acceleration + "mesh_aceleracion_%s_%s_%s_ytRS.csv" %(mode,mode_stars, name), sep = ",")
    else:
      mesh.to_csv(path_guardado + "mesh_aceleracion_%s_%s_ytRS.csv" %(mode,name), sep = ",")

    #------
    if plot == True:
      plot_acceleration_components(name, mode, mesh)


def plot_acceleration_components (name, mode, mesh_x, mesh_y, az, ar,aphi, rango_z=1e-14, rango_r=1e-13, rango_phi=1e-13):

    plt.style.use('dark_background')
    size = 5
    ancho = 25
    fig, ax = plt.subplots(1, 3, sharex=False, sharey=True,figsize = (15,4))
    az = ax[0].scatter(mesh_x, mesh_y, marker='s', c=az, 
                  cmap= "seismic", s = size, vmin =-rango_z, vmax = rango_z)

    cbar_az_ax = fig.add_axes([0.36, 0.1, 0.01,0.85 ])
    fig.colorbar(az,cbar_az_ax )

    ar = ax[1].scatter(mesh_x, mesh_y, marker='s', 
                  c=ar, cmap= "seismic", s =size, vmin =-rango_r, vmax =rango_r)

    cbar_ar_ax = fig.add_axes([0.63, 0.1, 0.01,0.85 ])
    fig.colorbar(ar,cbar_ar_ax )

    aphi = ax[2].scatter(mesh_x, mesh_y, marker='s',
                  c=aphi, cmap= "seismic", s = size,vmin =-rango_phi, vmax = rango_phi)

    cbar_aphi_ax = fig.add_axes([0.91,0.1, 0.01,0.85 ])
    fig.colorbar(ar,cbar_aphi_ax )

    for k in range(3):
      ax[k].set_xlabel("X [kpc]")
      ax[k].set_ylabel("Y [kpc]")
      ax[k].set_xlim(-ancho,ancho)
      ax[k].set_ylim(-ancho,ancho)

    plt.subplots_adjust(left=0.125,bottom=0.1, right=0.9, top=0.9, wspace=0.3,hspace=0.35)
    if mode_stars is not None:
        plt.savefig(path_guardado_figures + "acceleration_%s_%s_%s_ytRS.png"%(mode,mode_stars, name), format='png', dpi=150, bbox_inches='tight')
    else:
        plt.savefig(path_guardado_figures +"acceleration_%s_%s_ytRS.png"%(mode,name), format='png', dpi=150, bbox_inches='tight')
    
    gc.collect()

#-----------------------FOR SATELLITES AS POINTS -----------------------------------------

def acceleration_satellite_as_point (x_bin, y_bin, z_bin, sat, satelite, name, tidal = True):

    G = 1.3273e11
    kpc_to_km = 3.086e16
    sat = pd.read_csv(path_satelites + "satelite_%s_coordinates.csv"%satelite ,
                               sep = ",")
    sat_mass = pd.read_csv(path_satelites + "%s_complete_mass_table.csv"%satelite ,
                               sep = ",")
    x_sat0 = sat.loc[sat['Snapshot'] == name, 'X'].iloc[0]
    y_sat0 = sat.loc[sat['Snapshot'] == name, 'Y'].iloc[0]
    z_sat0 = sat.loc[sat['Snapshot'] == name, 'Z'].iloc[0]
    masa = sat_mass.loc[sat_mass['Snapshot'] == name, 'Mass'].iloc[0]

    #a = GM/(r² Npart)     <----
    ax_re, ay_re, az_re = force_in_mesh (x_bin, y_bin, z_bin, x_sat0, y_sat0, z_sat0, masa)
    if tidal == True:
        ax_0 = G* masa*x_sat0/((R**3)*(kpc_to_km**2))
        ay_0 =  G*masa*y_sat0/((R**3)*(kpc_to_km**2))
        az_0 = G*masa*z_sat0/((R**3)*(kpc_to_km**2))

        ax_re = ax_re - ax_0
        ay_re = ay_re - ay_0
        az_re = az_re - az_0
 
    
    a_r = ax_re*np.cos(mesh_phi) + ay_re*np.sin(mesh_phi)
    a_phi = +ax_re*(np.sin(mesh_phi)) - ay_re*np.cos(mesh_phi)
    
    return az_re,a_r, a_phi, x_sat0, y_sat0, z_sat0


def acceleration_satellites_all_plot (name):
    lookback = datos_edades.loc[datos_edades['Snapshot'] == name, 'Lookback'].iloc[0]
    fig, ax = plt.subplots(4, 3,
                           sharex=True, sharey=True,figsize = (10,12))
    fig.subplots_adjust(hspace=0.05, wspace=0.015)
    # axes are in a two-dimensional array, indexed by [row, col]
    plt.rc('font', family='serif')
    transparencia = 0.5
    tamano = 5
    satelites_lista = ["arania", "grillo", "mosquito"]
    a_z_tot = np.zeros(nbins**2)
    a_phi_tot = np.zeros(nbins**2)
    a_r_tot = np.zeros(nbins**2)
    for satelite in satelites_lista:
        if satelite == "arania":
            i = 0
            color_mark = "blue"
        if satelite == "grillo":
            i = 1
            color_mark = "red"
        if satelite == "mosquito":
            i = 2
            color_mark = "magenta"

        fuerza = acceleration_satellite_as_point(satelite, name)
        ax[0, 0].set_title("$a_{Z}$", fontsize = 15)
        ax[0, 1].set_title("$a_{R}$", fontsize = 15)
        ax[0, 2].set_title("$a_{\phi}$", fontsize = 15)
        ax[0, 1].text(-5, 85, "%.2f Gyr" %lookback, fontsize = 20)
         
        for j in range(0,3):
          #  print(j, i)
            im= ax[i, j].scatter(mesh_x, mesh_y, marker='o', c=fuerza[j],
                             cmap= "coolwarm", s = tamano,  vmin= -rango, vmax= rango,  alpha = transparencia)
            ax[i, j].plot(fuerza[3], fuerza[4], marker= "*", color = color_mark, ms = 15)
            #ax[i, j].clim(-1e-10, 1e-10)
            ax[i, j].set_xlim(-50,50)
            ax[i, j].set_ylim(-50,50)
            
   
        a_z_tot = a_z_tot + np.array(fuerza[0])
        a_r_tot = a_r_tot + np.array(fuerza[1])
        a_phi_tot = a_phi_tot + np.array(fuerza[2])

    cbar_ax_ar = fig.add_axes([0.66, 0.9, 0.235, 0.008])
    cbar_ar  = fig.colorbar(im, cax=cbar_ax_ar , orientation = "horizontal")          
    cbar_ar.ax.tick_params(labelsize= 10, top= True,bottom= False,
               labeltop=True,  labelbottom= False)
    
    ax[3, 0].scatter(mesh_x, mesh_y, marker='o', c=a_z_tot,
                     cmap= "coolwarm", s = tamano, vmin= -rango, vmax= rango, alpha = transparencia)

    ax[3, 1].scatter(mesh_x, mesh_y, marker='o', c=a_r_tot,
                     cmap= "coolwarm", s = tamano, vmin= -rango, vmax= rango,  alpha =transparencia)

    ax[3, 2].scatter(mesh_x, mesh_y, marker='o', c=a_phi_tot,
                 cmap= "coolwarm", s = tamano, vmin= -rango, vmax= rango, alpha = transparencia)

    for k in range(3):
      ax[k].set_xlabel("X [kpc]")
      ax[k].set_ylabel("Y [kpc]")
      ax[k].set_xlim(-ancho,ancho)
      ax[k].set_ylim(-ancho,ancho)
    
    plt.savefig(path_figuras + "acceleration_satelites_%s_ytRS_as_point.png"%name, format='png', dpi=150, bbox_inches='tight')

#----------------FOR SATELLITES ID------------------------------------

def separate_stream_remanent (sat, name, coord, mass):
    x_sat = coord.loc[coord['Snapshot'] == name, 'X'].iloc[0]
    y_sat = coord.loc[coord['Snapshot'] == name, 'Y'].iloc[0]
    z_sat = coord.loc[coord['Snapshot'] == name, 'Z'].iloc[0]
   # masa = mass.loc[mass['Snapshot'] == name, 'Mass'].iloc[0]
    tidal_radius = mass.loc[mass['Snapshot'] == name, 'Tidal_r'].iloc[0]
    DM= pd.read_csv(path_streams +"%s_%s_crossmatch_dm.csv" %(sat,name), sep = ",")
    stars= pd.read_csv(path_streams +"%s_%s_crossmatch_stars.csv" %(sat,name), sep = ",")
    DM = pd.concat ([DM, stars])
    R_sat = np.sqrt((DM["X"] - x_sat)**2 + (DM["Y"] - y_sat)**2 +(DM["Z"] - z_sat)**2)
    DM["R_sat"]= R_sat
    DM_core =  DM[(DM['R_sat']<= 5*tidal_radius)].copy()
    DM_stream =  DM[(DM['R_sat']> 5*tidal_radius)].copy()
    return DM_core, DM_stream


def satellites_acceleration_id (name, sat, plot = True):
    #satelites = ["arania", "grillo", "mosquito"]
    DM_core = pd.DataFrame()
    DM_stream = pd.DataFrame()
    lookback = datos_edades.loc[datos_edades['Snapshot'] == name, 'Lookback'].iloc[0]

    coord = pd.read_csv(path_datos + "satelite_%s_coordinates.csv"%sat ,sep = ",")
    mass = pd.read_csv(path_datos + "%s_complete_mass_table.csv"%sat,sep = ",")
    DM_core_i, DM_stream_i = separate_stream_remanent(sat, name,coord,mass)
    DM_core = pd.concat ([DM_core, DM_core_i])
    DM_stream = pd.concat ([DM_stream, DM_stream_i])

    print("Creating grid")
    #mesh_x, mesh_y, mesh_z, mesh_vphi, mesh_vr, mesh_vz, mesh_phi, mesh_mass, mesh_npart = snapshot_to_grid (name) 
    mesh_x, mesh_y, mesh_z, mesh_r, mesh_phi  = snapshot_to_grid(name) 

#---------------------PROGENITOR------------------------------------
    ax_core, ay_core, az_core, ar_core, aphi_core = calculate_force(DM_core_i,mesh_x, mesh_y, mesh_z, mesh_phi, tidal=True)
    #az_mean_prog = np.mean(np.abs(az_core))

#---------------------STREAMS------------------------------
    ax_stream, ay_stream, az_stream, ar_stream, aphi_stream = calculate_force(DM_stream_i,mesh_x, mesh_y, mesh_z, mesh_phi, tidal=True)
    # az_mean_stream = np.mean(np.abs(az_stream))
    data ={'X':mesh_x, 'Y':mesh_y,  'Z':mesh_z, 
                'ax_core':ax_core, 'ay_core':ay_core, 'az_core': az_core,
                'ar_core':ar_core , 'aphi_core':aphi_core,
                'ax_stream':ax_stream, 'ay_stream':ay_stream, 'az_stream': az_stream,
                'ar_stream':ar_stream , 'aphi_stream':aphi_stream }

    mesh_completa = pd.DataFrame(data)
    mesh_completa.to_csv(path_acceleration + "mesh_aceleracion_%s_%s_satellites_id_ytRS.csv" %(name, sat), sep = ",")
    gc.collect()

    if plot== True:
        #TODO Extra arguments for these plots?
        ancho = 30
        plt.style.use('dark_background')

        fig, ax = plt.subplots(2, 3, sharex=True, sharey=True,figsize = (13,8))
        fig.subplots_adjust(hspace=0.2, wspace=0.15)
        ax[0,0].scatter(mesh_x, mesh_y, marker='o', c=az_core, 
                    cmap= "coolwarm", s = 10, vmin =-rango, vmax = rango)

        ax[0,0].set_title("%.2f Gyr" %lookback, fontsize = 16)
        ar_ax = ax[0,1].scatter(mesh_x, mesh_y, marker='o', 
                    c=ar_core, cmap= "coolwarm", s = 10, vmin =-rango, vmax =rango)

        ax[0,1].set_title("Only core (part < 5 x Rtidal)", fontsize = 16)
        ax[0,2].scatter(mesh_x, mesh_y, marker='o',
                    c=aphi_core, cmap= "coolwarm", s = 10,vmin =-rango, vmax = rango)
    
        cbar_ax_ar = fig.add_axes([0.66, 0.9, 0.235, 0.008])
        cbar_ar  = fig.colorbar(ar_ax, cax=cbar_ax_ar , orientation = "horizontal")          
        cbar_ar.ax.tick_params(labelsize= 10, top= True,bottom= False,
                  labeltop=True,  labelbottom= False)

        az = ax[1,0].scatter(mesh_x, mesh_y, marker='o', c=az_stream, 
                    cmap= "coolwarm", s = 10, vmin =-rango, vmax = rango)

        ar = ax[1,1].scatter(mesh_x, mesh_y, marker='o', 
                    c=ar_stream, cmap= "coolwarm", s = 10, vmin =-rango, vmax =rango)


        ax[1,1].set_title("Only streams (part > 5 x Rtidal)", fontsize = 16)
        aphi = ax[1,2].scatter(mesh_x, mesh_y, marker='o',
                    c=aphi_stream, cmap= "coolwarm", s = 10,vmin =-rango, vmax = rango)

        cbar_ax_ar = fig.add_axes([0.66, 0.9, 0.235, 0.008])
        cbar_ar  = fig.colorbar(ar, cax=cbar_ax_ar , orientation = "horizontal")          
        cbar_ar.ax.tick_params(labelsize= 10, top= True,bottom= False,
                  labeltop=True,  labelbottom= False)


        for k in range(3):
          for m in range(2):
            ax[m,k].set_xlabel("X [kpc]")
            ax[m,k].set_ylabel("Y [kpc]")
            ax[m,k].set_xlim(-ancho,ancho)
            ax[m,k].set_ylim(-ancho,ancho)

        plt.savefig(path_figures +  "%s_core_streams_%s.png"%(sat,name), dpi = 100)
        plt.close()


#-----------------FOURIER -----------------------------------------------


def fourier(X,Y, peso=None, maximo = 22, minimo = 0,nbins = 22 ):
    maximo = 22
    minimo =0
    #paso = 1
    #nbins = (maximo-minimo)//paso
    AA =  np.zeros((3,nbins))
    armangle =  np.zeros((3,nbins))
    steparm = (maximo-minimo)/nbins
    radarm = np.arange(minimo,maximo + steparm, steparm )
    #nparticles = len(X)
    #peso_final = np.mean(peso)
    A = np.zeros(6)
    B = np.zeros(6)

    #dd = X**2 + Y**2
    dd = np.sqrt(X**2 + Y**2)
    indr=np.digitize(dd,radarm)-1
    nparticles = np.zeros(nbins)
    for m in range(0,5):
       for i in range(1,nbins):
            X_i=X[indr==i] 
            Y_i=Y[indr==i]

            A[m] = 0
            B[m] = 0
            a = []
            b = []
            a = np.arctan2(Y_i, X_i)
            b= np.arcsin(Y_i/np.sqrt(Y_i**2 + X_i**2))
            if peso is None:
                if m == 0:
                    A[m] = np.sum(np.cos(m*a))
                    B[m] = np.sum(np.sin(m*a))
                else :
                    A[m] = np.sum(2*np.cos(m*a))
                    B[m] = np.sum(2*np.sin(m*a))

            else:
                peso_i=peso[indr==i]
                if m ==0:
                    A[m] = np.sum(peso_i*np.cos(m*a))
                    B[m] = np.sum(peso_i*np.sin(m*a))
                else :
                    A[m] = np.sum(2*peso_i*np.cos(m*a))
                    B[m] = np.sum(2*peso_i*np.sin(m*a))
            
            AA[m,i] = np.sqrt(A[m]**2+ B[m]**2)
            if m > 0:
                armangle[m,i] = np.arctan2(B[m],A[m])
            elif m == 0:
                armangle[m,i] = 0
            if m ==0:
                nparticles[i]= len(a)


    Amp0 = AA[0,:]
    Amp1 = AA[1,:]
    Amp2 = AA[2,:]
    return radarm[:-1], nparticles, Amp0, Amp1, Amp2


def apply_fourier(modo, nbins = 22):
    datos_ = np.zeros((1,5))
    datos = datos_[1:]
    for t,name in enumerate(snapshots_analysis):
        df = pd.read_csv(path_acceleration + "mesh_aceleracion_%s_%s_ytRS.csv"%(modo, name) ,sep = ",")
        Rcenters, Amp0, Amp1, Amp2 =fourier(df["X"],df["Y"],peso=df["az"])
        for i in range(nbins):
            datos = np.append(datos,[[snapshots_analysis[t],Rcenters[i],Amp0[i],
            Amp1[i], Amp2[i]]], axis = 0)
    return datos

def apply_fourier_sat(nbins = 22):
    datos_c_ = np.zeros((1,5))
    datos_c = datos_c_[1:]
    datos_s_ = np.zeros((1,5))
    datos_s = datos_s_[1:]
    for t,name in enumerate(snapshots_analysis):
        df = pd.read_csv(path_acceleration + "mesh_aceleracion_%s_all_satellites_id_ytRS.csv"%(name) ,sep = ",")
        Rcenters_s, Amp0_s, Amp1_s, Amp2_s =fourier(df["X"],df["Y"],peso=df["az_stream"])
        Rcenters_c, Amp0_c, Amp1_c, Amp2_c=fourier(df["X"],df["Y"],peso=df["az_core"])
        for i in range(nbins):
        #print(len(Amp0[i]))
            datos_c = np.append(datos_c,[[snapshots_analysis[t],Rcenters_c[i],Amp0_c[i],
            Amp1_c[i], Amp2_c[i]]], axis = 0)
            datos_s = np.append(datos_s,[[snapshots_analysis[t],Rcenters_s[i],Amp0_s[i],
            Amp1_s[i], Amp2_s[i]]], axis = 0)
    return datos_c, datos_s

def acceleration_fourier (mode=1, plot_means=False, nbins = 22, maximo=22, minimo=0):
    stepr=(maximo-minimo)/nbins
    binsr=np.arange(minimo,maximo+stepr,stepr)
    Rcenters = (binsr[:-1] + binsr[1:]) / 2
    data =pd.read_csv(path_datos + "acceleration_mean_comparison_inter.csv", sep = ",", index_col = 0)
    modos = ["dm", "gas", "stars_disk", "stars_nodisk"]

    modos_gal = ["dm", "gas", "disk_stars", "nodisk_stars"]
    modos_sat = ["progenitor", "stream"]
    modos_completos = modos_gal + modos_sat

    plt.figure()
    fig, ax = plt.subplots(nrows=6, figsize=(12,18))
    fig.subplots_adjust(hspace=0.38)
    for index, modo in enumerate(modos):
        print(modo)
        FAMP= apply_fourier(modo)
        AMP_0 = FAMP[:,2].reshape(len(snapshots_analysis), nbins)
        AMP_1 = FAMP[:,3].reshape(len(snapshots_analysis), nbins)
        ax[index].pcolormesh(lookback,Rcenters, (AMP_1/AMP_0).T, 
                            vmin =  5e-17, vmax = 5e-15)
    FAMP_c, FAMP_s= apply_fourier_sat()
    AMP0_c = FAMP_c[:,2].reshape(len(snapshots_analysis), nbins)
    AMP0_s = FAMP_s[:,2].reshape(len(snapshots_analysis), nbins)
    if modo ==1:
        AMP1_c = FAMP_c[:,3].reshape(len(snapshots_analysis), nbins)
        AMP1_s = FAMP_s[:,3].reshape(len(snapshots_analysis), nbins)
    if modo ==2:
        AMP1_c = FAMP_c[:,4].reshape(len(snapshots_analysis), nbins)
        AMP1_s = FAMP_s[:,4].reshape(len(snapshots_analysis), nbins)

    ax[4].pcolormesh(lookback,Rcenters,(AMP1_c/AMP0_c).T, 
                    vmin = 3e-19, vmax = 5e-17)
    ax[5].pcolormesh(lookback,Rcenters,(AMP1_s/AMP0_s).T,
                    vmin = 3e-19, vmax = 5e-17)
    for index, aa in enumerate(ax):
        ax[index].set_xlim(6.3,0.1)
        ax[index].set_ylim(2,20)
        ax[index].set_title(modos_completos[index])
    if plot_means ==True:
        for index, modo in enumerate(modos_completos):
            ax2 = ax[index].twinx()
            ax2.plot( data["Lookback"],np.log10(data["az_%s_mean"%modo]),':',marker='o', 
                        markerfacecolor="white",
                    markeredgecolor='black', markeredgewidth=1, 
                        linewidth=2, ms = 6, color = 'cyan', label =modo)
            ax2.legend()
    ax[5].set_ylabel('Radius [kpc]',fontsize = 12)
    ax[5].set_xlabel('Lookback time [Gyr]',fontsize = 12)

    plt.savefig("Figuras/fourier_acceleration_components.png", dpi = 100, facecolor = "white",bbox_inches='tight')
    #plt.show()

def bending_fourier (mode=1, nbins = 22, maximo=22, minimo=0):
    for t,name in enumerate(snapshots_analysis):
        print(name)
        dfB = pd.read_csv(path_csv + "%s_stars_Rvir.csv"%name ,sep = ",")
        disco = pd.read_csv(path_disk + "Stars_disco_%s.csv" %name)
        dfA = dfB[dfB['ID'].isin(disco["ID"])]
        #dfC = dfA[(dfA['Age']>age1 )& (dfA['Age']< age2)].copy()
    # df=  dfC[(dfC['Z']>-z_lim)& (dfC['Z']< z_lim) &(dfC['R']< 25)].copy()
        df = dfA[(dfA['R']< 25)].copy()
        df["Phi"] = np.mod(np.arctan2(df["Y"], df["X"]), 2*np.pi)
        df["R"] = np.sqrt(df["X"]**2 + df["Y"]**2)

        df.loc[df['Phi']>np.pi,'Phi']=df.loc[df['Phi']>np.pi,'Phi']-2.*np.pi
        
    #  FAMP[t,:,:],FPHA[t,:,:],Rcenters=fourier(df["R"],df["Phi"],nbins=nbins,maxmod=maxmod,
    #                                           rmax=rmax,rmin=rmin,weight=df["Z"])
        Rcenters, nparticles, Amp0, Amp1, Amp2 =fourier(df["X"],df["Y"],df["Z"])
        
        
        for i in range(nbins):
            datos = np.append(datos,[[snapshots_analysis[t],lookback[t],Rcenters[i],nparticles[i],Amp0[i],
            Amp1[i], Amp2[i]]], axis = 0)

    column_names = ['snapshot_t','lookbacktime','Rcenters','Nparticles','amp0','amp1','amp2']
    datos_df = pd.DataFrame(datos, columns = column_names)

    datos_df.to_csv('fourier_z_11bins_all.csv', sep = ',', index = False)