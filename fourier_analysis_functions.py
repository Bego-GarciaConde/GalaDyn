from dataclasses import dataclass
import numpy as np
import pandas as pd
import gc
from multiprocessing import Pool
from config import *
import matplotlib.pyplot as plt


#-----------------FOURIER -----------------------------------------------

def fourier(X,Y, peso=None, maximo = 22, minimo = 0,nbins = 22, maxmode = 2):
    """  This function divides the disk in radial and operates fourier modes from 0 to maxmode 
        in galactic disks, minimo and maximo measured in kpc  """

    # Matrix of values of amplitudes for each mode and radial bin
    AA =  np.zeros((maxmode+1,nbins), dtype=np.float32)
    armangle =  np.zeros((maxmode+1,nbins))
    steparm = (maximo-minimo)/nbins
    radarm = np.arange(minimo,maximo + steparm, steparm )
    rcenter = (radarm[:-1] + radarm[1:]) / 2
    A = np.zeros(maxmode+1, dtype=np.float32)
    B = np.zeros(maxmode+1, dtype=np.float32)

    dd = np.sqrt(X**2 + Y**2)
    indr=np.digitize(dd,radarm)-1
    nparticles = np.zeros(nbins)
    #Iterating over fourier modes
    for m in range(0,maxmode+1):
       #Iterating over radial bins
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


    return rcenter, nparticles, AA


def apply_fourier_accelerations(modo, nbins = 22, maxmode = 2):
    datos = np.zeros((len(snapshots_analysis)*nbins, 5 +maxmode))
    index = 0
    for t,name in enumerate(snapshots_analysis):
        df = pd.read_csv(path_acceleration + "mesh_aceleracion_%s_%s_ytRS.csv"%(modo, name) ,sep = ",")
        Rcenters, npart, amplitudes =fourier(df["X"],df["Y"],peso=df["az"], maxmode = 2)
        for i in range(nbins):
            datos[index] = [snapshots_analysis[t],Rcenters[i],npart[i]]+ list(amplitudes[:,i])
            index = index + 1
    return datos

def apply_fourier_sat(nbins = 22, maxmode=2):
    datos_c = np.zeros((len(snapshots_analysis)*nbins, 5 +maxmode))
    datos_s = np.zeros((len(snapshots_analysis)*nbins, 5 +maxmode))

    index = 0
    for t,name in enumerate(snapshots_analysis):
        df = pd.read_csv(path_acceleration + "mesh_aceleracion_%s_all_satellites_id_ytRS.csv"%(name) ,sep = ",")
        Rcenters_s, npart_s, Amp_s =fourier(df["X"],df["Y"],peso=df["az_stream"])
        Rcenters_c, npart_c, Amp_c=fourier(df["X"],df["Y"],peso=df["az_core"])
        for i in range(nbins):
            datos_c[index] = [snapshots_analysis[t],lookback[t],Rcenters_c[i],npart_s[i]]+ list(Amp_c[:,i])
            datos_s[index] = [snapshots_analysis[t],lookback[t],Rcenters_s[i],npart_c[i]]+ list(Amp_s[:,i])
          #  datos_c = np.append(datos_c,[[snapshots_analysis[t],Rcenters_c[i],npart_s[i], Amp_c[:,i]]], axis = 0)
          #  datos_s = np.append(datos_s,[[snapshots_analysis[t],Rcenters_s[i],npart_c[i], Amp_s[:,i]]], axis = 0)
            index = index +1
    return datos_c, datos_s

def fourier_analysis_accelerations (mode=1, plot_means=False, nbins = 22, maximo=22, minimo=0):
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
        FAMP= apply_fourier_accelerations(modo)
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

#------------------Z and VZ-------------------

def filter_disk_particles(name):
    dfB = pd.read_csv(path_csv + f"{name}_stars_Rvir.csv",sep = ",")
    disco = pd.read_csv(path_disk + f"Stars_disco_{name}.csv")
    dfA = dfB[dfB['ID'].isin(disco["ID"])]
    df = dfA[(dfA['R']< 25)].copy()
    df["Phi"] = np.mod(np.arctan2(df["Y"], df["X"]), 2*np.pi)
    df["R"] = np.sqrt(df["X"]**2 + df["Y"]**2)
    return df


def fourier_analysis_bending (maxmode=3, nbins=22, maximo=22, minimo=0, plot=True):
    print("Searching for bending waves!")
    print("Analyzing fourierograms of Z and VZ")
    #Initializing data 
    datos_z = np.zeros((len(snapshots_analysis)*nbins, 5 +maxmode), dtype = np.float32)
    datos_vz = np.zeros((len(snapshots_analysis)*nbins, 5 +maxmode), dtype = np.float32)
 
    index = 0
    #Iterating over snapshots
    for t,name in enumerate(snapshots_analysis):
        print(name)
        etiqueta = f"{nbins}_prueba"
        df = filter_disk_particles(name)

        #Apply fourier
        Rcenters, nparticles, amplitudes_z =fourier(df["X"],df["Y"],peso=df["Z"], maxmode=maxmode)
        Rcenters, nparticles, amplitudes_vz =fourier(df["X"],df["Y"],peso=df["VZ"], maxmode=maxmode)

        #Preparing data to save into csv, iterating over bins
        for i in range(nbins):
            datos_z[index] = [snapshots_analysis[t], lookback[t], Rcenters[i], nparticles[i]] + list(amplitudes_z[:,i])
            datos_vz[index] = [snapshots_analysis[t], lookback[t], Rcenters[i], nparticles[i]] + list(amplitudes_vz[:,i])
            index = index +1

    #Saving data
    column_names = ['snapshot_t','lookbacktime','Rcenters','Nparticles']
    for numero_modo in range(0, maxmode +1):
        column_names.append(f"amp{numero_modo}")

    datos_z_df = pd.DataFrame(datos_z, columns=column_names)
    datos_vz_df = pd.DataFrame(datos_vz, columns=column_names)

    datos_z_df.to_csv(f'fourier_z_{etiqueta}.csv', sep = ',', index = False)
    datos_vz_df.to_csv(f'fourier_vz_{etiqueta}.csv', sep = ',', index = False)

    