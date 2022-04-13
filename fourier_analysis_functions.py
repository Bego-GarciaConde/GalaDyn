import numpy as np
import pandas as pd
import gc
from multiprocessing import Pool
from config import *
import matplotlib.pyplot as plt


#-----------------FOURIER -----------------------------------------------

def fourier(X,Y, peso=None, maximo = 22, minimo = 0,nbins = 22, maxmode = 2):
    maximo = 22
    minimo =0
    #paso = 1
    #nbins = (maximo-minimo)//paso
    AA =  np.zeros((maxmode+1,nbins))
    armangle =  np.zeros((maxmode+1,nbins))
    steparm = (maximo-minimo)/nbins
    radarm = np.arange(minimo,maximo + steparm, steparm )
    #nparticles = len(X)
    #peso_final = np.mean(peso)
    A = np.zeros(maxmode+1)
    B = np.zeros(maxmode+1)

    #dd = X**2 + Y**2
    dd = np.sqrt(X**2 + Y**2)
    indr=np.digitize(dd,radarm)-1
    nparticles = np.zeros(nbins)
    for m in range(0,maxmode):
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


def apply_fourier_accelerations(modo, nbins = 22):
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
    dfB = pd.read_csv(path_csv + "%s_stars_Rvir.csv"%name ,sep = ",")
    disco = pd.read_csv(path_disk + "Stars_disco_%s.csv" %name)
    dfA = dfB[dfB['ID'].isin(disco["ID"])]
    df = dfA[(dfA['R']< 25)].copy()
    df["Phi"] = np.mod(np.arctan2(df["Y"], df["X"]), 2*np.pi)
    df["R"] = np.sqrt(df["X"]**2 + df["Y"]**2)
    return df

def fourier_analysis_bending (mode=1, nbins = 22, maximo=22, minimo=0, plot=True):
    print("Searching for bending waves!")
    print("Analyzing fourierograms of Z and VZ")
    for t,name in enumerate(snapshots_analysis):
        print(name)
        etiqueta = "%s_prueba"%nbins

        df = filter_disk_particles(name)
        Rcenters, nparticles, Amp0_z, Amp1_z, Amp2_z =fourier(df["X"],df["Y"],df["Z"], maxmode = 2)
        Rcenters, nparticles, Amp0_vz, Amp1_vz, Amp2_vz =fourier(df["X"],df["Y"],df["VZ"], maxmode = 2)
        datos_ = np.zeros((1,7))
        datos_z, datos_vz = datos_[1:],datos_[1:]

        for i in range(nbins):
            datos_z = np.append(datos_z,[[snapshots_analysis[t],lookback[t],Rcenters[i],nparticles[i],Amp0_z[i],
            Amp1_z[i], Amp2_z[i]]], axis = 0)
            datos_vz = np.append(datos_vz,[[snapshots_analysis[t],lookback[t],Rcenters[i],nparticles[i],Amp0_vz[i],
            Amp1_vz[i], Amp2_vz[i]]], axis = 0)

    column_names = ['snapshot_t','lookbacktime','Rcenters','Nparticles','amp0','amp1','amp2']
    datos_z_df = pd.DataFrame(datos_z, columns = column_names)
    datos_vz_df = pd.DataFrame(datos_vz, columns = column_names)

    datos_z_df.to_csv('fourier_z_%s.csv'%etiqueta, sep = ',', index = False)
    datos_vz_df.to_csv('fourier_vz_%s.csv'%etiqueta, sep = ',', index = False)

    