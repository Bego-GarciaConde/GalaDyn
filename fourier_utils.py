from dataclasses import dataclass
import numpy as np
import pandas as pd
import gc
from multiprocessing import Pool
from config import *
import matplotlib.pyplot as plt
from snapshot_definition import Snapshot

class Fourier:
    def __init__(self, snapshots_analysis, lookback, maximo = 22, minimo = 0,nbins = 22, maxmode = 2):
        self.maximo = maximo
        self.minimo = minimo
        self.nbins = nbins
        self.maxmode = maxmode
        self.snapshots_analysis = snapshots_analysis
        self.lookback = lookback

#-----------------FOURIER -----------------------------------------------

    def fourier_method(self, X,Y, peso = None):
        """  This function divides the disk in radial and operates fourier modes from 0 to maxmode 
            in galactic disks, minimo and maximo measured in kpc  """

        # Matrix of values of amplitudes for each mode and radial bin
        AA =  np.zeros((self.maxmode+1,self.nbins))
        armangle =  np.zeros((self.maxmode+1,self.nbins))
        steparm = (self.maximo-self.minimo)/self.nbins
        radarm = np.arange(self.minimo,self.maximo + steparm, steparm )
        rcenter = (radarm[:-1] + radarm[1:]) / 2
        A = np.zeros(self.maxmode+1, dtype=np.float32)
        B = np.zeros(self.maxmode+1, dtype=np.float32)

        dd = np.sqrt(X**2 + Y**2)
        indr=np.digitize(dd,radarm)-1
        nparticles = np.zeros(self.nbins)
        #Iterating over fourier modes
        for m in range(0,self.maxmode+1):
            #Iterating over radial bins
            for i in range(1,self.nbins):
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
                        A[m] = np.sum(np.cos(m*a))
                        B[m] = np.sum(np.sin(m*a))
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


        return rcenter, nparticles, AA, armangle





    #------------------accelerations-------------------
    
    def apply_fourier_accelerations(self, comp, nbins):
        datos = np.zeros((len(snapshots_analysis)*self.nbins, 6 + 2*self.maxmode))
        #print(np.shape(datos))
        index = 0
        for t,name in enumerate(snapshots_analysis):
            print(name)
            df = pd.read_csv(path_acceleration + f"mesh_aceleracion_{comp}_{name}_ytRS_{nbins}.csv",sep = ",")
            limite = 1e-14
            for j,la in enumerate(df["az"]):
                if df["az"][j]>limite:
                    df["az"][j]=limite
                elif df["az"][j]< -limite:
                    df["az"][j]=limite
            Rcenters, npart, amplitudes, phases =self.fourier_method(df["X"],df["Y"],peso=df["az"])
            for i in range(self.nbins):
             #   print([snapshots_analysis[t],lookback[t], Rcenters[i],npart[i]] + list(amplitudes[:,i]) + list(phases[:,i]))
                datos[index] = [snapshots_analysis[t],lookback[t],Rcenters[i],npart[i]] + list(amplitudes[:,i]) + list(phases[:,i])
                index = index + 1
        self.save_fourierogram(datos, etiqueta=f"acceleration_{comp}", peso = "az")


    #------------------satellites-------------------
    
    def apply_fourier_sat(self):
        datos_c = np.zeros((len(snapshots_analysis)*self.nbins, 6 + 2*self.maxmode))
        datos_s = np.zeros((len(snapshots_analysis)*self.nbins, 6 + 2*self.maxmode))

        index = 0
        for t,name in enumerate(snapshots_analysis):
            df = pd.read_csv(path_acceleration + f"mesh_aceleracion_{name}_all_satellites_id_ytRS.csv" ,sep = ",")
            Rcenters_s, npart_s, Amp_s, phases_s = self.fourier_method(df["X"],df["Y"],peso=df["az_stream"])
            Rcenters_c, npart_c, Amp_c, phases_c = self.fourier_method(df["X"],df["Y"],peso=df["az_core"])
            for i in range(self.nbins):
                datos_c[index] = [snapshots_analysis[t],lookback[t],Rcenters_c[i],npart_s[i]]+ list(Amp_c[:,i]) + list(phases_c[:,i]) 
                datos_s[index] = [snapshots_analysis[t],lookback[t],Rcenters_s[i],npart_c[i]]+ list(Amp_s[:,i]) +  list(phases_s[:,i])
            #  datos_c = np.append(datos_c,[[snapshots_analysis[t],Rcenters_c[i],npart_s[i], Amp_c[:,i]]], axis = 0)
            #  datos_s = np.append(datos_s,[[snapshots_analysis[t],Rcenters_s[i],npart_c[i], Amp_s[:,i]]], axis = 0)
                index = index +1

        self.save_fourierogram(datos_c, etiqueta=f"sat_prog", peso = "az_core")
        self.save_fourierogram(datos_s, etiqueta=f"sat_streams", peso = "az_stream")
       # return datos_c, datos_s


    #------------------disk-------------------

    def apply_fourier_on_disk (self, peso=None):
        print(f"Analyzing fourierograms of {peso}")
        #Initializing data 
        datos = np.zeros((len(snapshots_analysis)*self.nbins, 6 +2*self.maxmode), dtype = np.float32)
    
        index = 0
        #Iterating over snapshots
        for t,name in enumerate(snapshots_analysis):
            print(name)
            etiqueta = "disc_5Gyr"
            snapshot = Snapshot(name)
            snapshot.load_stars()
            snapshot.load_disk()
           # df = snapshot.stars[(snapshot.stars["R"]< 25)&(snapshot.stars["Z"]< 2.7)&(snapshot.stars["Z"]> -2.7) &(snapshot.stars["Age"]<5000)]
           # dfA = snapshot.filter_disk_particles()
            #df = dfA[(dfA["R"]< 25)&(dfA["Z"]< 2.7)&(dfA["Z"]> -2.7) &(dfA["Age"]<5000)]
            df = snapshot.stars[(snapshot.stars["Age"]<5000)]
            print("Snapshot loaded!")

            #Apply fourier
            if peso == None:
                Rcenters, nparticles, amplitudes, phases = self.fourier_method(X= df["X"],Y = df["Y"])
            else:
                Rcenters, nparticles, amplitudes, phases = self.fourier_method(X = df["X"],Y= df["Y"], peso = df[f"{peso}"])
            for i in range(self.nbins):
                datos[index] = [snapshots_analysis[t], lookback[t], Rcenters[i], nparticles[i]] + list(amplitudes[:,i]) + list(phases[:,i])
                index = index +1

        self.save_fourierogram(datos, etiqueta=etiqueta, peso = peso)
        

    def apply_fourier_on_bar (self, peso=None):
        print(f"Analyzing fourierograms of {peso}")
        #Initializing data 
        datos = np.zeros((len(snapshots_analysis)*self.nbins, 6 +2*self.maxmode), dtype = np.float32)
        etiqueta = "bar_35_percentile_Gyr"
        index = 0
        #Iterating over snapshots
        for t,name in enumerate(snapshots_analysis):
            print(name)
            snapshot = Snapshot(name)
            snapshot.load_stars()
            snapshot.load_disk()
            df = snapshot.stars[~snapshot.stars['ID'].isin(snapshot.disk["ID"])]
            edad_a = np.percentile(df["Age"], 35)
            df = df[(df['Age']< edad_a)].copy()
            df = df[(df["R"]< 7)&(df["Z"]< 5)&(df["Z"]> -5)]
           # df = snapshot.stars[(snapshot.stars["R"]< 25)&(snapshot.stars["Z"]< 2.7)&(snapshot.stars["Z"]> -2.7) &(snapshot.stars["Age"]<5000)]
           # dfA = snapshot.filter_disk_particles()
          #  df = dfA[(dfA["R"]< 25)&(dfA["Z"]< 2.7)&(dfA["Z"]> -2.7) &(dfA["Age"]<5000)]
           # df = snapshot.stars[(snapshot.stars["Age"]<5000)]
            print("Bar loaded")
            #print(f"The bar has {}")

            #Apply fourier
            if peso == None:
                Rcenters, nparticles, amplitudes, phases = self.fourier_method(X= df["X"],Y = df["Y"])
            else:
                Rcenters, nparticles, amplitudes, phases = self.fourier_method(X = df["X"],Y= df["Y"], peso = df[f"{peso}"])
            for i in range(self.nbins):
                datos[index] = [snapshots_analysis[t], lookback[t], Rcenters[i], nparticles[i]] + list(amplitudes[:,i]) + list(phases[:,i])
                index = index +1

        self.save_fourierogram(datos, etiqueta=etiqueta, peso = peso)

        
    def save_fourierogram(self, datos, etiqueta, peso):

        column_names = ['snapshot_t','lookbacktime','Rcenters','Nparticles']
        for numero_modo in range(0, self.maxmode +1):
            column_names.append(f"amp{numero_modo}")
        for numero_modo in range(0, self.maxmode +1):
            column_names.append(f"phase{numero_modo}")

        
        datos = pd.DataFrame(datos, columns=column_names)
        
        print(f"Saving data as {path_results} fourier_{self.nbins}_{peso}_{etiqueta}.csv")
        datos.to_csv(path_results + f'fourier_{self.nbins}_{peso}_{etiqueta}.csv', sep = ',', index = False)