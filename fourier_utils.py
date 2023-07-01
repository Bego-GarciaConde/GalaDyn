from dataclasses import dataclass
import numpy as np
import pandas as pd
import gc
from multiprocessing import Pool
from config import *
import matplotlib.pyplot as plt
from snapshot_definition import Snapshot
from snapshot_definition import cartesian_to_spherical

class Fourier:

    def __init__(self, snapshots_analysis, lookback, maximo = 22, minimo = 0,nbins = 22, maxmode = 4):
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
                a = np.arctan2(Y_i, X_i)

                if peso is None:
                    factor = 1 if m == 0 else 2
                    A[m] = factor * np.sum(np.cos(m * a))
                    B[m] = factor * np.sum(np.sin(m * a))
                else:
                    peso_i = peso[indr == i]
                    factor = 1 if m == 0 else 2
                    A[m] = factor * np.sum(peso_i * np.cos(m * a))
                    B[m] = factor * np.sum(peso_i * np.sin(m * a))

                
                AA[m,i] = np.sqrt(A[m]**2+ B[m]**2) 
                armangle[m,i] = np.arctan2(B[m], A[m]) if m > 0 else 0
                nparticles[i] = len(a) if m == 0 else nparticles[i]
                    

        return rcenter, nparticles, AA, armangle


    #------------------accelerations-------------------
    
    def apply_fourier_accelerations(self, comp):
        datos = np.zeros((len(snapshots_analysis)*self.nbins, 6 + 2*self.maxmode))
        index = 0
        print("Componmente: ", comp)
        for t,name in enumerate(snapshots_analysis):
            print(name)
            if comp == "dm" or comp == "gas":
                print(comp)
                df = pd.read_csv(path_acceleration + f"mesh_aceleracion_{comp}_{name}_ytRS_{self.nbins}.csv",sep = ",")
            elif comp == "sum":
              df_dm = pd.read_csv(path_acceleration + f"mesh_aceleracion_dm_{name}_ytRS_{self.nbins}.csv",sep = ",")
              df_gas = pd.read_csv(path_acceleration + f"mesh_aceleracion_gas_{name}_ytRS_{self.nbins}.csv",sep = ",")
              df = df_dm.copy()
              df["az"] = np.array(df_dm["az"]) + np.array(df_gas["az"])
            elif comp == "stars_ellipsoid":
                df = pd.read_csv(path_acceleration + f"mesh_aceleracion_{comp}_{name}_ytRS_{self.nbins}.csv",sep = ",")
                limite = 4.8e-14
              
        #    limite = 4.8e-14 #This limit is calculated with MG/z**2, being M the maximum mass a 
            #star particle can have an z the minimum resolution of the simulation (108 pc, 54pc in absolute value)
            if comp == "gas":
                limite = 4.8e-14
            elif comp == "dm" or comp == "sum":
              limite = 4.8e-15

            for j, la in enumerate(df["az"]):
                if df["az"][j]>limite:
                    df["az"][j]=limite
                elif df["az"][j]< -limite:
                    df["az"][j]=-limite

            if torque==1:
                Rcenters, npart, amplitudes, phases =self.fourier_method(df["X"],df["Y"],peso=df["az"]*np.sqrt(df["X"]**2 +df["Y"]**2 ))
                etiqueta = f"acceleration_{comp}"
            else:
                Rcenters, npart, amplitudes, phases =self.fourier_method(df["X"],df["Y"],peso=df["az"])
                etiqueta = f"acceleration_{comp}_noise"

            for i in range(self.nbins):
                datos[index] = [snapshots_analysis[t],lookback[t],Rcenters[i],npart[i]] + list(amplitudes[:,i]) + list(phases[:,i])
                index = index + 1
        self.save_fourierogram(datos, etiqueta=etiqueta, peso = "az", subfolder = "accelerations")



    #------------------satellites-------------------
    
    def apply_fourier_sat(self, sat_name):
        datos_c = np.zeros((len(snapshots_analysis)*self.nbins, 6 + 2*self.maxmode))
        datos_s = np.zeros((len(snapshots_analysis)*self.nbins, 6 + 2*self.maxmode))

        index = 0
        for t,name in enumerate(snapshots_analysis):
            print(name)
            df = pd.read_csv(path_acceleration + f"mesh_aceleracion_{name}_{sat_name}_satellites_id_ytRS.csv" ,sep = ",")
            Rcenters_s, npart_s, Amp_s, phases_s = self.fourier_method(df["X"],df["Y"],peso=df["az_stream"])
            Rcenters_c, npart_c, Amp_c, phases_c = self.fourier_method(df["X"],df["Y"],peso=df["az_core"])
            for i in range(self.nbins):
                datos_c[index] = [snapshots_analysis[t],lookback[t],Rcenters_c[i],npart_s[i]]+ list(Amp_c[:,i]) + list(phases_c[:,i]) 
                datos_s[index] = [snapshots_analysis[t],lookback[t],Rcenters_s[i],npart_c[i]]+ list(Amp_s[:,i]) +  list(phases_s[:,i])
  
                index = index +1

        self.save_fourierogram(datos_c, etiqueta=f"sat_{sat_name}", peso = "az_core", subfolder = "satellites")
        self.save_fourierogram(datos_s, etiqueta=f"sat_{sat_name}", peso = "az_stream", subfolder = "satellites")
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
            etiqueta = f"disc_{tag}"
            snapshot = Snapshot(name)
            snapshot.load_stars()
            snapshot.load_disk()
            df =  filter_disk_particles_by_age()
          #  df = snapshot.filter_disk_particles()


            print("Snapshot loaded!")

            #Apply fourier
            if peso == None:
                Rcenters, nparticles, amplitudes, phases = self.fourier_method(X= df["X"],Y = df["Y"])
            else:
                Rcenters, nparticles, amplitudes, phases = self.fourier_method(X = df["X"],Y= df["Y"], peso = df[f"{peso}"])
            for i in range(self.nbins):
                datos[index] = [snapshots_analysis[t], lookback[t], Rcenters[i], nparticles[i]] + list(amplitudes[:,i]) + list(phases[:,i])
                index = index +1

        self.save_fourierogram(datos, etiqueta=etiqueta, peso = peso, subfolder = "disc")
        

    def apply_fourier_on_bar (self, stars_or_dm = "stars", peso=None):
        print(f"Analyzing fourierograms of {peso}")
        #Initializing data 
        datos = np.zeros((len(snapshots_analysis)*self.nbins, 6 +2*self.maxmode), dtype = np.float32)
        etiqueta = f"{stars_or_dm}_inner_structure_20kpc"
        index = 0
        #Iterating over snapshots
        for t,name in enumerate(snapshots_analysis):
            print(name)
            snapshot = Snapshot(name)
            if stars_or_dm == "stars":
                snapshot.load_stars()
                snapshot.load_disk()
                df = snapshot.stars[~snapshot.stars['ID'].isin(snapshot.disk["ID"])]
                edad_a = np.percentile(df["Age"], 35)
                df = df[(df['Age']< edad_a)].copy()
                df = df[(df["R"]< 7)&(df["Z"]< 5)&(df["Z"]> -5)]
            if stars_or_dm == "dm":
                snapshot.load_dm()
                df = cartesian_to_spherical(snapshot.dm)
               # df = df[df["R_sph"]<8].copy()
                limit_AB =20
                df = df[(df["R_sph"]<limit_AB)&(np.abs(df["Z"])<10)].copy()

            #Apply fourier
            if peso == None:
                Rcenters, nparticles, amplitudes, phases = self.fourier_method(X= df["X"],Y = df["Y"])
            else:
                Rcenters, nparticles, amplitudes, phases = self.fourier_method(X = df["X"],Y= df["Y"], peso = df[f"{peso}"])
            for i in range(self.nbins):
                datos[index] = [snapshots_analysis[t], lookback[t], Rcenters[i], nparticles[i]] + list(amplitudes[:,i]) + list(phases[:,i])
                index = index +1

        self.save_fourierogram(datos, etiqueta=etiqueta, peso = peso, subfolder = "bar")



    def apply_fourier_on_bending_breathing(self):
        print(f"Analyzing fourierograms of bending and breathing")
        #Initializing data 
        datos_bending  = np.zeros((len(snapshots_analysis)*self.nbins, 6 +2*self.maxmode), dtype = np.float32)
        datos_breathing  = np.zeros((len(snapshots_analysis)*self.nbins, 6 +2*self.maxmode), dtype = np.float32)
    
        index_bending = 0
        index_breathing = 0
        #Iterating over snapshots
        for t,name in enumerate(snapshots_analysis):
            print(name)

            snapshot = Snapshot(name)
            snapshot.load_stars()
            print("Snapshot loaded!")
            snapshot.calculate_bending_breathing()
            df = snapshot.bending_breathing_mode

            #Apply fourier

            Rcenters, nparticles, amplitudes, phases = self.fourier_method(X = df["X"],Y= df["Y"], peso = df["Bending"])
            for i in range(self.nbins):
                datos_bending[index_bending ] = [snapshots_analysis[t], lookback[t], Rcenters[i], nparticles[i]] + list(amplitudes[:,i]) + list(phases[:,i])
                index_bending  = index_bending  +1
            Rcenters, nparticles, amplitudes, phases = self.fourier_method(X = df["X"],Y= df["Y"], peso = df["Breathing"])
            for i in range(self.nbins):
                datos_breathing[index_breathing ] = [snapshots_analysis[t], lookback[t], Rcenters[i], nparticles[i]] + list(amplitudes[:,i]) + list(phases[:,i])
                index_breathing  = index_breathing  +1


        self.save_fourierogram(datos_bending, etiqueta=tag, peso = "bending", subfolder = "disc")
        self.save_fourierogram(datos_breathing, etiqueta=tag, peso = "breathing", subfolder = "disc")

        
    def save_fourierogram(self, datos, etiqueta, peso, subfolder):

        column_names = ['snapshot_t','lookbacktime','Rcenters','Nparticles']
        for numero_modo in range(0, self.maxmode +1):
            column_names.append(f"amp{numero_modo}")
        for numero_modo in range(0, self.maxmode +1):
            column_names.append(f"phase{numero_modo}")

        
        datos = pd.DataFrame(datos, columns=column_names)
        
        print(f"Saving data as {path_results} fourier_{self.nbins}_{peso}_{etiqueta}.csv")
        datos.to_csv(path_results + f'{subfolder}/fourier_{self.nbins}_{peso}_{etiqueta}.csv', sep = ',', index = False)