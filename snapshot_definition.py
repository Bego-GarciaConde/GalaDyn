
from os import path
from config import *
import yt
from yt import YTArray

import numpy as np
from yt.units import G
import array
import pandas as pd
import matplotlib.pylab as plt

from sklearn.linear_model import LinearRegression
from mpl_toolkits.axes_grid1 import make_axes_locatable
#global datos_edades
datos_edades = pd.read_csv(path_datos + "edades.csv", sep = ",",index_col = 0)
#global angular_momentum_ref
#angular_momentum_ref = YTArray([2.63425158e+29, 4.18786719e+28, -1.09858375e+29], "cm**2/s")

def cartesian_to_cylindrical (df):
    df["Phi"] =  np.mod(np.arctan2(df["Y"],df["X"]), 2*np.pi)
    df["R"] = np.sqrt(df["X"]**2 + df["Y"]**2)
    df["Vphi"] = (df["X"]*df["VY"] - df["Y"]*df["VX"])/df["R"] #todo revisar signos de phi y vphi
    df["VR"] = (df["X"]*df["VX"] + df["Y"]*df["VY"])/df["R"]
    return df


class Snapshot:
    def __init__(self, name):
        self.name = name
        self.path_snapshot = None
        self.lb = None
        self.center = None
        self.Rvir = None
        self.ds = None 
        self.dm = None
        self.gas = None
        self.stars = None
        self.disk = None
        self.disk_filt = None
        self.bending_breathing_mode = None
      


        def read_lb():
            self.lb = datos_edades.loc[datos_edades['Snapshot'] == self.name, 'Lookback'].iloc[0]

        def read_center_Rvir ():
            centro = np.loadtxt(path_datos +f'center_{self.name}.txt')
            center = YTArray([centro[0], centro[1], centro[2]], "cm")
        #   cx,cy,cz = center[0].in_units("cm"), center[1].in_units("cm"),  center[2].in_units("cm")
            Rvir = YTArray(centro[3], "kpc")
            self.center = center
            self.Rvir = Rvir

        def find_path_for_yt():
            # name = snapshots_analysis[i]
            if self.name < 425:
                path_snapshot = "/media/temp1/bego/GARROTXA_ART/"
            elif (self.name >= 425)&(name < 600):
                path_snapshot = "/srv/cab1/garrotxa/GARROTXA_ART/MW_003/RUN2.2/"
            elif (self.name >=600 )&(name < 800):
                path_snapshot = "/home/Garrotxa_ART/New_Run/"
            elif (self.name >= 800) & (name < 900) :
                path_snapshot = "/media/temp/bego/New_Resim/"
            elif self.name >= 900 :
                path_snapshot = "/media/temp1/GARROTXA_ART/MW_003/RUN2.2/"
            self.path_snapshot = path_snapshot
        
        print(f"Initializing snapshot {name}")
        find_path_for_yt()
        read_lb()
        print(f"Lookback time: {self.lb} Gyr")
        read_center_Rvir()
        
    def load_stars (self):
        self.stars = pd.read_csv(path_csv + f"{self.name}_stars_Rvir.csv",sep = ",")
        self.stars = cartesian_to_cylindrical(self.stars)

    def load_dm (self):
        self.dm = pd.read_csv(path_csv + f"{self.name}_dm_Rvir.csv",sep = ",")
        self.dm = cartesian_to_cylindrical(self.dm)

    def load_gas (self):
        self.gas = pd.read_csv(path_csv + f"Gas_{self.name}.csv",sep = ",")
        self.gas = cartesian_to_cylindrical(self.gas)

    def load_disk(self):
        self.disk = pd.read_csv(path_disk + f"Stars_disco_{self.name}.csv")

    def filter_disk_particles(self):
        dfA = self.stars[self.stars['ID'].isin(self.disk["ID"])]
        df = dfA[(dfA['R']< 25) &(dfA['Z']< 2.5)&(dfA['Z']>-2.5)].copy()
        df["Phi"] = np.mod(np.arctan2(df["Y"], df["X"]), 2*np.pi)
        df["R"] = np.sqrt(df["X"]**2 + df["Y"]**2)
        self.disk_filt = df
        return df

    def calculate_bending_breathing (self):
        xbins = np.linspace(-25, 25, 70)
        indx = np.digitize(self.disk_filt["X"], bins = xbins)
        binx_means = [self.disk_filt["X"][indx == i].mean() for i in range(1, len(xbins))]

        ybins = np.linspace(-25, 25, 70)
        indy = np.digitize(self.disk_filt["Y"], bins = ybins)
        biny_means = [self.disk_filt["Y"][indy == i].mean() for i in range(1, len(ybins))]
        pos_x = []
        pos_y = []
        A = []
        B = []

        Z = np.array(self.disk_filt["Z"])
        VZ = np.array(self.disk_filt["VZ"])

        X = np.array(self.disk_filt["X"])
        Y = np.array(self.disk_filt["Y"])
        for i in range(1, len(xbins)):
            for j in range(1, len(ybins)):
                
                indr = np.where((indx==i)&(indy==j))
                Z_i = Z[indr[0]]
                VZ_i = VZ[indr[0]]
                X_i = X[indr[0]]
                Y_i = Y[indr[0]]
                pos_x.append(np.mean(X_i))
                pos_y.append(np.mean(Y_i))
                if len(Z_i)<2:
                    B.append(0)
                    A.append(0)
                else:
                    Z_i = Z_i.reshape((-1, 1))
                    model = LinearRegression().fit(Z_i, VZ_i)
                    B.append(model.intercept_)
                    A.append(model.coef_)
        
        data = {"X": pos_x, "Y":pos_y,"Bending": np.array(B), "Breathing":np.array(A)}
        self.bending_breathing_mode=  pd.DataFrame(data)
        #return  self.bending_breathing_mode

    
    def load_accelerations(self):
        self.az_dm = pd.read_csv(path_acceleration + f"mesh_aceleracion_dm_{self.name}_ytRS_40.csv", sep = ",")
        self.az_gas =  pd.read_csv(path_acceleration + f"mesh_aceleracion_gas_{self.name}_ytRS_40.csv", sep = ",")
        
 