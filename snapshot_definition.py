
from os import path
from config import *
import yt
from yt import YTArray

import numpy as np
from yt.units import G
import array
import pandas as pd
import matplotlib.pylab as plt
from skspatial.objects import Point, Vector, Plane
from sklearn.linear_model import LinearRegression
from mpl_toolkits.axes_grid1 import make_axes_locatable

datos_edades = pd.read_csv(PATH_DATOS + "edades.csv", sep = ",",index_col = 0)

def cartesian_to_cylindrical (df):
    df["Phi"] =  np.mod(np.arctan2(df["Y"],df["X"]), 2*np.pi)
    df["R"] = np.sqrt(df["X"]**2 + df["Y"]**2)
    df["Vphi"] = (df["X"]*df["VY"] - df["Y"]*df["VX"])/df["R"] #todo revisar signos de phi y vphi
    df["VR"] = (df["X"]*df["VX"] + df["Y"]*df["VY"])/df["R"]
    return df

def cartesian_to_spherical (df):
    XsqPlusYsq = df["X"]**2 + df["Y"]**2
    df["Phi_sph"] =  np.arctan2(df["Y"],df["X"])
    df["R_sph"] = np.sqrt(XsqPlusYsq + df["Z"]**2) 
    df["Theta_sph"] = np.arctan2(df["Z"],np.sqrt(XsqPlusYsq))  
    df["Vr_sph"] = (df["X"]*df["VX"] + df["Y"]*df["VY"] +df["Z"]*df["VZ"])/df["R_sph"]
    return df

def centering(pos0,vel0,mass, ind,means=True,medians=False,L=True):
    pos=pos0.copy()
    vel=vel0.copy()
    
    pos2=pos0.copy()
    vel2=vel0.copy()
    
    #centering particles with means   
    if means==True:
        for j in range(3):
            mp=np.mean(pos0[ind,j])
            mv=np.mean(vel0[ind,j])
            print(j,mp,mv)
            pos[:,j]=pos0[:,j]-mp
            vel[:,j]=vel0[:,j]-mv

    #centering particles with medians   
    if medians==True:
        for j in range(3):
            pos[:,j]=pos0[:,j]-np.median(pos0[ind,j])
            vel[:,j]=vel0[:,j]-np.median(vel0[ind,j])    
    
    #reorienting particles so meadian L is oriented as Lz (Lx=Ly=0) 
    # the other free axis X is taken so as to produce the minimum variation of the old axis X
    if L==True:    
        Lx=pos[:,1]*vel[:,2]-pos[:,2]*vel[:,1]
        Ly=pos[:,2]*vel[:,0]-pos[:,0]*vel[:,2]
        Lz=pos[:,0]*vel[:,1]-pos[:,1]*vel[:,0]

       #finding median L vector only with selected particles
       # mLx,mLy,mLz= np.average(Lx[ind], weights=mass[ind]),np.average(Ly[ind], weights= mass[ind]),np.average(Lz[ind],weights= mass[ind])
        mLx,mLy,mLz= np.average(Lx[ind]),np.average(Ly[ind]),np.average(Lz[ind])
        m=np.sqrt(mLx*mLx+mLy*mLy+mLz*mLz)

        mLx,mLy,mLz=-mLx/m,-mLy/m,-mLz/m #normalization of the median L vector         
       # mLx,mLy,mLz=mLx/m,mLy/m,mLz/m
        v3=[mLx,mLy,mLz]

        plane = Plane(point=[0, 0, 0], normal=[mLx,mLy,mLz])

        #Projection of the X axis
        point = Point([1, 0, 0])
        v1 = plane.project_point(point) 
        m=np.sqrt(v1[0]**2+v1[1]**2+v1[2]**2)
        v1=v1/m

        #Third vector is the cross product of z' x'
        v2=np.cross(v3,v1)
        m=np.sqrt(v2[0]**2+v2[1]**2+v2[2]**2)
        v2=v2/m           
        A=np.matrix([np.array(v1),np.array(v2),np.array(v3)]).T

        B=np.linalg.inv(A)

    return B
 
    
def apply_transformation_matrix(matriz_trans, X,Y,Z):
    X_re= X*matriz_trans[0,0]+ Y*matriz_trans[0,1] + Z*matriz_trans[0,2]       
    Y_re= X*matriz_trans[1,0]+ Y*matriz_trans[1,1] + Z*matriz_trans[1,2]
    Z_re= X*matriz_trans[2,0]+ Y*matriz_trans[2,1] + Z*matriz_trans[2,2]
    return X_re, Y_re, Z_re

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
        self.mC = None
      


        def read_lb():
            self.lb = datos_edades.loc[datos_edades['Snapshot'] == self.name, 'Lookback'].iloc[0]

        def read_center_Rvir ():
            centro = np.loadtxt(PATH_DATOS +f'center_{self.name}.txt')
            center = YTArray([centro[0], centro[1], centro[2]], "cm")
        #   cx,cy,cz = center[0].in_units("cm"), center[1].in_units("cm"),  center[2].in_units("cm")
            Rvir = YTArray(centro[3], "kpc")
            self.center = center
            self.Rvir = Rvir

        def load_2nd_align():
             self.mC = np.loadtxt(PATH_2ND_ALIGNMENT + f"{name}_{ALIGN_RADIUS}kpc_matrix.txt")

        def find_path_for_yt():
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
        load_2nd_align()

        
    def load_stars (self):
        self.stars = pd.read_csv(path_csv + f"{self.name}_stars_Rvir.csv",sep = ",")
        if APPLY_ALIGNMENT == 1:
            self.stars = self.apply_align(self.stars)
        self.stars = cartesian_to_cylindrical(self.stars)

    def load_dm (self):
        self.dm = pd.read_csv(path_csv + f"{self.name}_dm_Rvir.csv",sep = ",")
        if APPLY_ALIGNMENT == 1:
            self.dm = self.apply_align(self.dm)
        self.dm = cartesian_to_cylindrical(self.dm)

    def load_gas (self):
        self.gas = pd.read_csv(path_csv + f"Gas_{self.name}.csv",sep = ",")
        if APPLY_ALIGNMENT == 1:
            self.gas = self.apply_align(self.gas)
        self.gas = cartesian_to_cylindrical(self.gas)

    def load_disk(self):
        #self.disk
        a = pd.read_csv(path_disk + f"cla_disco_{self.name}.csv")
        self.disk = a[a["cos_alpha"]>0.7]

    def apply_align(self, df):
        df["X"],df["Y"],df["Z"] = apply_transformation_matrix(self.mC,df["X"] ,df["Y"],df["Z"])
        df["VX"],df["VY"],df["VZ"] = apply_transformation_matrix(self.mC,df["VX"] ,df["VY"],df["VZ"])
        return df


    def filter_disk_particles(self):
        dfA = self.stars[self.stars['ID'].isin(self.disk["ID"])]
        age_since_merger = 9000 - (self.lb - lb_ref)*1000
        print("Age of stars since merger")
        dfA = dfA[dfA['Age']< age_since_merger]
        df = dfA[(dfA['R']< 25) &(dfA['Z']< 2.5)&(dfA['Z']>-2.5)].copy()
        self.disk_filt = df
        return df

    def filter_disk_particles_by_age(self):
        age_since_merger = 9000 - (self.lb - lb_ref)*1000
        print("Age of stars since merger")
        dfA = self.stars[self.stars['Age']< age_since_merger]
        df = dfA[(dfA['R']< 25) &(dfA['Z']< 2.5)&(dfA['Z']>-2.5)].copy()
        self.disk_filt = df
        return df

    def filter_young(self):
        dfA = self.stars[self.stars['Age']< 2000]
        df = dfA[(dfA['R']< 25) &(dfA['Z']< 2.5)&(dfA['Z']>-2.5)].copy()
        self.disk_filt = df
        return df

    def filter_intermediate(self):
        print("Age of stars since merger")
        dfA = self.stars[(self.stars['Age']> 2000)&(self.stars['Age']< 5000)]
        df = dfA[(dfA['R']< 25) &(dfA['Z']< 2.5)&(dfA['Z']>-2.5)].copy()
        self.disk_filt = df
        return df

    
    def filter_stellar_ellipsoid(self):
        age_since_merger = 9000 - (self.lb - lb_ref)*1000
        print("Age of stars since merger")
        dfA = self.stars[self.stars['Age']> age_since_merger]
        a = pd.read_csv(path_disk + f"cla_disco_{self.name}.csv")
        a = a[a["cos_alpha"]<0.7]
        dfA = dfA[dfA['ID'].isin(a["ID"])]
        return dfA
        
    def calculate_second_alignment(self):
        """
        The second alignment further aligns the galaxy with a selected R
        15 kpc --> better aligns intermediate regions (phase spirals)
        7 kpc  --> better aligns inner regions (bending modes)

        """

        print("Applying third alignment")
        df = self.filter_disk_particles_by_age()
        ind=(df['R']<ALIGN_RADIUS)&(df['Z']<3)&(df["Z"]>-3)
        pos0=df[['X','Y','Z']].to_numpy()
        vel0=df[['VX','VY','VZ']].to_numpy()
     
        np.savetxt(PATH_2ND_ALIGNMENT + f'{name}_{align_radius}kpcmC.txt', self.C)

      #  self.load_second_align()
            
    
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
                if len(Z_i)<10:
                    B.append(0)
                    A.append(0)
                else:
                    Z_i = Z_i.reshape((-1, 1))
                    model = LinearRegression().fit(Z_i, VZ_i)
                    B.append(model.intercept_)
                    A.append(model.coef_)
        
        data = {"X": pos_x, "Y":pos_y,"Bending": np.array(B), "Breathing":np.array(A)}
        self.bending_breathing_mode=  pd.DataFrame(data)


    
    def load_accelerations(self):
        self.az_dm = pd.read_csv(PATH_ACCELERATION + f"mesh_aceleracion_dm_{self.name}_ytRS_40.csv", sep = ",")
        self.az_gas =  pd.read_csv(PATH_ACCELERATION + f"mesh_aceleracion_gas_{self.name}_ytRS_40.csv", sep = ",")
        
 