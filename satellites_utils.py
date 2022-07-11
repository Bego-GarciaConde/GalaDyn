import numpy as np
import pandas as pd
from config import *



#global datos_edades
#datos_edades = pd.read_csv(path_datos + "edades.csv", sep = ",",index_col = 0)


class Satellite:
    def __init__(self, name, sat):
        self.name = name
        self.sat = sat
        self.lb = None
        self.coord = None
        self.mass = None
        self.tidalR = None
        self.stars = None
        self.dm = None
        self.core = None
        self.stream = None

        def read_lb():
            self.lb = datos_edades.loc[datos_edades['Snapshot'] == self.name, 'Lookback'].iloc[0]

        def read_coord_mass_tidal():
            tabla_coord = pd.read_csv(path_datos + f"satelite_{self.sat}_coordinates.csv",
                               sep = ",")
            tabla_mass = pd.read_csv(path_datos + f"{self.sat}_complete_mass_table.csv" ,
                               sep = ",")
            x = tabla_coord.loc[tabla_coord['Snapshot'] == name, 'X'].iloc[0]
            y = tabla_coord.loc[tabla_coord['Snapshot'] == name, 'Y'].iloc[0]
            z = tabla_coord.loc[tabla_coord['Snapshot'] == name, 'Z'].iloc[0]
            self.coord = [x,y,z]
            self.mass = tabla_mass.loc[tabla_mass['Snapshot'] == name, 'Mass'].iloc[0]
            self.tidalR = tabla_mass.loc[tabla_mass['Snapshot'] == name, 'Tidal_r'].iloc[0]

        def find_crossmatch ():
            self.stars = pd.read_csv(path_crossmatch + f"{sat}_{name}_crossmatch_stars.csv", sep = ",")
            self.dm = pd.read_csv(path_crossmatch + f"{sat}_{name}_crossmatch_dm.csv", sep = ",")

        print(f"Initializing snapshot {name}")
        
        read_lb()
        
        if sat != "all":
            read_coord_mass_tidal()
            find_crossmatch()
            print(self.tidalR)
        

    def separate_stream_remanent(self):
        all_data = pd.concat ([self.dm, self.stars])

        R_sat = np.sqrt((all_data["X"] - self.coord[0])**2 + (all_data["Y"] - self.coord[1])**2 + (all_data["Z"] - self.coord[2])**2)
        all_data["R_sat"]= R_sat
        DM_core =  all_data[(all_data['R_sat']<= 5*self.tidalR)].copy()
        DM_stream =  all_data[(all_data['R_sat']> 5*self.tidalR)].copy()
        self.core = DM_core
        self.stream = DM_stream

