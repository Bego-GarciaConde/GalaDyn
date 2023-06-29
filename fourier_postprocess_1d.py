import numpy as np
import pandas as pd
from config import *

        ##################################################
        #  Post-processing of fourierograms   #
        # This code takes the regions of interest and makes 1D variables
        #  of the fourierograms
        ##################################################


density = pd.read_csv(PATH_RESULTS + f"disc/fourier_22_Mass_disc_a2_since_merger.csv", sep = ",")
z = pd.read_csv(PATH_RESULTS + f"disc/fourier_22_Z_disk_factor2.csv", sep = ",")
vz = pd.read_csv(PATH_RESULTS + f"disc/fourier_22_VZ_disk_factor2.csv", sep = ",")
A = pd.read_csv(PATH_RESULTS + f"disc/fourier_22_breathing_since_merger_.csv", sep = ",")
B = pd.read_csv(PATH_RESULTS + f"disc/fourier_22_bending_since_merger_.csv", sep = ",")

ac_gas = pd.read_csv(PATH_RESULTS + f"accelerations/fourier_40_az_acceleration_gas.csv", sep = ",")
ac_dm = pd.read_csv(PATH_RESULTS + f"accelerations/fourier_40_az_acceleration_dm.csv", sep = ",")




variables = ["Z", "VZ", "Bending","Breathing","density", "BendingM2", "BreathingM2", "densityM2", "DM", "Gas"]
variables_m1 = ["Z", "VZ", "Bending","Breathing","density", "DM", "Gas"]
variables_acc = [ "DM", "Gas"]
rings = ["10", "20"]

estimator_10 = {}
estimator_20 = {}
#generation of dictionary

estimator_10["Lookback"] = []
estimator_20["Lookback"] = []

for var1 in variables:
    estimator_10[f"{var1}"] = []
    estimator_20[f"{var1}"] = []
    
data_frames = [z ,vz, B, A,density, ac_dm, ac_gas]

def create_estimator(df, var_name, mode):
    #df_5 = df[(df["Rcenters"]<5)&(df["Rcenters"]>0)]
    df_10 = df[(df["Rcenters"]<10)&(df["Rcenters"]>2)]
    #df_15 = df[(df["Rcenters"]<15)&(df["Rcenters"]>10)]
    df_20 = df[(df["Rcenters"]<20)&(df["Rcenters"]>10)]
    for i,snapshot in enumerate(snapshots_analysis):
        df_i = df_20[df_20["snapshot_t"]==snapshot]
        estimator_20[f"{var_name}"].append(np.mean(df_i[f"amp{mode}"]/df_i["Nparticles"]))

        df_i = df_10[df_10["snapshot_t"]==snapshot]
        estimator_10[f"{var_name}"].append(np.mean(df_i[f"amp{mode}"]/df_i["Nparticles"]))

        
for i, data in enumerate(data_frames):
    create_estimator(data,variables_m1[i], 1)

create_estimator(B,"BendingM2", 2)
create_estimator(A,"BreathingM2", 2)
create_estimator(density,"densityM2", 2)

for i,snapshot in enumerate(snapshots_analysis):
    lb = datos_edades.loc[datos_edades['Snapshot'] == snapshot, 'Lookback'].iloc[0]
    estimator_10["Lookback"].append(lb)
    estimator_20["Lookback"].append(lb)

for var in variables_acc:
    estimator_10[f"{var}"] = np.array(estimator_10[f"{var}"])*seconds_to_Myr
    estimator_20[f"{var}"] = np.array(estimator_20[f"{var}"])*seconds_to_Myr

df_dm_outer = ac_dm[(ac_dm["Rcenters"]>20)]
estimator_20["DM_out"] = []

for i,snapshot in enumerate(snapshots_analysis):
    df_i = df_dm_outer[df_dm_outer["snapshot_t"]==snapshot]
    estimator_20[f"DM_out"].append(np.mean(df_i["amp1"]/df_i["Nparticles"]))

estimator_20[f"DM_out"] = np.array(estimator_20[f"DM_out"])*seconds_to_Myr

df20_all = pd.DataFrame(data=estimator_20)
df10_all = pd.DataFrame(data=estimator_10)


df20_all.to_csv(PATH_RESULTS + "results_orignal/10-20kpc_dynamic_data_since_merger.csv", sep = ",")
df10_all.to_csv(PATH_RESULTS + "results_original/0-10kpc_dynamic_data_since_merger.csv", sep = ",")







# variables = ["Z", "VZ", "Bending","Breathing","density", "BendingM2", "BreathingM2", "densityM2", "DM", "Gas"]
# variables_m1 = ["Z", "VZ", "Bending","Breathing","density", "DM", "Gas"]
# variables_acc = [ "DM", "Gas"]
# rings = ["5", "10", "15", "20"]
# estimator_5 = {}
# estimator_10 = {}
# estimator_15 = {}
# estimator_20 = {}
# #generation of dictionary
# estimator_5["Lookback"] = []
# estimator_10["Lookback"] = []
# estimator_15["Lookback"] = []
# estimator_20["Lookback"] = []
# for var1 in variables:
#     estimator_5[f"{var1}"] = []
#     estimator_10[f"{var1}"] = []
#     estimator_15[f"{var1}"] = []
#     estimator_20[f"{var1}"] = []
    
# data_frames = [z ,vz, B, A,density, ac_dm, ac_gas]

# def create_estimator(df, var_name, mode):
#     df_5 = df[(df["Rcenters"]<5)&(df["Rcenters"]>0)]
#     df_10 = df[(df["Rcenters"]<10)&(df["Rcenters"]>5)]
#     df_15 = df[(df["Rcenters"]<15)&(df["Rcenters"]>10)]
#     df_20 = df[(df["Rcenters"]<20)&(df["Rcenters"]>15)]
#     for i,snapshot in enumerate(snapshots_analysis):
#         df_i = df_20[df_20["snapshot_t"]==snapshot]
#         estimator_20[f"{var_name}"].append(np.mean(df_i[f"amp{mode}"]/df_i["Nparticles"]))
#         df_i = df_15[df_15["snapshot_t"]==snapshot]
#         estimator_15[f"{var_name}"].append(np.mean(df_i[f"amp{mode}"]/df_i["Nparticles"]))
#         df_i = df_10[df_10["snapshot_t"]==snapshot]
#         estimator_10[f"{var_name}"].append(np.mean(df_i[f"amp{mode}"]/df_i["Nparticles"]))
#         df_i = df_5[df_5["snapshot_t"]==snapshot]
#         estimator_5[f"{var_name}"].append(np.mean(df_i[f"amp{mode}"]/df_i["Nparticles"]))
        
# for i, data in enumerate(data_frames):
#     create_estimator(data,variables_m1[i], 1)

# create_estimator(B,"BendingM2", 2)
# create_estimator(A,"BreathingM2", 2)
# create_estimator(density,"densityM2", 2)

# for i,snapshot in enumerate(snapshots_analysis):
#     lb = datos_edades.loc[datos_edades['Snapshot'] == snapshot, 'Lookback'].iloc[0]
#     estimator_5["Lookback"].append(lb)
#     estimator_10["Lookback"].append(lb)
#     estimator_15["Lookback"].append(lb)
#     estimator_20["Lookback"].append(lb)
# for var in variables_acc:
#     estimator_5[f"{var}"] = np.array(estimator_5[f"{var}"])*seconds_to_Myr
#     estimator_10[f"{var}"] = np.array(estimator_10[f"{var}"])*seconds_to_Myr
#     estimator_15[f"{var}"] = np.array(estimator_15[f"{var}"])*seconds_to_Myr
#     estimator_20[f"{var}"] = np.array(estimator_20[f"{var}"])*seconds_to_Myr




# df_dm_outer = ac_dm[(ac_dm["Rcenters"]>15)]
# df_dm_inner = ac_dm[(ac_dm["Rcenters"]<5)]
# estimator_15["DM_out"] = []
# estimator_15["DM_inner"] = []
# #estimator_inner["DM"] = []
# for i,snapshot in enumerate(snapshots_analysis):

#     df_i = df_dm_outer[df_dm_outer["snapshot_t"]==snapshot]
#     estimator_15[f"DM_out"].append(np.mean(df_i["amp1"]/df_i["Nparticles"]))

#     df_i = df_dm_inner[df_dm_inner["snapshot_t"]==snapshot]
#     estimator_15[f"DM_inner"].append(np.mean(df_i["amp1"]/df_i["Nparticles"]))
    



# estimator_15[f"DM_out"] = np.array(estimator_15[f"DM_out"])*seconds_to_Myr

# estimator_15[f"DM_inner"] = np.array(estimator_15[f"DM_inner"])*seconds_to_Myr

# df20_all = pd.DataFrame(data=estimator_20)
# df15_all = pd.DataFrame(data=estimator_15)
# df10_all = pd.DataFrame(data=estimator_10)
# df5_all = pd.DataFrame(data=estimator_5)  

# df20_all.to_csv("results/15-20kpc_dynamic_data_since_merger.csv", sep = ",")
# df15_all.to_csv("results/10-15kpc_dynamic_data_since_merger.csv", sep = ",")
# df10_all.to_csv("results/5-10kpc_dynamic_data_since_merger.csv", sep = ",")
# df5_all.to_csv("results/0-5kpc_dynamic_data_since_merger.csv", sep = ",")