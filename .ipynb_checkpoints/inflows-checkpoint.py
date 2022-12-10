import yt
import math
from yt import YTArray

from yt.utilities.cosmology import Cosmology

co = Cosmology(hubble_constant=0.7, omega_matter=0.3, 
               omega_lambda=0.7, omega_curvature=0.0)

import numpy as np 
from yt.units import G 
import matplotlib.pyplot as plt
import os
import array


import matplotlib
from scipy.interpolate import RectBivariateSpline
import pandas as pd
import matplotlib.colors as mcolors
import os
from scipy import stats
from matplotlib.colors import LogNorm
import matplotlib.colors as colors
from os.path import expanduser
home = expanduser("~")
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.gridspec


from matplotlib import rcParams

snapshots_analysis = [520,523,525, 527,530,532,535, 537,539,541,
543, 545,547, 550, 553, 555,557, 
560, 563, 565, 567,570,573, 575, 577, 580,
583, 585,587,590, 592,594,
596,598,600,
602, 604, 608, 610, 612, 614, 616, 618, 620, 622, 624, 626, 
629, 630, 632, 634, 636, 639, 640, 642, 644, 646, 648, 650, 652, 654, 656, 658, 660, 662, 
664, 666, 668,670, 672, 674, 676, 679, 681, 682, 684, 687, 689,
690, 692, 694, 698, 704,  706, 708,711, 712,714, 716,
718, 720, 722, 724, 726, 728, 731, 732, 734, 736, 739, 740, 742, 744, 746, 748, 751,752,
755, 756, 758, 761,763, 764, 766, 768, 770, 772, 774, 776, 778, 780, 
782, 784, 786, 788, 790, 792, 794, 797, 798, 802, 805, 806, 808, 810, 812, 814, 816,
818, 820, 822, 824, 826, 828, 830, 832, 834, 836, 839, 840, 842, 844, 846, 848, 850,
853, 855, 856, 858, 860, 862, 864, 867, 870,
 872, 875, 877, 879, 881, 883, 884, 888,
890, 892, 894, 898, 900, 902, 904, 907, 908, 910, 912, 915, 916, 918, 921, 922, 924, 927, 929, 
930, 932, 934, 937,
939, 941,942, 944, 946, 948, 950, 952, 954,956, 
958, 961, 963, 965, 966, 968, 970, 972, 974, 976, 979,
980, 982, 984, 989, 990, 993, 994, 996]
path_datos =  "/home/bego/GARROTXA_copia/datos_GARROTXA_resim/"
path_csv =  "/mnt/usb-TOSHIBA_EXTERNAL_USB_20220124010088F-0:0-part2/snapshots/"
datos_edades = pd.read_csv(path_datos + "edades.csv", index_col = 0, sep = ",")
path_results = ""

iniciar = 1

if iniciar == 1:
    dato_vacio = []
    df_vacio = pd.DataFrame(dato_vacio, columns=['Snapshot',  "Lookback", "Minf_15","Minf_20","Minf_25", "Minf_30" ])
    df_vacio.to_csv( "inflow_mass.csv", sep = ",")


for name in snapshots_analysis:
    print(name)
    lb = datos_edades.loc[datos_edades['Snapshot'] == name, 'Lookback'].iloc[0]
    gas = pd.read_csv(path_csv +"Gas_%s.csv" %name, sep = ",")
    gas["R_sph"] = np.sqrt(gas["X"]**2 + gas["Y"]**2 + gas["Z"]**2)
    gas_cold = gas[(gas['Temperature']<10000)].copy()


    gas_15 = gas[(gas["R_sph"]<15)&(gas["R_sph"]>14)&(gas["Vr"]<0)]
    mass_15 = np.sum(gas_15["Mass"])

    gas_20 = gas[(gas["R_sph"]<20)&(gas["R_sph"]>19)&(gas["Vr"]<0)]
    mass_20 = np.sum(gas_20["Mass"])

    gas_25 = gas[(gas["R_sph"]<25)&(gas["R_sph"]>24)&(gas["Vr"]<0)]
    mass_25 = np.sum(gas_25["Mass"])

    gas_30 = gas[(gas["R_sph"]<30)&(gas["R_sph"]>29)&(gas["Vr"]<0)]
    mass_30 = np.sum(gas_30["Mass"])


    datos_SFH = pd.read_csv(path_results + "inflow_mass.csv", index_col = 0, sep = ",")
    new_row = {'Snapshot':name,  "Lookback": lb,  "Minf_15": mass_15,"Minf_20":mass_20, "Minf_25":mass_25 , "Minf_30":mass_30  }#
    datos_SFH = datos_SFH.append(new_row, ignore_index = True)
    datos_SFH.to_csv(path_results + "inflow_mass.csv", sep = ",")