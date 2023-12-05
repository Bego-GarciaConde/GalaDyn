import pandas as pd
import numpy as np

        ##################################################
        #                     Settings                   #
        ##################################################

CALCULATE_2ND_ALIGNMENT = 0
APPLY_ALIGNMENT = 1
ALIGN_RADIUS = 7  #7 kpc or 15kpc

AC_CALCULATE = 0
AC_SATELLITES_AS_POINTS =     0 
#AC_SATELLITES_IDS =           1  
AC_DM =                        0
AC_GAS =                      0
AC_STARS =                    0
AC_STARS_ELLIPSOID =           0
AC_ALL_SAT =                  0

LIMIT = 40 #LIMIT FOR ACCELERATION MESH
G= 1.3273E11
KPC_TO_KM= 3.086E16
SECONDS_TO_MYR = 3.15576E+16
SOFTENING = 0

#-------
tag= f"align_{ALIGN_RADIUS}kpc"

FOURIER_ACCELERATION_DM=      0
FOURIER_ACCELERATION_GAS=     0
FOURIER_ACCELERATION_SUM =      0
FOURIER_ACCELERATION_TOTAL=    0
FOURIER_ACCELERATION_DM_INNER = 0
FOURIER_ACCELERATION_STARS=     0
FOURIER_ACCELERATION_STELLAR_ELLIPSOID=     0
FOURIER_ACCELERATION_SATELLITES = 0
TORQUE = 0

#-------
FOURIER_BENDING_BREATHING =   0
FOURIER_DENSITY =             0
FOURIER_MASS =                1
FOURIER_Z =                   1
FOURIER_VZ =                  1
FOURIER_VR =                  1 
FOURIER_VPHI =                1  

FOURIER_BAR  =                0
FOURIER_DM   =                 0
#-------





snapshots_analysis= [
 #       520,523,525, 527,530,532,535, 537,539,541,
 #543, 545,547, 550, 553, 555,557, 
 #560, 563, 565, 567,570,573, 575, 577, 580,
 #583, 585,587,590, 592,594,596,598,

 600, 602, 604, 608, 610, 612, 614, 616, 618,620,
 622, 624, 626, 629, 630, 632, 634, 636, 639, 640, 642, 644, 646, 648, 650, 652, 654, 656, 658, 660, 662, 
 664, 666,668,670, 672, 674, 676, 679, 681, 682, 684, 687, 689,
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
980, 982, 984, 989, 990, 993, 994, 996, 999]


LB_REF = 0.01321

# ---------------------------------------------------------------------------
#path_satellite_models = "/media/temp1/bego/snapshots/modelos_satelites/"
#PATH_CSV = "/mnt/usb-TOSHIBA_EXTERNAL_USB_20220124010088F-0:0-part2/snapshots_resim_new/"
PATH_CSV = "/mnt/usb-TOSHIBA_EXTERNAL_USB_20220124010088F-0:0-part2/snapshots_resim_1step/"
PATH_DATOS = "/home/bego/GARROTXA_copia/datos_GARROTXA_resim/"
PATH_2ND_ALIGNMENT = "/home/bego/GARROTXA/rotation_matrix_1step/"
PATH_CROSSMATCH = "/mnt/usb-TOSHIBA_EXTERNAL_USB_20220124010088F-0:0-part2/satelites_crossmatch/"
PATH_ACCELERATION = f"/home/bego/GARROTXA/acceleration_1step_{ALIGN_RADIUS}kpc/"
if ALIGN_RADIUS == 7:
        PATH_DISK = "/mnt/usb-TOSHIBA_EXTERNAL_USB_20220124010088F-0:0-part2/disk_1step/" #Aligned with 7 kpc
elif ALIGN_RADIUS == 15:
        PATH_DISK = "/home/bego/GARROTXA/disco/"#Aligned with 15 kpc
PATH_RESULTS = f"/home/bego/GARROTXA/GalaDyn/results/results_{ALIGN_RADIUS}kpc/"
PATH_FIGURES_BENDING = "/home/bego/GARROTXA/BendingBreathing/"

#----------------------------------------------------------------------------
satelites = ["arania", "grillo", "mosquito", "all"]


datos_edades = pd.read_csv(PATH_DATOS + "edades.csv", sep = ",",index_col = 0)
lookback = [datos_edades.loc[datos_edades['Snapshot'] == name, 'Lookback'].iloc[0] for name in snapshots_analysis]