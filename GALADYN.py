#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 24/04/2022
@author: bego
"""

from dynamics_functions import*
from accelerations import * 
from config import *
import ConfigParser


        ##################################################
        #             GALACTIC DYNAMICS                  #
        #   This code makes a complete analysis of       #
        #     the dynamics in GARROTXA model,
        #  measures accelerations of all components      #
        ##################################################

        ##################################################
        #                     Settings                   #
        ##################################################

# ac_satellites_as_points =     0 
# ac_satellites_ids =           0  
# ac_dm =                       0  
# ac_gas =                      0  
# ac_stars =                    0 
# #-------
# fourier_analysis_z  =         0
# fourier_analysis_vz =         0
# fourier_analysis_dm =         0
# fourier_analysis_gas =        0
# #-------
# comparison_plots =            0

# snapshots_analysis = [ 602,604, 608, 610, 612, 614, 616, 618, 620, 622, 624, 626, 
# 629, 630, 632, 634, 636, 639, 640, 642, 644, 646, 648, 650, 652, 654, 656, 658, 660, 662, 
# 664, 666, 668,670, 672, 674, 676, 679, 681, 682, 684, 687, 689,
# 690, 692, 694, 698, 704,  706, 708,711, 712,714, 716, 718, 720, 
# 722, 724, 726, 728, 731, 732, 734, 736, 739, 740, 742, 744, 746, 748, 751,752,
#  755, 756, 758, 761,763, 764, 766, 768, 770, 772, 774, 776, 778, 780, 
# 782, 784, 786, 788, 790, 792, 794, 797, 798, 802, 805, 806, 808, 810, 812, 814, 816,
#  818, 820, 822, 824, 826, 828, 
# 830, 832, 834, 836, 839, 840, 842, 844, 846, 848, 850,
# 853, 855, 856, 858, 860, 862, 864, 867, 870, 872, 875, 877, 879, 881, 883, 884, 886, 888,
# 890, 892, 894, 898, 900, 902, 904, 907, 908, 910, 912, 915, 
# 916, 918, 921, 922, 924, 927, 929, 
# 930, 932, 934, 937, 939, 941,942, 944, 946, 948, 950, 952, 954,956, 
# 958, 961, 963, 965, 966, 968, 970, 972, 974, 976, 979,
#  980, 982, 984, 986, 989, 990, 993, 994, 996 ]

# # ---------------------------------------------------------------------------
# #path_satellite_models = "/media/temp1/bego/snapshots/modelos_satelites/"
# path_csv = "/home/bego/GARROTXA/snapshots/"
# path_datos = "/home/bego/GARROTXA_copia/datos_GARROTXA_resim/"
# path_crossmatch = "/home/bego/GARROTXA/satelites_crossmatch/"
# path_figures = "/home/bego/GARROTXA/aceleration_figures/"
# path_mesh = "/home/bego/GARROTXA/acceleration/"
#----------------------------------------------------------------------------

def main():
    
    #TO DO: INICIALIZACION, COMPROBAR QUE EXISTEN TODOS LOS DIRECTORIOS 
    #TO DO: LEER LA CONFIGURACIÓN 
#    read_parameters()
    #COMPROBAR QUE EXISTE EL CROSSMATCHS Y ARCHIVOS DE CADA COMPONENTE
    for name in snapshots_analysis:
    #First steps accelerations
   # accelerations ()
        acceleration_in_mesh_comp (name,"dm",None, path_disk, plot = True, tidal = False)
        acceleration_in_mesh_comp (name,"gas",None, path_disk, plot = True, tidal = False)
        acceleration_in_mesh_comp (name,"stars","disk", path_disk, plot = True, tidal = False)
        acceleration_in_mesh_comp (name,"stars","nodisk", path_disk, plot = True, tidal = False)

    #Second steps: fourier


    #Third step: comparison


        
if __name__ == "__main__":
    main()