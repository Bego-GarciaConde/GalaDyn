U
    ��[c�/  �                   @   sb   d dl mZ d dlZd dlZd dlZd dlmZ d dl	T d dl
mZ d dlmZ G dd� d�ZdS )�    )�	dataclassN)�Pool)�*)�Snapshotc                   @   sT   e Zd Zddd�Zddd�Zd	d
� Zdd� Zddd�Zddd�Zdd� Z	dd� Z
dS )�Fourier�   r   �   c                 C   s(   || _ || _|| _|| _|| _|| _d S )N)�maximo�minimo�nbins�maxmode�snapshots_analysis�lookback)�selfr   r   r	   r
   r   r   � r   �,/home/bego/GARROTXA/GalaDyn/fourier_utils.py�__init__   s    zFourier.__init__Nc              
   C   s�  t �| jd | jf�}t �| jd | jf�}| j| j | j }t �| j| j| |�}|dd� |dd�  d }t j| jd t jd�}	t j| jd t jd�}
t �|d |d  �}t �	||�d }t �| j�}t
d| jd �D �]�}t
d| j�D �]�}|||k }|||k }t �||�}|dk�r�|dk�r`t �t �|| ��|	|< t �t �|| ��|
|< n8dt �t �|| �� |	|< dt �t �|| �� |
|< n�|||k }|dk�r�t �|t �|| � �|	|< t �|t �|| � �|
|< n@dt �|t �|| � � |	|< dt �|t �|| � � |
|< t �|	| d |
| d  �|||f< |dk�rxt �|
| |	| �|||f< q�|dkr�d|||f< t|�||< q�q�||||fS )z�  This function divides the disk in radial and operates fourier modes from 0 to maxmode 
            in galactic disks, minimo and maximo measured in kpc  �   N�����r   ��dtyper   )�np�zerosr   r   r	   r
   �arange�float32�sqrt�digitize�range�arctan2�sum�cos�sin�len)r   �X�Y�peso�AAZarmangleZsteparmZradarmZrcenter�A�B�dd�indr�
nparticles�m�i�X_i�Y_i�aZpeso_ir   r   r   �fourier_method   sD    


  &
zFourier.fourier_methodc              
   C   s^  t �tt�| j dd| j  f�}d}tt�D �]\}}t|� tj	t
d|� d|� d|� d� dd	�}d
}t|d �D ]F\}	}
|d |	 |kr�||d |	< qv|d |	 | k rv| |d |	< qv| j|d |d |d d�\}}}}t| j�D ]T}t| t| || || gt|d d �|f � t|d d �|f � ||< |d }q�q.| j|d|� �dd� d S )N�   r   r   �mesh_aceleracion_�_�_ytRS_�.csv�,��sepg [n�+=�azr#   r$   �r%   r   Zacceleration_��etiquetar%   )r   r   r"   r   r   r   �	enumerate�print�pd�read_csv�path_accelerationr1   r   r   �list�save_fourierogram)r   �compr   �datos�index�t�name�dfZlimite�j�la�Rcenters�npart�
amplitudes�phasesr-   r   r   r   �apply_fourier_accelerationsS   s     "&$Hz#Fourier.apply_fourier_accelerationsc                 C   s�  t �tt�| j dd| j  f�}t �tt�| j dd| j  f�}d}tt�D �]\}}tjt	d|� d|� d� dd�}| j
|d	 |d
 |d d�\}}	}
}| j
|d	 |d
 |d d�\}}}}t| j�D ]�}t| t| || |	| gt|d d �|f � t|d d �|f � ||< t| t| || || gt|
d d �|f � t|d d �|f � ||< |d }q�qP| j|d|� �dd� | j|d|� �dd� d S )Nr2   r   r   r3   r4   z_satellites_id_ytRS.csvr7   r8   r#   r$   �	az_streamr;   �az_corer   Z	sat_prog_r<   Zsat_streams_)r   r   r"   r   r   r   r>   r@   rA   rB   r1   r   r   rC   rD   )r   �sat_nameZdatos_cZdatos_srG   rH   rI   rJ   Z
Rcenters_sZnpart_sZAmp_sZphases_sZ
Rcenters_cZnpart_cZAmp_cZphases_cr-   r   r   r   �apply_fourier_satl   s    "" $$HHzFourier.apply_fourier_satFc                 C   s�  t d|� �� tjtt�| j dd| j  ftjd�}d}tt�D �]>\}}t |� d}t	|�}|�
�  |��  |�� }	|	|	d dk |	d dk@  }
t d	� |d kr�| j|
d
 |
d d�\}}}}n\|dkr�| j|
d
 |
d t�|
|�  �d�\}}}}n&| j|
d
 |
d |
|�  d�\}}}}t| j�D ]V}t| t| || || gt|d d �|f � t|d d �|f � ||< |d }�q*qB| j|||d� d S )N�Analyzing fourierograms of r2   r   r   r   Z"disk_filtered_All_age_2kpc_Z_limit�Z������Snapshot loaded!r#   r$   �r#   r$   T�r#   r$   r%   r   r<   )r?   r   r   r"   r   r   r   r   r>   r   �
load_stars�	load_disk�filter_disk_particlesr1   �absr   r   rC   rD   )r   r%   �	breathingrF   rG   rH   rI   r=   �snapshot�dfArJ   rM   r+   rO   rP   r-   r   r   r   �apply_fourier_on_disk�   s*    ( .&HzFourier.apply_fourier_on_disk�starsc                 C   s�  t d|� �� tjtt�| j dd| j  ftjd�}d}d}tt�D �]�\}}t |� t	|�}|dkr�|�
�  |��  |j|jd �|jd �  }	t�|	d	 d
�}
|	|	d	 |
k  �� }	|	|	d dk |	d dk @ |	d dk@  }	|dk�r$|��  |j|jd dk |jd dk @ |jd dk@  }	t d� |d k�rV| j|	d |	d d�\}}}}n&| j|	d |	d |	|�  d�\}}}}t| j�D ]V}t| t| || || gt|d d �|f � t|d d �|f � ||< |d }�q�qF| j|||d� d S )NrV   r2   r   r   Zdm_innerr   rd   �IDZAge�#   �R�   rW   �   ������dm�(   �
   i����z
Bar loadedr#   r$   rZ   r[   r   r<   )r?   r   r   r"   r   r   r   r   r>   r   r\   r]   rd   �isin�disk�
percentile�copy�load_dmrk   r1   r   r   rC   rD   )r   �stars_or_dmr%   rF   r=   rG   rH   rI   ra   rJ   Zedad_arM   r+   rO   rP   r-   r   r   r   �apply_fourier_on_bar�   s2    ((
0
 &HzFourier.apply_fourier_on_barc                 C   s�  t d� tjtt�| j dd| j  ftjd�}tjtt�| j dd| j  ftjd�}d}d}tt�D �]X\}}t |� t	|�}|�
�  |��  |��  |��  |j}|��  t d� | j|d |d |d	 d
�\}	}
}}t| j�D ]T}t| t| |	| |
| gt|d d �|f � t|d d �|f � ||< |d }q�| j|d |d |d d
�\}	}
}}t| j�D ]V}t| t| |	| |
| gt|d d �|f � t|d d �|f � ||< |d }�qjqh| j|ddd� | j|ddd� d S )Nz0Analyzing fourierograms of bending and breathingr2   r   r   r   rY   r#   r$   �Bendingr[   r   �	Breathing� Zbendingr<   r`   )r?   r   r   r"   r   r   r   r   r>   r   r\   r]   r^   �calculate_bending_breathing�bending_breathing_mode�plot_bending_breathingr1   r   r   rC   rD   )r   Zdatos_bendingZdatos_breathingZindex_bendingZindex_breathingrH   rI   ra   rJ   rM   r+   rO   rP   r-   r   r   r   �"apply_fourier_on_bending_breathing�   s2    (($H
$Hz*Fourier.apply_fourier_on_bending_breathingc              
   C   s�   ddddg}t d| jd �D ]}|�d|� �� qt d| jd �D ]}|�d|� �� qBtj||d	�}td
t� d| j� d|� d|� d�	� |jtd| j� d|� d|� d� ddd� d S )NZ
snapshot_tZlookbacktimerM   Z
Nparticlesr   r   �amp�phase)�columnszSaving data as z	 fourier_r4   r6   Zfourier_r7   F)r9   rG   )	r   r   �appendr@   �	DataFramer?   �path_resultsr   �to_csv)r   rF   r=   r%   �column_namesZnumero_modor   r   r   rD   �   s    $zFourier.save_fourierogram)r   r   r   r   )N)NF)rd   N)�__name__�
__module__�__qualname__r   r1   rQ   rU   rc   rt   r{   rD   r   r   r   r   r   
   s   


>
$
)&r   )�dataclassesr   �numpyr   �pandasr@   �gc�multiprocessingr   �config�matplotlib.pyplot�pyplot�plt�snapshot_definitionr   r   r   r   r   r   �<module>   s   