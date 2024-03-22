
ib3 = False; 
http_proxy_ib3 = 'http://proxy.fcen.uba.ar:8080'; 
https_proxy_ib3 = 'https://proxy.fcen.uba.ar:8080'; 
ftp_proxy_ib3 = 'ftp://proxy.fcen.uba.ar:8080'; 

# Generales
import os 
import time 
import math 
import sys 
import pandas as pd 
import numpy as np 

# Para debugging
from random import shuffle 

# Especificos bioinfo
from pyensembl import EnsemblRelease 
from Bio import Entrez, SeqIO, motifs, Seq 
Entrez.email = 'ekolomenski@gmail.com'; 
Entrez.api_key = '9bccbd4949f226a5930868df36c211e8b408'; 
import biomart 

# Importo analisischip
path_analisischip_casa = 'C:\\Users\\Admin\\Documents\\Scripts\\Git\\AnalisisSecuenciasChIP\\'; 
path_analisischip_ib3 = 'C:\\Users\\emili\\Documents\\Archivos nube\\Git\\analisissecuenciaschip\\'; 
if ib3:
    path_analisischip_main = path_analisischip_ib3; 
else:
    path_analisischip_main = path_analisischip_casa; 
sys.path.insert(1, path_analisischip_main); 
#from analisischip.seq import seq_data, seq_handler 

# Graficos
import matplotlib.pyplot as plt 
import seaborn as sns 

# Genomas de referencia
hg19 = EnsemblRelease(75, species='human'); 
hg38 = EnsemblRelease(102, species='human'); 
mm9 = EnsemblRelease(54, species='mouse'); 
mm10 = EnsemblRelease(102, species='mouse'); 


'''
### FALTA
# Abrir archivos .bed
# 
'''


#################################### VARIABLES ####################################

# Direcciones de archivos y outputs para ib3 y casa
path_git_casa = 'C:\\Users\\Admin\\Documents\\Scripts\\Git\\ChIP-seq-analysis-scripts\\v1.0_scripts_simplificados\\'; 
path_git_ib3 = 'C:\\Users\\emili\\Documents\\Archivos nube\\Git\\ChIP-seq-analysis-scripts\\v1.0_scripts_simplificados\\'; 
path_fasta_casa = 'D:\\ArchivosDoctorado\\Genomas\\'; 
path_fasta_ib3 = 'X:\\Genomas\\'; 
path_output_dump_casa = 'D:\\ArchivosDoctorado\\Output_dump\\'; 
path_output_dump_ib3 = 'X:\\Output_dump\\'; 

### Variables main()

# Path main depende de si estoy en ib3 o casa
if ib3:
    path_git_main = path_git_ib3; 
    path_fasta_main = path_fasta_ib3; 
    path_output_dump_main = path_output_dump_ib3; 
else:
    path_git_main = path_git_casa; 
    path_fasta_main = path_fasta_casa; 
    path_output_dump_main = path_output_dump_casa; 
path_in_main = path_git_main; ### Cambiar si es necesario
path_out_main = path_output_dump_main + 'output_git\\'; 
# Path que dependen de la especie
path_pwm_human = path_git_main + 'PWM_human\\'; 
path_pwm_mouse = path_git_main + 'PWM_mouse\\'; 


#################################### FUNCIONES ####################################


def _main():
    # Funcion para probar funciones en ejecucion del archivo
    return ''


#################################### RESULTADOS ###################################

output_dump = []; 

if __name__=='__main__':
    output_dump.append(_main()); 

