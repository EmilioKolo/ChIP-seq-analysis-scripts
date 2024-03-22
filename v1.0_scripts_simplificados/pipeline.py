
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
# Pasar a formato usado en .csv
    # chr_n, contig, pos_ini, pos_end
    # Sitios de union encontrados (por secuencia y por matriz)
        # Genes cercanos
        # Sitios de union otros TF cercanos
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

# Variables para pipeline
dist_max_main = 1000000; 
L_confirmados = ['GCAAGTG', 'GGAAGTG', 'GAAAGTG', 'ATAAGTG', 'GTAAGTG', 'CTAAGTG', 'TCAAGTG', 'TGAAGTG', 'TAAAGTG', 'TTAAGTG']; 



#################################### FUNCIONES ####################################


def _main():
    '''Funcion para probar scripts en ejecucion del archivo'''
    return ''


def pipeline_generador(nom_bed, genoma_ensembl, pssm, score_max, path_bed='', dist_max_gen=1000000, L_su=[]):
    '''Pipeline desde abrir archivos .bed y resultados RNA-seq hasta generar tablas de sitios y genes que se usan por otros scripts.
    genoma_ensembl es el elemento genoma de ensembl, usado para 
    dist_max_gen es la distancia maxima a la que se buscan genes cercanos.
    L_su es una lista de secuencias consideradas sitios de union.
    pssm es una matriz de puntaje para secuencias (FORMATO UNIFICADO POR UNA FUNCION)
    score_max es el score maximo aceptado para una secuencia correcta de pssm (uso 90% del puntaje maximo posible).'''

    # Abro archivo .bed
    M_peaks = abrir_bed(nom_bed, path_arch=path_bed); 
    # Inicializo la matriz de sitios de union
    M_su = []; 
    # Inicializo la matriz de genes
    M_genes = []; 
    # Recorro M_peaks para trabajar sobre cada peak
    for i in range(len(M_peaks)):
        curr_peak = M_peaks[i]; 
        # Uso genes_cercanos_peak() para buscar genes cercanos
        M_genes = M_genes + genes_cercanos_peak(curr_peak[0], curr_peak[1], curr_peak[2], dist_max_gen, genoma_ensembl); 
        # Uso secuencia_peak() para obtener la secuencia
        seq_peak = secuencia_peak(curr_peak[0], curr_peak[1], curr_peak[2]); 
        # Uso sitios_union_lista() para buscar sitios de union por lista
        M_su = M_su + sitios_union_lista(seq_peak, L_su); 
        # Uso sitios_union_pssm() para buscar sitios de union por pssm
        M_su = M_su + sitios_union_pssm(seq_peak, pssm, score_max); 
    ### FALTA
    # Recorrer M_su para chequeos (capaz es innecesario)
    # Guardar M_peaks y M_su
    ###
    return M_peaks, M_su, M_genes


def genes_cercanos_peak(chr_n, pos_ini, pos_end, dist_max, genoma_ensembl):
    '''Funcion para buscar los genes cercanos al rango definido por chr_n, pos_ini y pos_end en genoma_ensembl.
    Devuelve listas de genes, donde cada gen es una lista con ID, chr_n, pos0, tipo (prot_coding/otro) (ver si agrego algo mas).'''

    # Inicializo M_out
    M_out = []; 
    # Defino contig
    contig = chrn_to_contig(chr_n); 
    # Busco genes alrededor del peak con genoma_ensembl.genes_at_locus()
    L_genes_cerca_raw = genoma_ensembl.genes_at_locus(contig, pos_ini-dist_max, pos_end+dist_max); 
    # Recorro L_genes_cerca_raw
    for i in range(len(L_genes_cerca_raw)):
        curr_gen = L_genes_cerca_raw[i]; 
        # Inicializo L_gen, que se agrega a M_out
        L_gen = []; 
        # Defino pos0 con definir_pos0()
        pos0 = definir_pos0(curr_gen.start, curr_gen.end, curr_gen.strand); 
        # Cargo datos relevantes a L_gen
        L_gen.append(curr_gen.gene_id); # ID del gen para volver a buscarlo
        L_gen.append(chr_n); # chr_n para tenerlo en formato igual que el resto de scripts
        L_gen.append(pos0); # pos0 determina el +1 del gen
        L_gen.append(curr_gen.biotype); # Me importan mas los protein_coding que el resto
        # Cargo L_gen a M_out
        M_out.append(L_gen[:]); 
    return M_out


def secuencia_peak(chr_n, pos_ini, pos_end):
    '''Funcion para obtener la secuencia del rango definido por chr_n, pos_ini y pos_end.
    Devuelve un string con una secuencia de ADN.'''
    ### FALTA
    pass


def sitios_union_lista(seq_peak, L_sitios):
    '''Funcion para encontrar todas las ocurrencias de cada una de las secuencias en L_sitios dentro de seq_peak.
    Devuelve listas de sitios de union, donde cada sitio incluye chr_n, pos_ini, pos_end, seq y fuente="lista" (ver si agrego algo mas).'''
    ### FALTA
    pass


def sitios_union_pssm(seq_peak, pssm, score_max):
    '''Funcion para encontrar todas las secuencias con puntaje mayor a score_max dentro de seq_peak.
    Devuelve listas de sitios de union, donde cada sitio incluye chr_n, pos_ini, pos_end, seq y fuente="pssm" (ver si agrego algo mas).'''
    ### FALTA
    pass


def chrn_to_contig(chr_n):
    '''Funcion para pasar de chr_n a contig.
    Recibe un string de formato "chr"+N, donde N puede ser un numero o letras (X, Y, M).
    Devuelve un string. con solo N'''
    return chr_n[3:]


def definir_pos0(start, end, strand):
    '''Funcion que recibe el start, end y strand de un elemento gen de EnsemblRelease y devuelve el +1.'''

    # Veo si strand es + o -
    if strand == '+':
        ret = start; 
    elif strand == '-':
        ret = end; 
    # Si no es ninguno de los dos, devuelvo gene_start por defecto
    else:
        print('WARNING: Strand ' + str(strand) + ' no reconocido. Se devuelve start como +1.')
        ret = start; 
    return ret


def abrir_bed(nom_arch, path_arch='', sep='\t', ext='.bed', col_num=3):
    '''Funcion que abre archivos .bed y devuelve una matriz con una lista de peaks.
    col_num determina el numero de columnas que se devuelven. Por defecto, solo recopilo chr_n, pos_ini y pos_end.'''
    
    # Inicializo la matriz que se devuelve
    M_out = []; 
    # Defino la direccion del archivo con nom_arch y dir_arch
    if path_arch=='':
        filepath = nom_arch + ext; 
    else:
        filepath = os.path.join(path_arch, nom_arch + ext); 
    # Abro el archivo filepath
    with open(filepath, 'r') as f_out:
        print('Archivo ' + nom_arch + ' abierto.')
        # Recorro cada fila
        for curr_line in f_out.readlines():
            # Transformo la linea en lista
            L_line = curr_line.rstrip().split(sep=sep); 
            # Cargo L_line en M_out
            M_out.append(L_line[:]); 
    return M_out


#################################### RESULTADOS ###################################

output_dump = []; 

if __name__=='__main__':
    output_dump.append(_main()); 

