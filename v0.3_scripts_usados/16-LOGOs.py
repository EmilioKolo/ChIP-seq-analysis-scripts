
ib3 = False; 
http_proxy_ib3 = 'http://proxy.fcen.uba.ar:8080'; 
https_proxy_ib3 = 'https://proxy.fcen.uba.ar:8080'; 
ftp_proxy_ib3 = 'ftp://proxy.fcen.uba.ar:8080'; 

# Generales
import os
import time
import copy
from random import shuffle
from pyensembl import EnsemblRelease

# Analisis de secuencias y genomas
from Bio import Entrez, SeqIO, motifs, Seq
Entrez.email = 'ekolomenski@gmail.com'; 
Entrez.api_key = '9bccbd4949f226a5930868df36c211e8b408'; 
import biomart

# Generacion de figuras
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import logomaker as lm

# Importo analisischip
import sys
path_analisischip_casa = 'C:\\Users\\Admin\\Documents\\Scripts\\Git\\AnalisisSecuenciasChIP\\'; 
path_analisischip_ib3 = 'C:\\Users\\emili\\Documents\\Archivos nube\\Git\\analisissecuenciaschip\\'; 
if ib3:
    path_analisischip_main = path_analisischip_ib3; 
else:
    path_analisischip_main = path_analisischip_casa; 
sys.path.append(path_analisischip_main); 
from analisischip.seq import seq_data, seq_handler

# Genomas de referencia
hg19 = EnsemblRelease(75, species='human'); 
hg38 = EnsemblRelease(102, species='human'); 
mm9 = EnsemblRelease(54, species='mouse'); 
mm10 = EnsemblRelease(102, species='mouse'); 


'''
Funciones para hacer LOGOs con secuencias de ADN
'''

#################################### VARIABLES ####################################

# Direcciones de archivos y outputs para ib3 y casa
path_dropbox_casa = 'D:\\Users\\Admin\\Dropbox\\'; 
path_dropbox_ib3 = 'C:\\Users\\emili\\Dropbox\\'; 
path_fasta_casa = 'D:\\ArchivosDoctorado\\Genomas\\'; 
path_fasta_ib3 = 'X:\\Genomas\\'; 
path_output_dump_casa = 'D:\\ArchivosDoctorado\\Output_dump\\'; 
path_output_dump_ib3 = 'X:\\Output_dump\\'; 



### Variables main()

# Path main depende de si estoy en ib3 o casa
if ib3:
    path_dropbox_main = path_dropbox_ib3; 
    path_fasta_main = path_fasta_ib3; 
    path_output_dump_main = path_output_dump_ib3; 
else:
    path_dropbox_main = path_dropbox_casa; 
    path_fasta_main = path_fasta_casa; 
    path_output_dump_main = path_output_dump_casa; 
path_bed_main = path_dropbox_main + 'Doctorado\\3-Genes transactivados rio abajo\\0-Fuentes\\Papers ChIP-seq\\'; 
path_trabajo_main = path_dropbox_main + 'Doctorado\\3-Genes transactivados rio abajo\\1-Interpretacion de datos ChIP-seq\\'; 
path_in_main = path_trabajo_main + ''; # Uso path_pwm de human y mouse
path_out_main = path_output_dump_main + 'LOGOs\\'; 
# Path que dependen de la especie
path_pwm_human = path_dropbox_main + 'Doctorado\\3-Genes transactivados rio abajo\\1-Interpretacion de datos ChIP-seq\\PWM_human\\'; 
path_pwm_mouse = path_dropbox_main + 'Doctorado\\3-Genes transactivados rio abajo\\1-Interpretacion de datos ChIP-seq\\PWM_mouse\\'; 

L_sitios = ['GCAAGTG', 'GGAAGTG', 'GAAAGTG', 'ATAAGTG', 'GTAAGTG', 'CTAAGTG', 'TCAAGTG', 'TGAAGTG', 'TAAAGTG', 'TTAAGTG']; 

# Nombres de archivos de matrices de peso y factores de transcripcion asociados
L_arch_pwm_raton = ['NKX25_MOUSE.H11MO.0.A.pcm', 'TBX20_MOUSE.H11MO.0.C.pcm', 'MEIS1_MOUSE.H11MO.1.A.pcm', 'TGIF1_MOUSE.H11MO.0.A.pcm', 
                    'HAND1_MOUSE.H11MO.0.C.pcm', 'MAF_MOUSE.H11MO.1.A.pcm', 'GATA1_MOUSE.H11MO.1.A.pcm', 'GATA6_MOUSE.H11MO.0.A.pcm', 'GATA4_MOUSE.H11MO.0.A.pcm']; 
L_arch_pwm_humano = ['NKX25_HUMAN.H11MO.0.B.pcm', 'TBX20_HUMAN.H11MO.0.D.pcm', 'MEIS1_HUMAN.H11MO.1.B.pcm', 'TGIF1_HUMAN.H11MO.0.A.pcm', 'HAND1_HUMAN.H11MO.0.D.pcm', 
                     'MAF_HUMAN.H11MO.1.B.pcm', 'GATA1_HUMAN.H11MO.1.A.pcm', 'GATA6_HUMAN.H11MO.0.A.pcm', 'GATA4_HUMAN.H11MO.0.A.pcm']; 
L_nombres_pssm_raton = ['Nkx25', 'Tbx20', 'Meis1', 'Tgif1', 'Hand1', 'Maf', 'Gata1', 'Gata6', 'Gata4']; 
L_nombres_pssm_humano = ['NKX25', 'TBX20', 'MEIS1', 'TGIF1', 'HAND1', 'MAF', 'GATA1', 'GATA6', 'GATA4']; 



#################################### FUNCIONES ####################################


def _main():
    # Funcion para probar funciones en ejecucion del archivo

    # Prueba de generar logos
    if False:
        # Prueba de generar logos
        m_human, pwm_human, pssm_human = abrir_pssm(L_arch_pwm_humano[0], path_arch=path_pwm_human, solo_pssm=False, pseudocounts=0.5); 
        m_mouse, pwm_mouse, pssm_mouse = abrir_pssm(L_arch_pwm_raton[0], path_arch=path_pwm_mouse, solo_pssm=False, pseudocounts=0.5); 
        # Paso m, pwm o pssm a data frame con Pandas
        m_m_df = pd.DataFrame(m_mouse.counts); 
        #m_m_df = pd.DataFrame(pwm_mouse); 
        #m_m_df = pd.DataFrame(pssm_mouse); 
        info_mat = lm.transform_matrix(m_m_df, from_type='counts', to_type='information'); 
        # Inicializo pyplot
        plt.ion(); 
        # Uso lm para hacer LOGO con el dataframe
        logo = lm.Logo(info_mat); 
        plt.show(block=True); 
    
    # Genero y guardo logos de lista
    if False: # Funciones para generar logos de PWM
        # Abro las matrices de peso de cada lista
        for i in range(len(L_arch_pwm_humano)):
            nom_arch_pwm_h = L_arch_pwm_humano[i]; 
            nom_arch_pwm_m = L_arch_pwm_raton[i]; 
            nom_out_pwm_h = L_nombres_pssm_humano[i]; 
            nom_out_pwm_m = L_nombres_pssm_raton[i]; 
            # Consigo las matrices de raton y humano
            m_human, pwm_human, pssm_human = abrir_pssm(nom_arch_pwm_h, path_arch=path_pwm_human, solo_pssm=False, pseudocounts=0.5); 
            m_mouse, pwm_mouse, pssm_mouse = abrir_pssm(nom_arch_pwm_m, path_arch=path_pwm_mouse, solo_pssm=False, pseudocounts=0.5); 
            # Creo los dataframes
            m_h_df = pd.DataFrame(m_human.counts); 
            m_m_df = pd.DataFrame(m_mouse.counts); 
            # Uso generar_logo() para generar los logos de counts
            m_h_df = generar_logo(m_h_df, nom_out_pwm_h+'_human_counts', path_out=path_out_main+'counts\\'); 
            m_m_df = generar_logo(m_m_df, nom_out_pwm_m+'_mouse_counts', path_out=path_out_main+'counts\\'); 
            # Creo las matrices de informacion
            info_mat_h = lm.transform_matrix(m_h_df, from_type='counts', to_type='information'); 
            info_mat_m = lm.transform_matrix(m_m_df, from_type='counts', to_type='information'); 
            # Uso generar_logo() para generar los logos de informacion
            info_mat_h = generar_logo(info_mat_h, nom_out_pwm_h+'_human_info', path_out=path_out_main+'info\\'); 
            info_mat_m = generar_logo(info_mat_m, nom_out_pwm_m+'_mouse_info', path_out=path_out_main+'info\\'); 

    # Genero el logo en base a sitios confirmados
    if False: # Uso pipeline_logo_seq()
        # Defino el nombre de los archivos creados
        nom_out = 'NKX25_sitios_confirmados'; 
        pipeline_logo_seq(L_sitios, nom_out, path_out=path_out_main); 

    # Genero el logo en base a secuencias cerca de genes up o downregulados
    if True: # Uso pipeline_logo_seq() y extraer_secuencias_filtro()
        ### HUMAN
        print('>>> EMPEZANDO ANALISIS DE DATOS EN HUMANO')
        # Genero L_sitios con extraer_secuencias_filtro()
        L_sitios_up = extraer_secuencias_filtro('datos_pca_full_human', 6, [(8,'1')], path_arch=path_in_main); 
        print('Largo L_sitios_up: ' + str(len(L_sitios_up)))
        # Defino el nombre de los archivos creados
        nom_out = 'human_NKX25_sitios_upreg'; 
        pipeline_logo_seq(L_sitios_up, nom_out, path_out=path_out_main); 
        # Genero L_sitios con extraer_secuencias_filtro()
        L_sitios_down = extraer_secuencias_filtro('datos_pca_full_human', 6, [(8,'-1')], path_arch=path_in_main); 
        print('Largo L_sitios_down: ' + str(len(L_sitios_down)))
        # Defino el nombre de los archivos creados
        nom_out = 'human_NKX25_sitios_downreg'; 
        pipeline_logo_seq(L_sitios_down, nom_out, path_out=path_out_main); 
        ### MOUSE
        print('>>> EMPEZANDO ANALISIS DE DATOS EN RATON')
        # Genero L_sitios con extraer_secuencias_filtro()
        L_sitios_up = extraer_secuencias_filtro('datos_pca_full_mouse', 6, [(8,'1')], path_arch=path_in_main); 
        print('Largo L_sitios_up: ' + str(len(L_sitios_up)))
        # Defino el nombre de los archivos creados
        nom_out = 'mouse_NKX25_sitios_upreg'; 
        pipeline_logo_seq(L_sitios_up, nom_out, path_out=path_out_main); 
        # Genero L_sitios con extraer_secuencias_filtro()
        L_sitios_down = extraer_secuencias_filtro('datos_pca_full_mouse', 6, [(8,'-1')], path_arch=path_in_main); 
        print('Largo L_sitios_down: ' + str(len(L_sitios_down)))
        # Defino el nombre de los archivos creados
        nom_out = 'mouse_NKX25_sitios_downreg'; 
        pipeline_logo_seq(L_sitios_down, nom_out, path_out=path_out_main); 
    return ''


def abrir_csv(nom_arch, path_arch='', ext='.csv', sep=';', ignore_headers=True):
    # Abre un archivo en formato de matriz con filas por lineas y columnas separadas por sep
    # Devuelve una matriz con una lista de filas, cada una con una lista de columnas
    # ignore_headers ignora la primer fila

    # Inicializo la matriz que se devuelve
    M_out = []; 
    # Defino la direccion del archivo con nom_arch y path_arch
    if path_arch=='':
        filepath = nom_arch + ext; 
    else:
        filepath = os.path.join(path_arch, nom_arch + ext); 
    # Booleano para revisar que se ignore la primera fila
    header_ignorado = not ignore_headers; 
    # Abro el archivo en modo reading
    with open(filepath, 'r') as F_open:
        print('Archivo ' + nom_arch + ' abierto.')
        # Recorro cada fila
        for curr_line in F_open.readlines():
            if header_ignorado:
                # Transformo la linea en lista
                L_line = curr_line.rstrip().split(sep=sep); 
                # Cargo L_line en M_out
                M_out.append(L_line[:]); 
            else:
                # Si header_ignorado empieza sin ser True, se ignora la primera linea
                header_ignorado = True; 
    return M_out


def abrir_L_pssm(L_nom_arch, path_arch='', pseudocounts=0.5):
    # Extrae una lista de matrices de counts y devuelve sus PSSM correspondientes

    # Inicializo la lista que se devuelve
    L_pssm = []; 
    # Recorro L_nom_arch
    for nom_arch in L_nom_arch:
        # Uso abrir_pssm() con solo_pssm=True para extraer los PSSM
        L_pssm.append(abrir_pssm(nom_arch, path_arch=path_arch, solo_pssm=True, pseudocounts=pseudocounts)); 
    return L_pssm


def abrir_pssm(nom_arch, path_arch='', solo_pssm=True, pseudocounts=0.5):
    # Extrae la matriz de counts de nom_arch y la pasa a PSSM

    # Extraigo el motif con extraer_motif_pcm()
    m = extraer_motif_pcm(nom_arch, path_arch=path_arch); 
    # Reviso que m no sea lista vacia y genero pwm y pssm
    if m != '':
        pwm = m.counts.normalize(pseudocounts=pseudocounts); 
        pssm = pwm.log_odds(); 
    else:
        pwm = ''; 
        pssm = ''; 
    # Devuelvo solo pssm por default
    if solo_pssm:
        return pssm
    # Si solo_pssm==False, devuelvo m y pwm tambien
    else:
        return m, pwm, pssm
    

def extraer_motif_pcm(nom_arch, path_arch=''):
    # Abre archivo con motif en formato pfm-four-columns y devuelve el primer motif

    # Defino dir_arch
    if path_arch != '':
        dir_arch = os.path.join(path_arch, nom_arch); 
    else:
        dir_arch = str(nom_arch); 
    # Abro pcm con motifs.parse
    handle = open(dir_arch); 
    record = motifs.parse(handle, 'pfm-four-columns'); 
    handle.close(); 
    # Veo que record tenga un solo motif
    if len(record) == 1:
        m = record[0]; 
    elif len(record) > 1:
        print('WARNING: ' + nom_arch + ' tiene ' + str(len(record)) + ' motifs. Se devuelve solo el primero.')
        m = record[0]; 
        print('Motifs en record:')
        for i in record:
            print(i.name)
            print(i.counts)
    else:
        print('ERROR: No se encontro ningun motif en ' + nom_arch + '. Se devuelve string vacio.')
        m = ''; 
    return m


def extraer_secuencias_filtro(nom_arch, col_seq, M_filt=[], path_arch='', ext_arch='.csv', sep_arch=';', ignore_headers=True):
    # Funcion para extraer secuencias de un archivo .csv
    # Selecciona secuencias en columna col_seq
    # Filtra filas de acuerdo a M_filt, que contiene tuplas de formato (col_filt, val_filt)
    # Solo agarra secuencias de filas que tengan el valor val_filt en la columna col_filt

    # Inicializo la lista que se devuelve
    L_out = []; 
    # Abro la matriz en nom_arch
    M_csv = abrir_csv(nom_arch, path_arch=path_arch, ext=ext_arch, sep=sep_arch, ignore_headers=ignore_headers); 
    # Recorro M_csv
    for i in range(len(M_csv)):
        curr_row = M_csv[i]; 
        # Inicializo el booleano que define si me quedo con la secuencia de curr_row
        pasa_filtro = True; 
        # Recorro M_filt
        for j in range(len(M_filt)):
            # Defino columna y valor
            curr_col = int(M_filt[j][0]); 
            curr_val = str(M_filt[j][1]); 
            # Veo que curr_row[curr_col] sea igual a curr_val
            if str(curr_row[curr_col])!=curr_val:
                # Si cualquiera de los filtros no da, el filtro queda False
                pasa_filtro=False; 
        # Veo si la fila pasa el filtro
        if pasa_filtro:
            # Si pasa el filtro, anoto la secuencia en L_out
            L_out.append(str(curr_row[col_seq])); 
    return L_out


def generar_logo(df, nom_out, path_out='.\\', ext='.png'):
    # Funcion para generar LOGO en base a una matriz de pesos en formato Bio.motif

    # Defino dirarch
    dirarch = os.path.join(path_out, nom_out + ext); 
    # Inicializo plt
    plt.ion(); 
    # Uso lm para hacer LOGO con m_df
    logo = lm.Logo(df); 
    # Guardo el plot
    plt.savefig(dirarch); 
    # Cierro la figura (no se si sirve)
    plt.close(); 
    return df


def pipeline_logo_seq(L_sitios_usado, nom_out, path_out='.\\'):
    # Funcion para correr los scripts necesarios para crear LOGOs en base a una lista de secuencias
    # Crea dos LOGOs, uno con counts y otro con informacion

    # Creo un dataframe con L_sitios_usado
    counts_mat_conf = lm.alignment_to_matrix(L_sitios_usado); 
    # Uso generar_logo() para generar los logos de counts
    counts_mat_conf = generar_logo(counts_mat_conf, nom_out+'_counts', path_out=path_out); 
    # Creo las matrices de informacion
    info_mat_conf = lm.transform_matrix(counts_mat_conf, from_type='counts', to_type='information'); 
    # Uso generar_logo() para generar los logos de informacion
    info_mat_conf = generar_logo(info_mat_conf, nom_out+'_info', path_out=path_out); 


#################################### RESULTADOS ###################################

output_dump = []; 

if __name__=='__main__':
    output_dump.append(_main()); 

