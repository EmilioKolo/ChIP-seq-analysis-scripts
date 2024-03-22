
ib3 = False; 
http_proxy_ib3 = 'http://proxy.fcen.uba.ar:8080'; 
https_proxy_ib3 = 'https://proxy.fcen.uba.ar:8080'; 
ftp_proxy_ib3 = 'ftp://proxy.fcen.uba.ar:8080'; 

# Generales
import os 
import time 
import copy 
import math 
from random import shuffle 
from pyensembl import EnsemblRelease 
import pandas as pd 
import numpy as np 

# Analisis de secuencias y genomas
from Bio import Entrez, SeqIO, motifs, Seq 
Entrez.email = 'ekolomenski@gmail.com'; 
Entrez.api_key = '9bccbd4949f226a5930868df36c211e8b408'; 
import biomart 

# Graficos
import matplotlib.pyplot as plt 
import seaborn as sns 

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
Funciones para clasificar sitios de union de NKX2-5 segun si tienen sitios cercanos de otros factores de transcripcion

seleccion_sitios_peaks() usa las funciones mas abajo para conseguir la lista de genes con sitios cercanos a NKX2-5

seleccion_genes_no_cercanos() usa estructura similar a seleccion_sitios_peaks() para conseguir la lista de genes sin otros sitios cercanos a NKX2-5

* Hay que correr ShinyGO con los outputs 
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
path_in_main = path_trabajo_main; 
path_out_main = path_output_dump_main + 'GeneOntology\\'; 
# Path que dependen de la especie
path_pwm_human = path_dropbox_main + 'Doctorado\\3-Genes transactivados rio abajo\\1-Interpretacion de datos ChIP-seq\\PWM_human\\'; 
path_pwm_mouse = path_dropbox_main + 'Doctorado\\3-Genes transactivados rio abajo\\1-Interpretacion de datos ChIP-seq\\PWM_mouse\\'; 

# Nombre de archivo con sitios nkx25
nom_sitios_main = '17-su_human_500k_rnaseq'; 
# Nombre de archivo con peaks de ChIP-seq
nom_peaks_main = '13-Anderson100k'; 
# Nombre de archivo que guarda los sitios seleccionados y los cercanos
nom_su_cercanos_main = '18-sitios_nkx25_y_cercanos'; 


#################################### FUNCIONES ####################################


def _main():
    # Funcion para probar funciones en ejecucion del archivo

    '''# Corro seleccion_sitios_peaks()
    M_seleccion = seleccion_sitios_peaks(nom_sitios_main, nom_peaks_main, path_sitios=path_in_main, path_peaks=path_in_main); 
    # Defino la matriz de headers
    M_headers = ['chr_n', 'pos_ini', 'pos_end', 'seq', 'gene', 'log_fc'] + ['tf', 'pos_ini_tf', 'pos_end_tf', 'dist']*9; 
    # Guardo M_seleccion
    M_seleccion = guardar_csv(M_seleccion, nom_su_cercanos_main, path_out=path_out_main, sep_out=';', L_head=M_headers); 
    '''

    '''# Corro seleccion_genes_grupos() para cada uno de los factores de transcripcion
    M_tf = [['nkx25', 9], ['tbx20', 13], ['meis1', 17], ['tgif1', 21], ['hand1', 25], 
            ['maf', 29], ['gata1', 33], ['gata6', 37], ['gata4', 41]]; 
    # Recorro cada uno de los factores de transcripcion
    for curr_tf in M_tf:
        # Corro seleccion_genes_grupos() usando curr_tf[1] para col_dist_tf
        dict_tf = seleccion_genes_grupos(nom_su_cercanos_main, path_arch=path_in_main, col_dist_tf=curr_tf[1]); 
        # Guardo dict_tf con nombre dependiendo de curr_tf[0]
        dict_tf = guardar_dict(dict_tf, 'genes_nkx25_'+curr_tf[0], path_out=path_out_main); 
    '''
    
    '''# Defino L_col_dist_tf_main para correr seleccion_genes_no_cercanos()
    L_col_dist_tf_main = [13, 17, 21, 25, 29, 33, 37, 41]; 
    # Corro seleccion_genes_no_cercanos()
    dict_no_tf = seleccion_genes_no_cercanos(nom_su_cercanos_main, path_arch=path_in_main, L_col_dist_tf=L_col_dist_tf_main); 
    # Guardo dict_no_tf con nombre dependiendo de curr_tf[0]
    dict_no_tf = guardar_dict(dict_no_tf, 'genes_nkx25_solo', path_out=path_out_main); 
    '''

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


def guardar_csv(M_csv, nom_out, path_out='.\\', sep_out=';', L_head=[]):
    # Funcion para guardar una matriz en formato csv 

    # Defino extension
    ext = '.csv'
    # Defino la direccion del archivo en base a dir_out y nom_out
    dirarch = os.path.join(path_out, nom_out + ext); 
    # Creo el archivo
    with open(dirarch, 'w') as F_out:
        print('Archivo ' + nom_out + ext + ' creado.')
    # Lo vuelvo a abrir en modo append
    with open(dirarch, 'a') as F_out:
        # Creo el titulo en base a L_head
        if len(L_head) > 0:
            str_head = ''; 
            # Recorro L_head
            for h in L_head:
                str_head = str_head + str(h) + sep_out; 
            # Elimino la ultima ocurrencia de sep
            str_head=str_head.rstrip(sep_out); 
            # Agrego end of line y guardo en F_out
            F_out.write(str_head + '\n'); 
        # Recorro M_out
        for L_out in M_csv:
            curr_str = ''; 
            # Recorro L_out
            for i in L_out:
                curr_str = curr_str + str(i) + sep_out; 
            # Elimino la ultima ocurrencia de sep
            curr_str=curr_str.rstrip(sep_out); 
            # Agrego end of line y guardo en F_out
            F_out.write(curr_str + '\n'); 
    return M_csv


def guardar_dict(dict_out, nom_out, path_out='.\\', sep_out=';'):
    # Funcion para guardar un diccionario en formato csv 

    # Defino extension
    ext = '.csv'; 
    # Defino la direccion del archivo en base a dir_out y nom_out
    dirarch = os.path.join(path_out, nom_out + ext); 
    # Defino lista de keys
    L_keys = list(dict_out.keys()); 
    # Creo el archivo
    with open(dirarch, 'w') as F_out:
        print('Archivo ' + nom_out + ext + ' creado.')
    # Lo vuelvo a abrir en modo append
    with open(dirarch, 'a') as F_out:
        # Recorro L_keys
        for curr_key in L_keys:
            # Agrego curr_key y dict_out[curr_key] a F_out
            F_out.write(str(curr_key) + sep_out + str(dict_out[curr_key]) + '\n'); 
    return dict_out


def seleccion_genes_no_cercanos(nom_arch, path_arch='', col_gene=4, col_fc=5, L_col_dist_tf=[]):
    # Funcion para devolver listas de genes sin sitios que no sean nkx2-5 cerca

    # Abro el archivo nom_arch en path_arch
    M_arch = abrir_csv(nom_arch, path_arch=path_arch, ignore_headers=True); 
    # Inicializo la lista que se devuelve
    dict_out = {}
    # Recorro M_arch
    for i in range(len(M_arch)):
        curr_row = M_arch[i]; 
        # Defino los valores que me interesan de acuerdo a columnas dadas a la funcion
        curr_gene = curr_row[col_gene]; 
        curr_updown = fc_to_updown(float(curr_row[col_fc])); 
        # Inicializo un booleano para determinar que no haya ningun sitio cerca
        sin_sitio_cerca = True; 
        # Recorro L_col_dist_tf 
        for curr_col_dist_tf in L_col_dist_tf:
            # Veo que ninguna de las columnas en L_col_dist_tf sea 'NA'
            if curr_row[curr_col_dist_tf]!='NA':
                try:
                    # Trato de pasar la distancia a int() antes de confirmar que hay sitio cerca
                    curr_dist = int(curr_row[curr_col_dist_tf]); 
                    sin_sitio_cerca = False; 
                except:
                    print('WARNING: curr_row[curr_col_dist_tf]=' + str(curr_row[curr_col_dist_tf]) + ' da error sin ser igual a NA.')
        # Una vez recorrido L_col_dist_tf, veo si sin_sitio_cerca sigue siendo True
        if sin_sitio_cerca:
            # Reviso que curr_gene no este en dict_out.keys()
            if not (curr_gene in dict_out.keys()):
                # Cargo el valor de updown al diccionario con el gen como key
                dict_out[curr_gene] = int(curr_updown); 
            # Reviso que dict_out[curr_gene] tenga el mismo valor para todas las instancias de curr_gene
            elif (dict_out[curr_gene]!=0) and (dict_out[curr_gene]!=int(curr_updown)):
                print('WARNING: Gen ' + str(curr_gene) + ' tiene mas de un valor de updown. Se carga 0.')
                # Pongo valor 0 y no vuelvo a tirar el warning para curr_gene
                dict_out[curr_gene] = 0; 
    return dict_out


def seleccion_genes_grupos(nom_arch, path_arch='', col_gene=4, col_fc=5, col_dist_tf=9):
    # Funcion para devolver listas de genes cerca de sitios de nkx2-5 con otros sitios cerca

    # Abro el archivo nom_arch en path_arch
    M_arch = abrir_csv(nom_arch, path_arch=path_arch, ignore_headers=True); 
    # Inicializo la lista que se devuelve
    dict_out = {}
    # Recorro M_arch
    for i in range(len(M_arch)):
        curr_row = M_arch[i]; 
        # Defino los valores que me interesan de acuerdo a columnas dadas a la funcion
        curr_gene = curr_row[col_gene]; 
        curr_updown = fc_to_updown(float(curr_row[col_fc])); 
        curr_dist_tf = curr_row[col_dist_tf]; 
        # Cargo solo los genes donde curr_dist_tf es distinto de 'NA'
        if curr_dist_tf!='NA':
            # Si curr_dist_tf no es un numero (si algo falla), esto tira error
            curr_dist_tf = int(curr_row[col_dist_tf]); 
            # Reviso que curr_gene no este en dict_out.keys()
            if not (curr_gene in dict_out.keys()):
                # Cargo el valor de updown al diccionario con el gen como key
                dict_out[curr_gene] = int(curr_updown); 
            # Reviso que dict_out[curr_gene] tenga el mismo valor para todas las instancias de curr_gene
            elif (dict_out[curr_gene]!=0) and (dict_out[curr_gene]!=int(curr_updown)):
                print('WARNING: Gen ' + str(curr_gene) + ' tiene mas de un valor de updown. Se carga 0.')
                # Pongo valor 0 y no vuelvo a tirar el warning para curr_gene
                dict_out[curr_gene] = 0; 
    return dict_out


def fc_to_updown(n):
    # Funcion para convertir fold change (float de -inf a +inf, sin pasar por 0) a updown (-1 o 1)
    if n<0:
        ret = -1; 
    elif n>0:
        ret = 1; 
    else:
        print('ERROR: n es igual a 0.')
        ret = 0; 
    return ret


def seleccion_sitios_peaks(nom_sitios, nom_peaks, path_sitios='', path_peaks='', dist_max=1000):
    # Funcion que agrega info de sitios cercanos a sitios de union en archivo nom_sitios

    # Inicializo la matriz que se devuelve
    M_out = []; 
    # Abro el archivo con sitios de union y extraigo los datos 
    M_sitios = abrir_csv(nom_sitios, path_arch=path_sitios, ignore_headers=True); 
    # Abro el archivo con peaks y extraigo los datos 
    M_peaks = abrir_csv(nom_peaks, path_arch=path_peaks, ignore_headers=True); 
    # Recorro M_sitios
    for i in range(len(M_sitios)):
        curr_sitio = M_sitios[i]; 
        # Veo que tenga largo mas de 4 y gen cercano
        if (len(curr_sitio) > 4) and (curr_sitio[4] != ''):
            # Uso buscar_sitios_cercanos() para buscar curr_sitio en M_peaks
            mult_out = buscar_sitios_cercanos(curr_sitio, M_peaks, dist_max=dist_max); 
            # Recorro mult_out y agrego cada fila a M_out
            for o in range(len(mult_out)):
                curr_out = mult_out[o]; 
                M_out.append(curr_out[:]); 
        elif len(curr_sitio) > 4:
            print('WARNING: sitio ' + str(curr_sitio) + ' de largo mayor a 4 y sin gen cercano.')
    return M_out


def buscar_sitios_cercanos(L_sitio, M_peaks, dist_max=1000):
    # Funcion que busca sitios cercanos a L_sitios en M_peaks
    # Busca sitios a una distancia maxima de dist_max

    # Matriz de columnas importantes por sitio
    M_otros_tf = [['nkx25', 7, 23], ['tbx20', 8, 25], ['meis1', 9, 27], ['tgif1', 10, 29], ['hand1', 11, 31], 
                  ['maf', 12, 33], ['gata1', 13, 35], ['gata6', 14, 37], ['gata4', 15, 39]]; 
    # Inicializo lista con datos de otros sitios de union
    L_out = L_sitio[:]; 
    # Inicializo una matriz por si se devuelve mas de una lista
    M_out = []; 
    # Defino chr_n, pos_ini y pos_end de sitio de union
    chr_n_su = L_sitio[0]; 
    pos_ini_su = int(L_sitio[1]); 
    pos_end_su = int(L_sitio[2]); 
    # Booleano para cuando ya se encontro L_sitio en M_peaks
    sitio_encontrado = False; 
    # Recorro M_peaks 
    for i in range(len(M_peaks)):
        curr_peak = M_peaks[i]; 
        # Defino chr_n, pos_ini y pos_end de peak
        chr_n_peak = curr_peak[0]; 
        pos_ini_peak = int(curr_peak[2]); 
        pos_end_peak = int(curr_peak[3]); 
        # Reviso que chr_n sea igual en sitio y peak
        if chr_n_peak==chr_n_su:
            # Veo si curr_peak contiene a L_sitio
            if ((pos_ini_su<pos_end_peak) and (pos_ini_su>pos_ini_peak)) or ((pos_end_su<pos_end_peak) and (pos_end_su>pos_ini_peak)):
                # Veo si ya se habia encontrado L_sitio
                if sitio_encontrado:
                    print('WARNING: Sitio ' + str(L_sitio[:3]) + ' ya encontrado. Se agrega otra fila.')
                # Uso extraer_sitios_cercanos() para buscar sitios cercanos con M_otros_tf y pasarlos a formato
                L_out = extraer_sitios_cercanos(L_sitio, curr_peak, M_otros_tf, dist_max=dist_max); 
                # Agrego L_out a M_out
                M_out.append(L_out[:]); 
    return M_out


def extraer_sitios_cercanos(L_sitio, L_peak, M_otros_tf, dist_max=1000, sep_su=',', sep_pos='_'):
    # Funcion para buscar sitios cercanos a L_sitios en L_peak
    # Usa M_otros_tf para ver las columnas con info importante para otros sitios
    # Selecciona sitios a distancia menor o igual a dist_max

    # Inicializo la lista que se devuelve
    L_out = L_sitio[:]; 
    # Recorro M_otros_tf
    for i in range(len(M_otros_tf)):
        L_curr_tf = M_otros_tf[i]; 
        # Para cada TF, la lista contiene el nombre del tf, la columna con X y la columna con los sitios
        col_nom = L_curr_tf[0]; 
        col_x = L_curr_tf[1]; 
        col_su = L_curr_tf[2]; 
        # Inicializo lista que se agrega a L_out
        curr_su_out = [str(col_nom), 'NA', 'NA', 'NA']; 
        # Veo si col_x tiene X y si col_su contiene texto
        if L_peak[col_x] == 'X' and L_peak[col_su] != '':
            # Extraigo los sitios de union, separandolos por sep_su
            L_curr_su = L_peak[col_su].split(sep_su); 
            # Inicializo distancia minima a sitio
            min_dist = dist_max+1; 
            # Recorro los sitios de union
            for j in range(len(L_curr_su)):
                # Defino curr_su
                curr_su = L_curr_su[j].split(sep_pos); 
                # Defino los rangos de sitio y curr_su
                r_sitio = [min(int(L_sitio[1]), int(L_sitio[2])), max(int(L_sitio[1]), int(L_sitio[2]))]; 
                r_curr_su = [min(int(curr_su[0]),int(curr_su[1])), max(int(curr_su[0]),int(curr_su[1]))];  
                # Calculo la distancia entre sitio y curr_su
                dist_curr_su = calcular_dist_rangos(r_sitio, r_curr_su); 
                # Veo si curr_su esta a distancia menor a min_dist
                if dist_curr_su < min_dist:
                    # Reviso que dist_curr_su no sea 0 si estoy viendo nkx2-5
                    if not ((dist_curr_su==0) and (str(col_nom)=='nkx25')):
                        # Reviso que dist_curr_su no sea menor a 0
                        if dist_curr_su<0:
                            print('ERROR: Sitio encontrado con valor de distancia negativo')
                        # Redefino curr_su_out[1], curr_su_out[2] y curr_su_out[3]
                        curr_su_out[1] = r_curr_su[0]; 
                        curr_su_out[2] = r_curr_su[1]; 
                        curr_su_out[3] = int(dist_curr_su); 
                        # Actualizo min_dist
                        min_dist = int(dist_curr_su); 
        # Confirmo que no haya problemas con uno de los dos no alineando con el otro
        elif L_peak[col_x] == 'X' or L_peak[col_su] != '':
            print('WARNING: La columna de confirmacion no coincide con la columna de sitio de union.')
        # Una vez confirmo que no hay sitios de union o los recorro todos, agrego curr_su_out al final de L_out
        L_out = L_out + curr_su_out; 
    return L_out


def calcular_dist_rangos(r1, r2):
    # Funcion que devuelve la distancia entre dos rangos
    # Si tienen overlap, devuelve 0
    # Asumo que r1 y r2 tienen la posicion inicial menor a la final

    # Inicializo el valor que se devuelve
    ret = 0; 
    # Asigno el valor mas grande a x y el mas chico a y
    x, y = sorted((r1, r2)); 
    # Si x[1] esta entre x[0] y y[0], entonces los rangos no tienen overlap
    # all() revisa que todos los rangos esten ordenados bien
    if x[0] <= x[1] < y[0] and all(y[0] <= y[1] for y in (r1,r2)):
        ret = y[0] - x[1]; 
    return ret


#################################### RESULTADOS ###################################

output_dump = []; 

if __name__=='__main__':
    output_dump.append(_main()); 

