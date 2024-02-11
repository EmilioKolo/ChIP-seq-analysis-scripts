
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
from Bio import Entrez, SeqIO
Entrez.email = 'ekolomenski@gmail.com'; 
Entrez.api_key = '9bccbd4949f226a5930868df36c211e8b408'; 
import biomart

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
Funcion para contar cuantos genes identificados por GREAT aparecen modificados en RNA-seq

FALTA:
-Abrir GREAT para agarrar lista de genes
-Abrir RNA-seq para agarrar lista de IDs
-Probar comparar (naive) listas de nombres
-Agregar mas complejidad de ser necesario
'''

#################################### VARIABLES ####################################

# Direcciones de archivos y outputs
path_bed_casa = 'D:\\Users\\Admin\\Dropbox\\Doctorado\\3-Genes transactivados rio abajo\\0-Fuentes\\Papers ChIP-seq\\'; 
path_bed_ib3 = 'C:\\Users\\emili\\Dropbox\\Doctorado\\3-Genes transactivados rio abajo\\0-Fuentes\\Papers ChIP-seq\\'; 
path_fasta_casa = 'D:\\Archivos doctorado\\Genomas\\'; 
path_fasta_ib3 = 'X:\\Genomas\\'; 
path_out_casa = 'D:\\Archivos doctorado\\Output_dump\\PeaksClasificados\\'; 
path_out_ib3 = 'X:\\Output_dump\\PeaksClasificados\\'; 
path_rnaseq_casa = 'D:\\Users\\Admin\\Dropbox\\Doctorado\\3-Genes transactivados rio abajo\\1-Interpretacion de datos ChIP-seq\\'; 
path_rnaseq_ib3 = 'C:\\Users\\emili\\Dropbox\\Doctorado\\3-Genes transactivados rio abajo\\1-Interpretacion de datos ChIP-seq\\'; 

# Nombres archivos usados
nombre_GREAT = '3312genesGREATdupays'; 
nombre_rnaseq = '6-Dupays2015_RNAseq'; 


# Variables main()

if ib3:
    path_bed_main = path_bed_ib3; 
    path_fasta_main = path_fasta_ib3; 
    path_out_main = path_out_ib3; 
    path_rnaseq_main = path_rnaseq_ib3; 
else:
    path_bed_main = path_bed_casa; 
    path_fasta_main = path_fasta_casa; 
    path_out_main = path_out_casa; 
    path_rnaseq_main = path_rnaseq_casa; 



#################################### FUNCIONES ####################################


def abrir_csv(nom_arch, dir_arch='', ext='.csv', sep=';', ignore_headers=True):
    # Abre un archivo en formato de matriz con filas por lineas y columnas separadas por sep
    # Devuelve una matriz con una lista de filas, cada una con una lista de columnas
    # ignore_headers ignora la primer fila

    # Inicializo la matriz que se devuelve
    M_out = []; 
    # Defino la direccion del archivo con nom_arch y dir_arch
    if dir_arch=='':
        filepath = nom_arch + ext; 
    else:
        filepath = os.path.join(dir_arch, nom_arch + ext); 
    # Booleano para revisar que se ignore la primera fila
    header_ignorado = not ignore_headers; 
    # Abro el archivo en modo reading
    with open(filepath, 'r') as F:
        print('Archivo ' + nom_arch + ' abierto.')
        # Recorro cada fila
        for curr_line in F.readlines():
            if header_ignorado:
                # Transformo la linea en lista
                L_line = curr_line.rstrip().split(sep=sep); 
                # Cargo L_line en M_out
                M_out.append(L_line[:]); 
            else:
                # Si header_ignorado empieza sin ser True, se ignora la primera linea
                header_ignorado = True; 
    return M_out


def abrir_csv_headers(nom_arch, dir_arch='', ext='.csv', sep=';'):
    # Abre un archivo en formato de matriz con filas por lineas y columnas separadas por sep
    # Devuelve una lista con los headers (primera fila) y una matriz con una lista de filas, cada una con una lista de columnas
    # Usa abrir_csv() con ignore_headers=False para generar la matriz

    # Cargo M_csv con abrir_csv()
    M_csv = abrir_csv(nom_arch, dir_arch=dir_arch, ext=ext, sep=sep, ignore_headers=True); 
    # Inicializo la lista de headers
    L_head = []; 
    # Leo la primer fila del csv para extraer los headers
    # Defino la direccion del archivo con nom_arch y dir_arch
    if dir_arch=='':
        filepath = nom_arch + ext; 
    else:
        filepath = os.path.join(dir_arch, nom_arch + ext); 
    # Abro el archivo en modo reading
    with open(filepath, 'r') as F:
        str_head = F.readlines()[0]; 
        L_head = str_head.rstrip().split(sep=sep); 
    return M_csv, L_head


def abrir_lista_genes_rnaseq(nom_rnaseq, path_rnaseq='', incluir_aliases=True): 
    # Devuelve una lista de genes de una tabla .csv con resultados de RNA-seq

    # Abro el archivo RNA-seq
    M_rnaseq, L_headers_rnaseq = abrir_csv_headers(nom_rnaseq, dir_arch=path_rnaseq); 
    # Agarro los ids de las columnas que me interesan
    gene_symbol_col = L_headers_rnaseq.index('Gene Symbol'); 
    aliases_col = L_headers_rnaseq.index('Aliases'); 
    # Inicializo la lista de genes que se devuelven
    L_genes = []; 
    # Recorro M_rnaseq
    for i in range(len(M_rnaseq)):
        curr_L = M_rnaseq[i]; 
        # Extraigo la info de las columnas que me interesan
        curr_gene_symbol = curr_L[gene_symbol_col]; 
        curr_aliases = curr_L[aliases_col].split('|'); 
        # Agrego los genes a L_genes
        L_genes.append(str(curr_gene_symbol)); 
        # Agrego aliases si incluir_aliases es True
        if incluir_aliases:
            for alias in curr_aliases:
                if alias != '':
                    L_genes.append(str(alias)); 
    return L_genes


def abrir_txt(nom_arch, path_arch=''):
    # Abre un txt y lo devuelve como una lista de filas

    # Inicializo la lista que se devuelve
    L_out = []; 
    # Veo si nom_arch termina en .txt
    if nom_arch[-4:]=='.txt':
        nom_usado = nom_arch; 
    else:
        nom_usado = nom_arch + '.txt'; 
    # Defino la direccion del archivo con nom_arch y dir_arch
    if path_arch=='':
        filepath = nom_usado; 
    else:
        filepath = os.path.join(path_arch, nom_usado); 
    # Abro el archivo en modo reading
    with open(filepath, 'r') as F:
        print('Archivo ' + nom_arch + ' abierto.')
        # Recorro cada fila
        for curr_line in F.readlines():
            L_out.append(curr_line.rstrip()); 
    return L_out


def buscar_comunes_listas(lista1, lista2):
    # Busca los elementos en comun entre dos listas

    # Inicializo la lista que se devuelve
    L_out = []; 
    # Recorro lista1
    for elemento1 in lista1:
        # Si un elemento de lista1 esta en lista2, lo agrego a L_out
        if elemento1 in lista2:
            # Veo que no este repetido en L_out
            if elemento1 in L_out:
                print(elemento1 + ' repetido en lista1')
            else:
                L_out.append(str(elemento1)); 
    return L_out


def dict_regulation_rnaseq(nom_rnaseq, path_rnaseq=''): 
    # Devuelve un diccionario con los valores de up/down regulacion para los aliases de resultados de RNA-seq

    # Inicializo el diccionario que se devuelve
    dict_regulation = {}; 
    # Abro el archivo RNA-seq
    M_rnaseq, L_headers_rnaseq = abrir_csv_headers(nom_rnaseq, dir_arch=path_rnaseq); 
    # Agarro los ids de las columnas que me interesan
    gene_symbol_col = L_headers_rnaseq.index('Gene Symbol'); 
    aliases_col = L_headers_rnaseq.index('Aliases'); 
    regulation_col = L_headers_rnaseq.index('Regulation ([hypo] vs [co])'); 
    # Inicializo la lista de genes que se devuelven
    L_genes = []; 
    # Recorro M_rnaseq
    for i in range(len(M_rnaseq)):
        curr_L = M_rnaseq[i]; 
        # Extraigo la info de las columnas que me interesan
        curr_gene_symbol = curr_L[gene_symbol_col]; 
        curr_aliases = curr_L[aliases_col].split('|'); 
        curr_regulation = curr_L[regulation_col]; 
        # Agrego curr_regulation a dict_regulation
        if curr_gene_symbol in dict_regulation:
            print('ERROR: ' + curr_gene_symbol + ' repetido.')
            if dict_regulation[curr_gene_symbol] == curr_regulation:
                print('\tARREGLADO')
        else:
            dict_regulation[curr_gene_symbol] = curr_regulation; 
    return dict_regulation


def _main():
    # Funcion para probar funciones en ejecucion del archivo

    # Genero las listas de genes y la lista de combinados
    L_great = abrir_txt(nombre_GREAT, path_arch=path_rnaseq_main); 
    L_genes_rnaseq = abrir_lista_genes_rnaseq(nombre_rnaseq, path_rnaseq=path_rnaseq_main, incluir_aliases=False); 
    L_comun = buscar_comunes_listas(L_genes_rnaseq, L_great); # 309 genes
    # 
    dict_regulation = dict_regulation_rnaseq(nombre_rnaseq, path_rnaseq=path_rnaseq_main); 
    # Inicializo contadores de up/down regulacion
    cont_up = 0; 
    cont_down = 0; 
    # Recorro L_comun
    for gen in L_comun:
        if dict_regulation[gen] == 'down':
            cont_down += 1; 
        elif dict_regulation[gen] == 'up':
            cont_up += 1; 
        else:
            print('ERROR: ' + gen + ' da ' + dict_regulation[gen] + ' en dict_regulation.')
    print('cont_up=' + str(cont_up)) # 137 upregulados
    print('cont_down=' + str(cont_down)) # 172 downregulados
    return L_comun


#################################### RESULTADOS ###################################

output_dump = []; 

if __name__=='__main__':
    output_dump.append(_main()); 

