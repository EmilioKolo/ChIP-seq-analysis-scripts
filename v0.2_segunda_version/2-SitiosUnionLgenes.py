
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

# Importo analisischip
import sys
path_analisischip_casa = 'C:\\Users\\Admin\\Documents\\Scripts\\Git\\AnalisisSecuenciasChIP\\'; 
path_analisischip_ib3 = ''; 
sys.path.append(path_analisischip_casa); 
from analisischip.seq import seq_data, seq_handler

# Genomas de referencia
hg19 = EnsemblRelease(75, species='human'); 
hg38 = EnsemblRelease(102, species='human'); 
mm9 = EnsemblRelease(54, species='mouse'); 
mm10 = EnsemblRelease(102, species='mouse'); 


'''
Funciones para busqueda de sitios de union alrededor de genes

Usadas para buscar sitios de union en listas de genes sacadas de GREAT

Funciones para sacar secuencias en promotores (usadas para MEME-ChIP)

### FALTA
# Funciones para sacar secuencias para MEME-ChIP con seq_data
    X Funcion para conseguir secuencia: seq_data.secuencias_rangos_fasta()
        # Probar seq_data.secuencias_rangos_fasta()
    X Funcion para recortar rangos a tamaño estandar: seq_data.rangos_mismo_largo()
        # Probar seq_data._rango_largo_definido()
        # Probar seq_data.rangos_mismo_largo()
    # Funcion para seleccionar rangos de seq_data que esten cerca de genes
    # Funcion para seleccionar rangos de seq_data que tengan sitios de union
# Funciones para sacar secuencias para MEME-ChIP con GREAT
###
'''

#################################### VARIABLES ####################################

# Defino la direccion del .fasta
path_fasta_casa = 'D:\\Archivos doctorado\\Genomas\\'; 
path_fasta_ib3 = 'X:\\Genomas\\'; 
# Defino las direcciones de output
path_out_casa = 'D:\\Archivos doctorado\\Output_dump\\SitiosUnionLgenes\\'; 
path_out_graficos_casa = 'D:\\Archivos doctorado\\Output_dump\\Graficos\\SitiosUnionLgenes\\'; 
path_out_ib3 = 'X:\\Output_dump\\SitiosUnionLgenes\\'; 
path_out_graficos_ib3 = 'X:\\Output_dump\\Graficos\\SitiosUnionLgenes\\'; 

path_ib3 = 'C:\\Users\\emili\\Dropbox\\Doctorado\\3-Genes transactivados rio abajo\\1-Interpretacion de datos ChIP-seq\\'; 
path_casa = 'D:\\Dropbox\\Doctorado\\3-Genes transactivados rio abajo\\1-Interpretacion de datos ChIP-seq\\'; 

nombre_rnaseq = 'Dupays2015_RNAseq'; 
ext_rnaseq = '.csv'; 
sep_rnaseq = ';'; 
gene_id_pos = 12; 
gene_name_pos = 9; #10 contiene aliases y 11 contiene Entrez ID

nombre_great_dupays_1mb = 'GREAT_Dupays_1Mb'; 
nombre_great_dupays_100kb = 'GREAT_Dupays_100kb'; 
nombre_great_anderson_1mb = 'GREAT_Anderson_1Mb'; 
nombre_great_anderson_100kb = 'GREAT_Anderson_100kb'; 
ext_great = '.txt'; 
sep_great = '\t'; 

nombre_rnaseq_ceci_5kb_genebody = 'rnaseq_ceci_5kb_genebody'; 
nombre_rnaseq_ceci_10kb_genebody = 'rnaseq_ceci_10kb_genebody'; 
nombre_rnaseq_ceci_25kb_genebody = 'rnaseq_ceci_25kb_genebody'; 
nombre_rnaseq_ceci_10kb_TSS = 'rnaseq_ceci_10kb_TSS'; 
ext_rnaseq_ceci = '.txt'; 
sep_rnaseq_ceci = ';'; 

# Variables main()

path_main = path_casa; 
path_fasta_main = path_fasta_casa; 
path_out_main = path_out_casa; 
path_out_graficos_main = path_out_graficos_casa; 


#################################### FUNCIONES ####################################


def abrir_archivo_great_genes(nom_arch, dir_arch='', ext='.txt', sep='\t', ignore_headers=True):
    # Abre un archivo GREAT con lista de genes y devuelve la lista de genes

    # Inicializo la lista que se devuelve
    L_out = []; 
    # Recopilo la matriz del output de GREAT
    M_csv = abrir_csv(nom_arch, dir_arch=dir_arch, ext=ext, sep=sep, ignore_headers=ignore_headers); 
    # Recorro M_csv
    for i in range(len(M_csv)):
        # Solo agarro las filas con elementos (por si hay errores)
        if len(M_csv[i])>0:
            # Agrego los nombres de los genes (en posicion 0) a L_out
            L_out.append(M_csv[i][0]); 
    return L_out


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


def abrir_rnaseq_genes(nom_arch, col_genes, dir_arch='', ext='.csv', sep=';', ignore_headers=True):
    # Abre un archivo .csv con output de RNA-seq y extrae genes de una columna dada por col_genes

    # Inicializo la lista que se devuelve
    L_out = []; 
    # Recopilo la matriz del output de RNA-seq
    M_csv = abrir_csv(nom_arch, dir_arch=dir_arch, ext=ext, sep=sep, ignore_headers=ignore_headers); 
    # Recorro M_csv
    for i in range(len(M_csv)):
        # Solo agarro las filas con suficientes elementos para llegar a col_genes (por si hay errores)
        if len(M_csv[i])>col_genes+1:
            # Agrego los nombres de los genes (en posicion 0) a L_out
            L_out.append(M_csv[i][col_genes]); 
        elif len(M_csv[i])>0:
            print('ERROR. Fila demasiado corta: ' + str(M_csv[i]))
    return L_out


def abrir_txt(nom_arch, dir_arch='', ext='.txt', ignore_header=False):
    # Abre un archivo de texto
    # Devuelve una lista con cada linea, sin el fin de linea

    # Inicializo la lista que se devuelve
    L_out = []; 
    # Defino la direccion del archivo con nom_arch y dir_arch
    if dir_arch=='':
        filepath = nom_arch + ext; 
    else:
        filepath = os.path.join(dir_arch, nom_arch + ext); 
    # Booleano para revisar que se ignore la primera fila
    header_ignorado = not ignore_header; 
    # Abro el archivo en modo reading
    with open(filepath, 'r') as F:
        print('Archivo ' + nom_arch + ' abierto.')
        # Recorro cada fila
        for curr_line in F.readlines():
            if header_ignorado:
                # Elimino el fin de linea
                line = curr_line.rstrip(); 
                # Cargo line en L_out si tiene largo mayor a 0
                if len(line) > 0:
                    L_out.append(str(line)); 
            else:
                # Si header_ignorado empieza sin ser True, se ignora la primera linea
                header_ignorado = True; 
    return L_out


def gene_name_to_id(gene_name, genome):
    # Busca todos los genes del genoma correspondientes a gene_name
    # Los devuelve como una lista
    try:
        L_out = genome.genes_by_name(gene_name); 
    except:
        print('Gen ' + str(gene_name) + ' no encontrado en genoma. Se devuelve lista vacia.')
        L_out = []; 
    return L_out


def _main():
    # Funcion para probar funciones en ejecucion del archivo

    # Inicializo una matriz para devolver
    M_out = []; 

    # Abro archivos great y extraigo lista de genes
    #L_great_dupays_1mb = abrir_archivo_great_genes(nombre_great_dupays_1mb, dir_arch=path_main); 
    #L_great_dupays_100kb = abrir_archivo_great_genes(nombre_great_dupays_100kb, dir_arch=path_main); 
    #L_great_anderson_1mb = abrir_archivo_great_genes(nombre_great_anderson_1mb, dir_arch=path_main); 
    #L_great_anderson_100kb = abrir_archivo_great_genes(nombre_great_anderson_100kb, dir_arch=path_main); 
    #shuffle(L_great_dupays_100kb)
    #print(L_great_dupays_100kb[:10])

    # Abro RNA-seq y extraigo lista de genes
    #L_rnaseq_dupays = abrir_rnaseq_genes(nombre_rnaseq, gene_name_pos, dir_arch=path_main); 
    #print(L_rnaseq_dupays[:10])

    # Abro las listas de genes que Ceci Ballares proceso de RNA-seq
    #L_rnaseq_ceci_5kb_genebody = abrir_txt(nombre_rnaseq_ceci_5kb_genebody, dir_arch=path_main); 
    #L_rnaseq_ceci_10kb_genebody = abrir_txt(nombre_rnaseq_ceci_10kb_genebody, dir_arch=path_main); 
    #L_rnaseq_ceci_25kb_genebody = abrir_txt(nombre_rnaseq_ceci_25kb_genebody, dir_arch=path_main); 
    #L_rnaseq_ceci_10kb_TSS = abrir_txt(nombre_rnaseq_ceci_10kb_TSS, dir_arch=path_main); 

    # Inicializo seq_handler
    handle = seq_handler('hg19', path_fasta=path_fasta_main, path_archivos=path_out_main); 
    '''# Defino genomas a descargar
    #L_genomas_descargados = ['hg19', 'mm9', 'hg38', 'mm10'];
    #L_genomas_descargados = ['hg38', 'mm10']; 
    # Descargo los genomas
    #handle.generar_archivos_base([], L_genomas_descargados, [], descargar_genomas=True, path_fasta=path_fasta_main, path_archivos=path_out_main); '''
    # Defino variables para archivos .bed
    path_bed_casa = 'D:\\Dropbox\\Doctorado\\3-Genes transactivados rio abajo\\0-Fuentes\\Papers ChIP-seq\\';
    path_bed_ib3 = ''; ### BUSCAR PATH
    path_bed = path_bed_casa; 
    bed_dupays = ['Dupays2015', 'mm9']; 
    bed_anderson = ['Anderson2018-GSE89457consensus', 'hg19']; 
    #L_bed = [bed_dupays, bed_anderson]; 
    L_bed = []; 
    # Defino genomas para buscar promotores
    L_genomas_promotores = ['mm9', 'hg19']; 
    # Defino rangos usados para promotores
    L_rangos = [(-1500,1500), (-10000,10000), (-50000,50000)]; 
    # Creo archivos promotores y .bed
    #handle.generar_archivos_base(L_rangos, L_genomas_promotores, L_bed, path_fasta=path_fasta_main, path_archivos=path_out_main, path_bed=path_bed); 
    
    return M_out


#################################### RESULTADOS ###################################

output_dump = []; 

if __name__=='__main__':
    output_dump.append(_main()); 
