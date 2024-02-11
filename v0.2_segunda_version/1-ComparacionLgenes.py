
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

# Genomas de referencia
hg19 = EnsemblRelease(75, species='human'); 
hg38 = EnsemblRelease(102, species='human'); 
mm9 = EnsemblRelease(54, species='mouse'); 
mm10 = EnsemblRelease(102, species='mouse'); 


'''
Funciones para analisis de listas de genes

Usadas para comparar listas de genes sacadas de GREAT con listas de genes de RNA-seq (Dupays)
'''

#################################### VARIABLES ####################################

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

nombre_main_test = nombre_great_dupays_1mb; 

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


def comparar_listas(L1, L2):
    # Compara los elementos 

    # Inicializo las listas que se devuelven
    L_com = []; 
    L_diff1 = []; 
    L_diff2 = []; 
    # Recorro L1
    for i in range(len(L1)):
        curr_gene = L1[i]; 
        # Recorro L2
        for j in range(len(L2)):
            comp_gene = L2[j]; 
            # Comparo dos genes
            if curr_gene == comp_gene:
                # Si son iguales, los anoto en L_com
                L_com.append(str(curr_gene)); 
    # Recorro L1 y L2 de nuevo para determinar L_diff1 y L_diff2
    for i in range(len(L1)):
        # Veo si L1[i] esta en L_com
        if not (L1[i] in L_com): 
            # Si no esta en L_com, lo agrego a L_diff1
            L_diff1.append(str(L1[i])); 
    for i in range(len(L2)):
        # Veo si L2[i] esta en L_com
        if not (L2[i] in L_com): 
            # Si no esta en L_com, lo agrego a L_diff2
            L_diff2.append(str(L2[i])); 
    return L_com, L_diff1, L_diff2


def comparar_listas_display(L1, L2, titulo_L1='L1', titulo_L2='L2', print_L1=False, print_L2=False):
    # Usa comparar_listas y printea resultados importantes

    # Display
    print('>Comparando ' + titulo_L1 + ' con ' + titulo_L2 + '.')
    # Comparo listas
    L_com, L_diff1, L_diff2 = comparar_listas(L1, L2); 
    # Display
    print('En comun: ' + str(len(L_com)))
    print('En ' + titulo_L1 + ' no en comun: ' + str(len(L_diff1)))
    if print_L1:
        print(L_diff1)
    print('En ' + titulo_L2 + ' no en comun: ' + str(len(L_diff2)))
    if print_L2:
        print(L_diff2)
    return L_com


def _main():
    # Funcion para probar funciones en ejecucion del archivo

    # Inicializo una matriz para devolver
    M_out = []; 

    # Abro archivos great y extraigo lista de genes
    L_great_dupays_1mb = abrir_archivo_great_genes(nombre_great_dupays_1mb, dir_arch=path_main); 
    L_great_dupays_100kb = abrir_archivo_great_genes(nombre_great_dupays_100kb, dir_arch=path_main); 
    #L_great_anderson_1mb = abrir_archivo_great_genes(nombre_great_anderson_1mb, dir_arch=path_main); 
    #L_great_anderson_100kb = abrir_archivo_great_genes(nombre_great_anderson_100kb, dir_arch=path_main); 
    #shuffle(L_great_dupays_100kb)
    #print(L_great_dupays_100kb[:10])

    # Abro RNA-seq y extraigo lista de genes
    L_rnaseq_dupays = abrir_rnaseq_genes(nombre_rnaseq, gene_name_pos, dir_arch=path_main); 
    #print(L_rnaseq_dupays[:10])

    # Abro las listas de genes que Ceci Ballares proceso de RNA-seq
    L_rnaseq_ceci_5kb_genebody = abrir_txt(nombre_rnaseq_ceci_5kb_genebody, dir_arch=path_main); 
    L_rnaseq_ceci_10kb_genebody = abrir_txt(nombre_rnaseq_ceci_10kb_genebody, dir_arch=path_main); 
    L_rnaseq_ceci_25kb_genebody = abrir_txt(nombre_rnaseq_ceci_25kb_genebody, dir_arch=path_main); 
    L_rnaseq_ceci_10kb_TSS = abrir_txt(nombre_rnaseq_ceci_10kb_TSS, dir_arch=path_main); 

    # Comparo listas
    # Dupays 1mb vs RNA-seq total
    #M_out.append(comparar_listas_display(L_great_dupays_1mb, L_rnaseq_dupays, titulo_L1='Dupays 1mb', titulo_L2='RNA-seq')); 
    # Dupays 100kb vs RNA-seq total
    #M_out.append(comparar_listas_display(L_great_dupays_100kb, L_rnaseq_dupays, titulo_L1='Dupays 100kb', titulo_L2='RNA-seq'));
    print()
    print('### RNA-seq ceci 5kb genebody')
    # Dupays 1mb vs RNA-seq ceci 5kb genebody
    M_out.append(comparar_listas_display(L_great_dupays_1mb, L_rnaseq_ceci_5kb_genebody, titulo_L1='Dupays 1mb', titulo_L2='RNA-seq ceci 5kb genebody', print_L2=True)); 
    # Dupays 100kb vs RNA-seq ceci 5kb genebody
    M_out.append(comparar_listas_display(L_great_dupays_100kb, L_rnaseq_ceci_5kb_genebody, titulo_L1='Dupays 100kb', titulo_L2='RNA-seq ceci 5kb genebody', print_L2=True)); 
    print()
    print('### RNA-seq ceci 10kb genebody')
    # Dupays 1mb vs RNA-seq ceci 10kb genebody
    M_out.append(comparar_listas_display(L_great_dupays_1mb, L_rnaseq_ceci_10kb_genebody, titulo_L1='Dupays 1mb', titulo_L2='RNA-seq ceci 10kb genebody', print_L2=True)); 
    # Dupays 100kb vs RNA-seq ceci 10kb genebody
    M_out.append(comparar_listas_display(L_great_dupays_100kb, L_rnaseq_ceci_10kb_genebody, titulo_L1='Dupays 100kb', titulo_L2='RNA-seq ceci 10kb genebody', print_L2=True)); 
    print()
    print('### RNA-seq ceci 25kb genebody')
    # Dupays 1mb vs RNA-seq ceci 25kb genebody
    M_out.append(comparar_listas_display(L_great_dupays_1mb, L_rnaseq_ceci_25kb_genebody, titulo_L1='Dupays 1mb', titulo_L2='RNA-seq ceci 25kb genebody', print_L2=True)); 
    # Dupays 100kb vs RNA-seq ceci 25kb genebody
    M_out.append(comparar_listas_display(L_great_dupays_100kb, L_rnaseq_ceci_25kb_genebody, titulo_L1='Dupays 100kb', titulo_L2='RNA-seq ceci 25kb genebody', print_L2=True)); 
    print()
    print('### RNA-seq ceci 10kb TSS')
    # Dupays 1mb vs RNA-seq ceci 10kb TSS
    M_out.append(comparar_listas_display(L_great_dupays_1mb, L_rnaseq_ceci_10kb_TSS, titulo_L1='Dupays 1mb', titulo_L2='RNA-seq ceci 10kb TSS', print_L2=True)); 
    # Dupays 100kb vs RNA-seq ceci 10kb TSS
    M_out.append(comparar_listas_display(L_great_dupays_100kb, L_rnaseq_ceci_10kb_TSS, titulo_L1='Dupays 100kb', titulo_L2='RNA-seq ceci 10kb TSS', print_L2=True)); 

    ### Chequear distancias de genes GREAT a peak
    # ej. AACS, Peak_1739 (+211448); Peak_1739: chr12, 125761031, 125761714
    #print(hg19.genes_by_name('AACS'))
    # [Gene(gene_id='ENSG00000081760', gene_name='AACS', biotype='protein_coding', contig='12', start=125549925, end=125627873, strand='+', genome='GRCh37')]
    # [Gene(gene_id='ENSG00000081760', gene_name='AACS', biotype='protein_coding', contig='12', start=125065434, end=125143333, strand='+', genome='GRCh38')]
    # Gen: (125549925, 125627873) // Peak: (125761031, 125761714)
    # ej. NPPA, Peak_36 (-228); Peak_36: chr1, 11907794, 11908342
    #print(hg19.genes_by_name('NPPA'))
    # [Gene(gene_id='ENSG00000175206', gene_name='NPPA', biotype='protein_coding', contig='1', start=11905766, end=11908402, strand='-', genome='GRCh37')]
    # [Gene(gene_id='ENSG00000175206', gene_name='NPPA', biotype='protein_coding', contig='1', start=11845709, end=11848345, strand='-', genome='GRCh38')]
    # Gen: (11905766, 11908402) // Peak: (11907794, 11908342)
    
    return M_out


#################################### RESULTADOS ###################################

output_dump = []; 

if __name__=='__main__':
    output_dump.append(_main()); 

