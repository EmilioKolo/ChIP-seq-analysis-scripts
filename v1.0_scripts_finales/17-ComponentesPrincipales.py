
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

# Analisis de secuencias y genomas
from Bio import Entrez, SeqIO, motifs, Seq 
from numpy import array 
from Bio.Cluster import pca 
from sklearn.preprocessing import StandardScaler 
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
Scripts para hacer analisis de componentes principales y/o preparar datos para ese analisis

pipeline_graficos_finales() tiene hardcodeada la generacion de varios graficos importantes
    * Usa generar_grafico_2d() y generar_grafico_3d()

pipeline_generar_lista_sitios() genera tablas con sitios de union de NKX2-5, sus secuencias y todos los genes cercanos

pipeline_agregar_info_rnaseq() agrega filtrado por RNA-seq de genes cercanos (9-ClasificacionPeaks3.py genera uno de los archivos que necesita)

pipeline_pca() termina de procesar las filas del csv generado por pipeline_agregar_info_rnaseq() y genera tablas y graficos relevantes

generar_pca_secuencias() hace PCA en base a solo informacion de secuencia

generar_pca_seq_one_hot() genera una tabla con todos los datos necesarios para (???)

### FALTA
- Mejorar eficiencia de generacion de graficos
    * Ver generar_grafico_2d() y generar_grafico_3d()
    - Separar datos por label, crear figura recorriendo labels en vez de cada punto
###
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
path_out_main = path_output_dump_main + 'PCA\\'; 
# Path que dependen de la especie
path_pwm_human = path_dropbox_main + 'Doctorado\\3-Genes transactivados rio abajo\\1-Interpretacion de datos ChIP-seq\\PWM_human\\'; 
path_pwm_mouse = path_dropbox_main + 'Doctorado\\3-Genes transactivados rio abajo\\1-Interpretacion de datos ChIP-seq\\PWM_mouse\\'; 


nom_arch_mouse = '13-Dupays100k'; 
genoma_mouse = mm9; 
nom_genoma_mouse = 'mm9'; 

nom_arch_human = '13-Anderson100k'; 
genoma_human = hg19; 
nom_genoma_human = 'hg19'; 

nom_su_arch_mouse = '17-su_mouse_500k'; 
nom_su_arch_human = '17-su_human_500k'; 

nom_rnaseq_mouse = '17-lista_genes_rnaseq_raton'; 
nom_rnaseq_human = '17-lista_genes_rnaseq_humano'; 

nom_out_mouse = '17-su_mouse_500k_rnaseq'; 
nom_out_human = '17-su_human_500k_rnaseq'; 

lista_col = [4]; 



#################################### FUNCIONES ####################################


def _main():
    # Funcion para probar funciones en ejecucion del archivo

    '''# Funciones para generar archivos .csv con sitios de union de NKX2-5 usando pipeline_generar_lista_sitios()
    M_sitios_mouse = pipeline_generar_lista_sitios(nom_arch_mouse, genoma_mouse, nom_genoma_mouse, nom_su_arch_mouse, path_csv=path_in_main, 
                                                   path_fasta=path_fasta_main, path_out=path_out_main, sep_csv=';', sep_out=';'); 
    shuffle(M_sitios_mouse); 
    for s in M_sitios_mouse[:20]:
        print(s)
    M_sitios_human = pipeline_generar_lista_sitios(nom_arch_human, genoma_human, nom_genoma_human, nom_su_arch_human, path_csv=path_in_main, 
                                                   path_fasta=path_fasta_main, path_out=path_out_main, sep_csv=';', sep_out=';'); 
    shuffle(M_sitios_human); 
    for s in M_sitios_human[:20]:
        print(s)
    '''
    
    '''# Funciones para agregar confirmacion por RNA-seq usando pipeline_agregar_info_rnaseq()
    M_sitios_mouse = pipeline_agregar_info_rnaseq(nom_su_arch_mouse, nom_rnaseq_mouse, nom_out_mouse, path_su=path_in_main, path_rnaseq=path_in_main, 
                                                  path_out=path_in_main, col_listas=lista_col, col_genes=4); 
    M_sitios_human = pipeline_agregar_info_rnaseq(nom_su_arch_human, nom_rnaseq_human, nom_out_human, path_su=path_in_main, path_rnaseq=path_in_main, 
                                                  path_out=path_in_main, col_listas=lista_col, col_genes=4); 
    '''

    '''# Pruebas de PCA
    matrix = array([[ 0.,  0.,  0.],
                    [ 1.,  0.,  0.],
                    [ 7.,  3.,  0.],
                    [ 4.,  2.,  6.]])
    col_mean, coord, pc, ev = pca(matrix)
    m = matrix - (col_mean + dot(coord, pc))
    #print('m: ')
    #print(m)
    #print('col means: ')
    #print(col_mean)
    #print('coord: ')
    #print(coord)
    #print('coord n=2:')
    #print(coord[:,:2])
    print('pc: ')
    print(pc)
    print('eigenvalues:')
    print(ev)
    print('sum(ev): ')
    print(sum(ev))
    print('var_exp: ')
    print([x/sum(ev) for x in ev])'''

    # Pruebas de invertir_seq()
    #print(invertir_seq('AAGTG'))

    '''# Pruebas de procesar_seq()
    M_csv_mouse = abrir_csv(nom_out_mouse, path_arch=path_in_main, ignore_headers=True); 
    M_csv_human = abrir_csv(nom_out_human, path_arch=path_in_main, ignore_headers=True); 
    shuffle(M_csv_mouse); 
    shuffle(M_csv_human); 
    for i in M_csv_human[:20]:
        print(i)
        print(procesar_seq(i[3]))
    print()
    print('#####################################################################')
    print()
    for j in M_csv_mouse[:20]:
        print(j)
        print(procesar_seq(j[3]))'''
    
    '''# Pruebas de pipeline_pca()
    M_pca_human, pc_var_human, pc_human = pipeline_pca(nom_out_human, hg19, 'pca_human_dist_inv', path_out=path_out_main, path_arch=path_in_main, 
                                                       n_pca=2, curr_dist_inv=True, curr_dist_rangos=False); 
    display_pca(M_pca_human, pc_var_human, pc_human, titulo='human', n_pc=5, n_ej=5); 
    M_pca_mouse, pc_var_mouse, pc_mouse = pipeline_pca(nom_out_mouse, mm9, 'pca_mouse_dist_inv', path_out=path_out_main, path_arch=path_in_main, 
                                                       n_pca=2, curr_dist_inv=True, curr_dist_rangos=False); 
    display_pca(M_pca_mouse, pc_var_mouse, pc_mouse, titulo='mouse', n_pc=5, n_ej=5); 
    '''
    
    # Prueba de generar_pca_secuencias()
    M_pca_seq_human, pc_var_seq_human, pc_seq_human = generar_pca_secuencias(nom_out_human, 'pca_seq_n3_human', path_out=path_out_main, path_arch=path_in_main, 
                                                                             sep_in=';', ext_in='.csv', n_pca=3); 
    display_pca(M_pca_seq_human, pc_var_seq_human, pc_seq_human, titulo='seq human', n_pc=5, n_ej=5); 
    '''M_pca_seq_mouse, pc_var_seq_mouse, pc_seq_mouse = generar_pca_secuencias(nom_out_mouse, 'pca_seq_n3_mouse', path_out=path_out_main, path_arch=path_in_main, 
                                                                             sep_in=';', ext_in='.csv', n_pca=3); 
    display_pca(M_pca_seq_mouse, pc_var_seq_mouse, pc_seq_mouse, titulo='seq mouse', n_pc=5, n_ej=5); 
    '''

    # Pruebas de std_len_M_pca()
    #M_test = [[1,2,3], ['a'], [2,4,5,6]]; 
    #print(std_len_M_pca(M_test))

    '''# Pruebas de generar_pca_seq_one_hot()
    generar_pca_seq_one_hot(nom_out_human, 'procesado_pca_seq_n3_human_dist_norm', hg19, path_out=path_out_main, path_arch=path_in_main, 
                            sep_in=';', ext_in='.csv', n_pca=3, dist_t='n', onehot_full_seq=True); 
    generar_pca_seq_one_hot(nom_out_mouse, 'procesado_pca_seq_n3_mouse_dist_norm', mm9, path_out=path_out_main, path_arch=path_in_main, 
                            sep_in=';', ext_in='.csv', n_pca=3, dist_t='n', onehot_full_seq=True); 

    generar_pca_seq_one_hot(nom_out_human, 'procesado_pca_seq_n3_human_onehot_dist_norm', hg19, path_out=path_out_main, path_arch=path_in_main, 
                            sep_in=';', ext_in='.csv', n_pca=3, dist_t='n', onehot_full_seq=False); 
    generar_pca_seq_one_hot(nom_out_mouse, 'procesado_pca_seq_n3_mouse_onehot_dist_norm', mm9, path_out=path_out_main, path_arch=path_in_main, 
                            sep_in=';', ext_in='.csv', n_pca=3, dist_t='n', onehot_full_seq=False); 
    '''

    # Pruebas de generacion de graficos con pipeline_graficos_finales()
    #pipeline_graficos_finales('datos_pca_full_human', nom_mod='human', path_full=path_in_main, path_out=path_out_main); 
    #pipeline_graficos_finales('datos_pca_full_mouse', nom_mod='mouse', path_full=path_in_main, path_out=path_out_main); 

    # Uso pipeline_pca() para generar 
    #M_pca_human, pc_var_human, pc_human = pipeline_pca(nom_out_human, hg19, 'pca_human_nodist', path_out=path_out_main, path_arch=path_in_main, 
    #                                                   n_pca=2, dist_t='0'); 
    #display_pca(M_pca_human, pc_var_human, pc_human, titulo='human_congreso', n_pc=5, n_ej=5); 

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


def abrir_csv_headers(nom_arch, dir_arch='', ext='.csv', sep=';'):
    # Abre un archivo en formato de matriz con filas por lineas y columnas separadas por sep
    # Devuelve una lista con los elementos de la primera fila

    # Inicializo la lista que se devuelve
    L_out = []; 
    # Defino la direccion del archivo con nom_arch y dir_arch
    if dir_arch=='':
        filepath = nom_arch + ext; 
    else:
        filepath = os.path.join(dir_arch, nom_arch + ext); 
    # Abro el archivo en modo reading
    with open(filepath, 'r') as F_open:
        # Agarro el primer elemento en F_open.readlines()
        headers = F_open.readlines()[0]; 
        # Defino L_out con headers
        L_out = headers.rstrip().split(sep=sep); 
    return L_out


def abrir_csv_listas(nom_arch, col_lista=[], path_arch='', ext='.csv', sep=';', sub_sep=',', ignore_headers=True):
    # Funcion para abrir un csv que incluye columnas con listas de elementos

    # Primero abro el .csv normalmente con abrir_csv()
    M_csv = abrir_csv(nom_arch, path_arch=path_arch, ext=ext, sep=sep, ignore_headers=ignore_headers); 
    # Recorro M_csv
    for i in range(len(M_csv)):
        curr_row = M_csv[i]; 
        # Recorro col_lista para ver a que columnas hay que hacerle split
        for j in col_lista:
            # Veo que curr_row tenga elemento j y que curr_row[j] sea distinto de string vacio antes de hacer split
            if (len(curr_row)>j) and curr_row[j]!='':
                # Si no es string vacio, uso split con sub_sep para transformarlo en lista
                M_csv[i][j] = curr_row[j].split(sub_sep); 
            elif len(curr_row)>j:
                # Si es string vacio, paso a lista vacia
                M_csv[i][j] = []; 
            else:
                # Si no tiene elemento j hago append de lista vacia
                M_csv[i].append([]); 
    return M_csv


def display_pca(M_pca, pc_var, pc, titulo='human', n_pc=5, n_ej=5):
    print()
    print('>>> PCA ' + titulo)
    print()
    print('Varianza representada por primeros ' + str(n_pc) + ' PC:')
    print(pc_var[:n_pc])
    print('Primeros ' + str(n_pc) + ' PC:')
    for i in range(n_pc):
        print(list(pc[i,]))
    print('Ejemplos de M_pca:')
    shuffle(M_pca)
    for j in M_pca[:n_ej]:
        print(j)


def dist_gene(pos_ini, pos_end, gene_id, genoma):
    # Funcion para calcular la distancia del sitio (pos_ini,pos_end) a gene.start de gene_id

    # Busco gene_id en genoma
    curr_gen = genoma.gene_by_id(gene_id); 
    # Defino posicion inicial en base a strand
    if curr_gen.strand=='+':
        gen_start = curr_gen.start; 
        gen_forward = True; 
    elif curr_gen.strand=='-':
        gen_start = curr_gen.end; 
        gen_forward = False; 
    else:
        print('WARNING: Strand no se pudo procesar en gen ' + str(curr_gen) + ', asumo forward')
        gen_start = curr_gen.start; 
        gen_forward = True; 
    # Veo que gen_start no este entre pos_ini y pos_end
    if (gen_start<pos_end) and (gen_start>pos_ini):
        # Si esta entre pos_ini y pos_end, la distancia es 0
        d_out = 0; 
    elif gen_start<pos_ini:
        # Si gen_start es menor a pos_ini, calculo distancia a pos_ini
        if gen_forward:
            # Si el gen es forward y pos_ini es mayor a gen_start, la distancia es negativa
            d_out = gen_start-pos_ini; 
        else:
            # Si el gen es reverse y pos_ini es mayor a gen_start, la distancia es positiva
            d_out = pos_ini-gen_start; 
    elif gen_start>pos_end:
        # Si gen_start es mayor a pos_end, calculo distancia a pos_end
        if gen_forward:
            # Si el gen es forward y pos_end es menor a gen_start, la distancia es positiva
            d_out = gen_start-pos_end; 
        else:
            # Si el gen es reverse y pos_end es menor a gen_start, la distancia es negativa
            d_out = pos_end-gen_start; 
    else:
        print('ERROR: Esto no deberia pasar')
    return d_out


def extraer_sitios(nom_csv, path_csv='', sep_csv=';'):
    # Abre un .csv y extrae informacion relevante a sitios de union
    # Devuelve una matriz con sitios de union e informacion relevante para procesamiento

    # Extraigo la matriz de nom_csv en path_csv
    M_csv = abrir_csv(nom_csv, path_arch=path_csv, sep=sep_csv, ignore_headers=True); 
    # Extraigo headers del archivo nom_csv, aplico .lower() para estandarizar titulos de columnas
    L_head_csv = [x.lower() for x in abrir_csv_headers(nom_csv, dir_arch=path_csv, sep=sep_csv)]; 
    # Inicializo la matriz de sitios de union
    M_sitios = []; 
    # Defino posicion de columnas que me interesan
    col_chr_n = L_head_csv.index('chr_n'); 
    col_contig = L_head_csv.index('contig'); 
    col_su_L = L_head_csv.index('pos_su_lista'); 
    col_seq_su_L = L_head_csv.index('seq_su_lista'); 
    col_su_pssm_nkx25 = L_head_csv.index('pos_pssm_nkx25'); 
    col_genes_rnaseq = L_head_csv.index('genes_rnaseq'); 
    # Recorro M_csv
    for i in range(len(M_csv)):
        curr_peak = M_csv[i]; 
        # Extraigo la informacion de las columnas que me interesan
        su_L = pasar_a_lista(curr_peak[col_su_L], ','); 
        seq_su_L = pasar_a_lista(curr_peak[col_seq_su_L], ','); 
        su_pssm_nkx25 = pasar_a_lista(curr_peak[col_su_pssm_nkx25], ','); 
        genes_rnaseq = pasar_a_lista(curr_peak[col_genes_rnaseq], ','); 
        chr_n = curr_peak[col_chr_n]; 
        contig = curr_peak[col_contig]; 
        ### Prueba
        if chr_n != ('chr'+str(contig)):
            print('WARNING: contig ' + str(contig) + ' no se corresponde a ' + str(chr_n))
        ###
        # Veo si su_L y seq_su_L tienen un solo elemento
        if (len(su_L) == 1) and (len(seq_su_L) == 1):
            # Si tienen un solo elemento, lo cargo a M_sitios con secuencia definida
            M_sitios.append([str(contig)] + su_L + seq_su_L + [genes_rnaseq]); 
        else:
            # En el resto de los casos, recorro su_L
            for curr_su_L in su_L:
                # Si hay mas de un sitio, cargo cada sitio de union por separado sin secuencia definida
                M_sitios.append([str(contig), str(curr_su_L), '', genes_rnaseq]); 
        # Recorro su_pssm_nkx25
        for curr_su_pssm in su_pssm_nkx25:
            # Cargo cada sitio de union por separado sin secuencia definida
            M_sitios.append([str(contig), str(curr_su_pssm), '', genes_rnaseq]); 
    return M_sitios


def generar_datos_grafico_pca(M_pca, M_ref, col_ref, n_pca):
    # Funcion que maneja la generacion de datos un grafico de componentes principales (acepta n_pca=2 y n_pca=3, para n_pca>3 tira warning y usa n_pca=3)
    # M_pca es la matriz con los valores de componentes principales correspondientes a M_ref
    # col_ref es la/s columna/s usada/s dentro de M_ref para usar de referencia para el grafico

    # Inicializo la matriz del grafico
    M_graf = []; 
    # Veo que M_pca y M_ref tengan el mismo largo
    if len(M_pca)!=len(M_ref): 
        print('WARNING: M_pca y M_ref tienen largos distintos. Se usa el largo minimo entre ambos.')
        len_m = min(len(M_pca), len(M_ref)); 
    else:
        len_m = len(M_pca); 
    # Recorro M_pca y M_ref
    for i in range(len_m):
        curr_pca = M_pca[i]; 
        curr_ref_row = M_ref[i]; 
        curr_ref = curr_ref_row[col_ref]; 
        # Agrego elementos acordes a M_graf
        if n_pca <= 2: 
            if n_pca < 2:
                print('WARNING: n_pca menor a 2, se usa n_pca=2.')
            M_graf.append([curr_ref, curr_pca[0], curr_pca[1]]); 
        elif n_pca >= 3:
            if n_pca > 3:
                print('WARNING: n_pca mayor a 3, se usa n_pca=3.')
            M_graf.append([curr_ref, curr_pca[0], curr_pca[1], curr_pca[2]]); 
    return M_graf


def generar_grafico_2d(M_data, col_x=0, col_y=1, col_label=2, L_col_color=[3,4,5], plot_title='PCA', nom_out='', path_out='.\\', ext_out='.png', show_plot=True, L_axis=[]):
    # Funcion para generar un grafico 2d con los datos de una matriz
    # col_x y col_y determinan las columnas de los datos correspondientes a esos ejes en M_data
    # L_col_color determina las columnas que definen el color de cada punto, contiene 3 numeros entre 0 y 1
    # col_label determina la columna que define la clasificacion de cada punto
    
    # Uso matplotlib para hacer el grafico
    fig, ax = plt.subplots(); 
    fig.set_size_inches(16, 16); 
    # Defino labels
    if len(L_axis)>=2:
        plt.xlabel(L_axis[0]); 
        plt.ylabel(L_axis[1]); 
    else:
        plt.xlabel('x'); 
        plt.ylabel('y'); 
    # Recorro M_data
    for i in range(len(M_data)):
        curr_row = M_data[i]; 
        # Defino si uso color
        if L_col_color:
            # Defino la tupla de color
            curr_color = (float(curr_row[L_col_color[0]]),float(curr_row[L_col_color[1]]),float(curr_row[L_col_color[2]])); 
            # Agrego los datos con ax.scatter()
            ax.scatter(float(curr_row[col_x]), float(curr_row[col_y]), c=curr_color, label=curr_row[col_label]); 
        else:
            # Agrego los datos con ax.scatter()
            ax.scatter(float(curr_row[col_x]), float(curr_row[col_y]), c='black', label=curr_row[col_label]); 
    # Agrego ax.legend() para ver labels
    #ax.legend(); 
    # Defino titulo
    plt.title(plot_title); 
    # Si nom_out es distinto de string vacio, guardo el grafico
    if nom_out!='':
        # Defino full_out con path_out, nom_out y ext_out
        full_out = os.path.join(path_out, nom_out+ext_out); 
        # Guardo la figura con plt.savefig()
        plt.savefig(full_out, bbox_inches='tight'); 
    # Si show_plot es True, muestro el plot
    if show_plot:
        plt.show(); 
    else:
        plt.close(); 
    return M_data


def generar_grafico_3d(M_data, col_x=0, col_y=1, col_z=2, col_label=3, L_col_color=[4,5,6], plot_title='PCA', nom_out='', path_out='.\\', ext_out='.png', show_plot=True, L_axis=[]):
    # Funcion para generar un grafico 3d con los datos de una matriz
    # col_x, col_y y col_z determinan las columnas de los datos correspondientes a esos ejes en M_data
    # L_col_color determina las columnas que definen el color de cada punto, contiene 3 numeros entre 0 y 1
    # col_label determina la columna que define la clasificacion de cada punto

    # Uso matplotlib para hacer el grafico
    fig = plt.figure(figsize=(16, 16)); 
    ax = fig.add_subplot(projection='3d'); 
    # Defino labels
    if len(L_axis)>=3:
        ax.set_xlabel(L_axis[0]); 
        ax.set_ylabel(L_axis[1]); 
        ax.set_zlabel(L_axis[2]); 
    else:
        ax.set_xlabel('x'); 
        ax.set_ylabel('y'); 
        ax.set_zlabel('z'); 
    # Recorro M_data
    for i in range(len(M_data)):
        curr_row = M_data[i]; 
        # Defino si uso color
        if L_col_color:
            # Defino la tupla de color
            curr_color = (float(curr_row[L_col_color[0]]),float(curr_row[L_col_color[1]]),float(curr_row[L_col_color[2]])); 
            # Agrego los datos con ax.scatter()
            ax.scatter(float(curr_row[col_x]), float(curr_row[col_y]), float(curr_row[col_z]), c=curr_color, label=str(curr_row[col_label])); 
        else:
            # Agrego los datos con ax.scatter()
            ax.scatter(float(curr_row[col_x]), float(curr_row[col_y]), float(curr_row[col_z]), c='black', label=str(curr_row[col_label])); 
    # Agrego ax.legend() para ver labels
    #ax.legend(); 
    # Defino titulo
    plt.title(plot_title); 
    # Si nom_out es distinto de string vacio, guardo el grafico
    if nom_out!='':
        # Defino full_out con path_out, nom_out y ext_out
        full_out = os.path.join(path_out, nom_out+ext_out); 
        # Guardo la figura con plt.savefig()
        plt.savefig(full_out, bbox_inches='tight'); 
    # Si show_plot es True, muestro el plot
    if show_plot:
        plt.show(); 
    else:
        plt.close(); 
    return M_data


def generar_grafico_pca2(M_data, nom_out, col_data=['seq', 'PC1', 'PC2'], col_num=[0,1,2], path_out='', ext_out='.png'):
    # Funcion para generar scatterplot 2D con M_data
    # col_data es una lista con las 3 columnas de M_data que correspondan en col_num
    # col_num determina col_data segun el formato [clasificacion, x, y]

    # Creo dataframe con M_data
    pca_df = pd.DataFrame(data=M_data, columns=[col_data[col_num[0]], col_data[col_num[1]], col_data[col_num[2]]]); 
    # Uso seaborn para hacer el grafico
    sns.lmplot(x=col_data[col_num[1]], y=col_data[col_num[2]], data=pca_df, hue=col_data[col_num[0]], fit_reg=False, legend=True); 
    # Defino el titulo
    plt.title('PCA ' + col_data[col_num[0]]); 
    # Defino full_out con path_out, nom_out y ext_out
    full_out = os.path.join(path_out, nom_out+ext_out); 
    plt.savefig(full_out, bbox_inches='tight'); 
    plt.show(); 
    return M_data


def generar_grafico_pca3(M_data, nom_out, col_data=['seq', 'PC1', 'PC2', 'PC3'], col_num=[0,1,2,3], path_out='', ext_out='.png'):
    # Funcion para generar scatterplot 3D con M_data
    # col_data es una lista con las 4 columnas de M_data que correspondan en col_num
    # col_num determina col_data segun el formato [clasificacion, x, y, z]

    # Inicializo lista de colores para el grafico
    L_color = []; 
    # Creo matriz para el grafico
    M_graf = []; 
    # Recorro M_data
    for i in range(len(M_data)):
        curr_row = M_data[i]; 
    ### FALTA
    # Pasar lo de abajo a formato de arriba
    ###
    # Creo dataframe con M_data
    pca_df = pd.DataFrame(data=M_data, columns=[col_data[col_num[0]], col_data[col_num[1]], col_data[col_num[2]], col_data[col_num[3]]]); 
    # Uso matplotlib para hacer el grafico
    fig = plt.figure(figsize=(12, 12)); 
    ax = fig.add_subplot(projection='3d'); 
    # Defino x, y, z, y color
    x = pca_df.iloc[:,1]; 
    y = pca_df.iloc[:,2]; 
    z = pca_df.iloc[:,3]; 
    # Defino la lista de colores por posicion
    color = list_to_color(x,y,z); 
    # Paso x, y, z y seq a lista
    L_seq = list(pca_df.iloc[:,0]); 
    L_x = list(x); 
    L_y = list(y); 
    L_z = list(z); 
    # Inicializo una matriz para reordenar x, y, z, seq y color
    M_graf = []; 
    # Recorro x, y, z, seq y color para pasarlos a matriz
    for m in range(len(L_seq)):
        # Agarro solo los elementos unicos de M_graf
        if not ([L_x[m],L_y[m],L_z[m], color[m], L_seq[m]] in M_graf):
            M_graf.append([L_x[m],L_y[m],L_z[m], color[m], L_seq[m]]); 
    # Recorro M_graf para hacer el plot
    for m in range(len(M_graf)):
        ax.scatter(M_graf[m][0],M_graf[m][1],M_graf[m][2], c=M_graf[m][3], label=M_graf[m][4]); 
    #ax.legend(M_leg, L_seq); 
    #ax.legend(loc='upper left'); 
    ax.legend(); 
    plt.title('PCA secuencias'); 
    # Defino full_out con path_out, nom_out y ext_out
    full_out = os.path.join(path_out, nom_out+ext_out); 
    plt.savefig(full_out, bbox_inches='tight'); 
    plt.show(); 
    return M_data, M_graf


def generar_M_log(M_pca_in, M_pca_std, M_pca_out, M_pca_graph, pc_var, pc):
    # Funcion que registra todas las columnas importantes de las distintas matrices y listas generadas y las guarda en una (o dos?) variables

    # Inicializo la matriz que se devuelve
    M_log = []; 
    # Agrego la info de varianza de componentes principales
    M_log.append(['Varianza de componentes principales']); 
    M_log.append(list(pc_var)); 
    # Agrego los componentes principales
    M_log.append(['Componentes principales']); 
    for curr_pc in pc:
        M_log.append(list(curr_pc)); 
    # Verifico que las distintas matrices tengan todas el mismo largo
    len_pca_in = len(M_pca_in); 
    len_pca_std = len(M_pca_std); 
    len_pca_out = len(M_pca_out); 
    len_pca_graph = len(M_pca_graph); 
    if len_pca_in==len_pca_std==len_pca_out==len_pca_graph:
        len_usada = int(len_pca_in); 
    else:
        print('WARNING: Largos de matrices difieren entre si. Se usa el valor minimo.')
        print('M_pca_in: ' + str(len_pca_in))
        print('M_pca_std: ' + str(len_pca_std))
        print('M_pca_out: ' + str(len_pca_out))
        print('M_pca_graph: ' + str(len_pca_graph))
        len_usada = min(len_pca_in, len_pca_std, len_pca_out, len_pca_graph); 
    # Defino titulos de la tabla
    head_pca_in = ['chr_n', 'pos_ini', 'pos_end', 'seq_full', 'gene_id', 'log_fc', 'seq_cut', 'dist', 'updown', 'dist_transform', 'seq_onehot']; 
    head_pca_std = ['std_updown', 'std_dist', 'std_seq']; 
    head_pca_out = ['PC1', 'PC2', 'PC3']; 
    head_pca_graph = ['ref', 'x', 'y', 'z']; 
    head_pca_in = std_len_head(head_pca_in, len(M_pca_in[0]), ext=''); 
    head_pca_std = std_len_head(head_pca_std, len(M_pca_std[0]), ext=''); 
    head_pca_out = std_len_head(head_pca_out, len(M_pca_out[0]), ext=''); 
    head_pca_graph = std_len_head(head_pca_graph, len(M_pca_graph[0]), ext=''); 
    # Agrego el titulo a M_log
    M_log.append(head_pca_in + head_pca_std + head_pca_out + head_pca_graph); 
    # Recorro las 4 matrices
    for i in range(len_usada): 
        curr_pca_in = M_pca_in[i]; 
        curr_pca_std = list(M_pca_std[i,:]); 
        curr_pca_out = list(M_pca_out[i,:]); 
        curr_pca_graph = M_pca_graph[i]; 
        # Agrego cada elemento de cada una de las filas a M_log
        M_log.append(curr_pca_in + curr_pca_std + curr_pca_out + curr_pca_graph); 
    return M_log


def generar_PCA(M_in, n_pc=2):
    # Funcion para generar datos relevantes a analisis de componentes principales

    # Transformo M_in en array de numpy
    arr_in = StandardScaler().fit_transform(array(M_in)); 
    # Uso pca() para generar datos importantes
    _col_mean, coord, pc, ev = pca(arr_in); 
    # Selecciono las n_pc primeras columnas de coord
    M_pc = coord[:,:n_pc]; 
    # Uso eigenvalues para calcular varianza explicada por cada pc
    pc_var = [x/sum(ev) for x in ev]; 
    # Devuelvo la matriz seleccionada (M_pc), la varianza explicada (pc_var) y los componentes principales (pc)
    return M_pc, pc_var, pc, arr_in


def generar_pca_secuencias(nom_arch, nom_out, path_out='', path_arch='', sep_in=';', ext_in='.csv', n_pca=3, ext_out='.png'):
    # Funcion para hacer analisis de componentes principales solo sobre las secuencias en nom_arch
    # Codifica cada base de las secuencias como 4 booleanos (codificacion one-hot por nucleotido)

    # Abro el csv nom_arch en path_arch
    M_csv = abrir_csv(nom_arch, path_arch=path_arch, ext=ext_in, sep=sep_in, ignore_headers=True); 
    # Inicializo la matriz que devuelve la funcion
    M_out = []; 
    # Inicializo la matriz que se va a usar para PCA
    M_pca_in, M_out = procesar_seq_para_pca(M_csv, min_row_len=4, seq_col=3); 
    # Uso generar_pca() con M_pca_in
    M_pca_out, pc_var, pc, _arr_in = generar_PCA(M_pca_in, n_pc=n_pca); 
    # Inicializo matriz con info para grafico
    M_pca_graph = []; 
    # Recorro M_pca_out
    for i in range(len(M_pca_out)): 
        curr_pca_out = M_pca_out[i]; 
        curr_out = M_out[i]; 
        # Agrego los valores relevantes a M_pca_graph
        M_pca_graph.append([curr_out[6]] + list(curr_pca_out)); 
    # Si n_pca es 2
    if n_pca == 2:
        # Uso generar_grafico_pca2(), que usa seaborn para un grafico 2d
        M_pca_graph = generar_grafico_pca2(M_pca_graph, nom_out, col_data=['seq', 'PC1', 'PC2'], col_num=[0,1,2], path_out=path_out, ext_out=ext_out); 
    # Si n_pca es 3 o mas
    elif n_pca >= 3:
        # Si n_pca es mayor a 3, tiro advertencia
        if n_pca > 3:
            print('WARNING: El numero de componentes principales seleccionado es mayor a 3. Se grafican los primeros 3.')
        ### CUIDADO
        # Esta funcion no esta terminada
        ###
        # Uso generar_grafico_pca3(), que usa matplotlib para un grafico 3d
        M_pca_graph, M_graf_out = generar_grafico_pca3(M_pca_graph, nom_out, col_data=['seq', 'PC1', 'PC2', 'PC3'], col_num=[0,1,2,3], 
                                                       path_out=path_out, ext_out=ext_out); 
    # Hardcodeo titulos de las matrices
    L_titulos = [['chr_n', 'pos_ini', 'pos_end', 'seq_raw', 'gene', 'log_fc', 'seq_cut', 'seq_onehot'], ['PCA_in'], ['PC1', 'PC2', 'PC3'], ['ref', 'x', 'y', 'z']]; 
    # Junto las matrices con datos relevantes en M_log
    M_log = juntar_M_log([M_out, M_pca_in, M_pca_out, M_pca_graph], L_titulos); 
    # Guardo M_log
    M_log = guardar_csv(M_log, nom_out, path_out=path_out); 
    return M_out, pc_var, pc


def generar_pca_seq_one_hot(nom_arch, nom_out, genoma, path_out='.\\', path_arch='', sep_in=';', ext_in='.csv', n_pca=2, dist_t='i', onehot_full_seq=True):
    # Funcion para hacer analisis de componentes principales solo sobre las secuencias en nom_arch
    # Guarda el output en un archivo .csv para ser 
    # Codifica cada secuencia posible como un booleano de una lista (codificacion one-hot por secuencia)
    # Puede usar distancia como variable en z para diferenciar puntos apelmazados
    # dist_t (dist transform) determina el tipo de transformacion que se le hace a la distancia
        # dist_t='i' se usa 1/dist
        # dist_t='r' se usa simplificar_dist() para clasificar en rangos
        # dist_t='l' se usa log_para_grafico()
        # transform='li' o transform='il' se usa 1/log_para_grafico()
        # dist_t='n' o dist_t='' se usa la distancia normalmente
        # Cualquier otro valor tira warning y se usa la distancia normalmente

    # Abro el csv nom_arch en path_arch
    M_csv = abrir_csv(nom_arch, path_arch=path_arch, ext=ext_in, sep=sep_in, ignore_headers=True); 
    # Inicializo la matriz que registra datos de M_pca_in y M_csv
    M_reg = []; 
    # Inicializo la matriz que se usa para PCA
    M_pca_in = []; 
    # Inicializo una matriz de registro de datos
    M_log = []; 
    # Inicializo diccionario de secuencias para hacer one-hot con secuencia completa
    onehot_dict = {}; 
    # Recorro M_csv
    for i in range(len(M_csv)):
        curr_row = M_csv[i]; 
        # Solo agarro curr_row si tiene info de gen cercano
        if len(curr_row) > 4:
            if onehot_full_seq:
                # Defino secuencia con procesar_seq() y confirmo si esta en onehot_dict.keys() o si hay que agregarla
                _, seq_recortada = procesar_seq(curr_row[3]); 
                if not(str(seq_recortada) in onehot_dict.keys()):
                    onehot_dict[str(seq_recortada)] = [0]*len(onehot_dict.keys()) + [1]; 
                seq_usada = onehot_dict[str(seq_recortada)]; 
            else:
                seq_usada, seq_recortada = procesar_seq(curr_row[3]); 
            # Defino distancia usada con dist_gene() y procesar_dist()
            curr_dist = dist_gene(int(curr_row[1]), int(curr_row[2]), curr_row[4], genoma); 
            curr_dist_usada = procesar_dist(curr_dist, transform=dist_t); 
            # Defino si el gen esta upregulado o downregulado
            curr_fc = float(curr_row[5]); 
            if curr_fc > 0:
                updownreg = 1; 
            elif curr_fc < 0:
                updownreg = -1; 
            # Agrego las columnas relevantes a M_pca_in
            M_pca_in.append([updownreg, curr_dist_usada] + seq_usada); 
            # Agrego todas las columnas posibles a M_reg
            M_reg.append(curr_row + [seq_recortada, curr_dist, updownreg, curr_dist_usada] + seq_usada); 
    # Estandarizo el largo de M_pca_in y M_reg (algunas seq_usada terminan antes)
    if onehot_full_seq:
        M_pca_in = std_len_M_pca(M_pca_in); 
        M_reg = std_len_M_pca(M_reg); 
    # Uso generar_pca() con M_pca_in
    M_pca_out, pc_var, pc, M_pca_std = generar_PCA(M_pca_in, n_pc=n_pca); 
    # Defino la posicion del valor de referencia en M_reg
    col_ref = 1; 
    # Inicializo matriz con info para grafico
    M_pca_graph = generar_datos_grafico_pca(M_pca_out, M_pca_in, col_ref, n_pca); 
    # Genero M_log con las matrices y datos que fui armando
    M_log = generar_M_log(M_reg, M_pca_std, M_pca_out, M_pca_graph, pc_var, pc); 
    # Guardo M_log en el archivo nom_out dentro de path_out
    M_log = guardar_csv(M_log, nom_out, path_out=path_out, sep_out=sep_in); 
    return M_log


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


def invertir_seq(seq):
    # Funcion para generar el reverso complementario de seq

    # Diccionario de nucleotidos 
    dict_rev = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}; 
    # Inicializo la secuencia que se devuelve
    seq_out = ''; 
    # Recorro seq
    for i in range(len(seq)):
        # Recorro seq desde el final hasta el principio
        curr_n = seq[-i-1]; 
        # Agrego el reverso segun dict_rev a seq_out
        seq_out = seq_out + dict_rev[curr_n]; 
    return seq_out


def juntar_M_log(L_matrices, L_titulos):
    # Funcion para juntar una lista de matrices como una sola matriz
    # Asume que tienen la misma cantidad de filas

    # Inicializo la matriz que se devuelve
    M_out = []; 
    # Inicializo listas de largos de filas para cada matriz
    L_row_len = []; 
    # Inicializo len_m
    len_m = len(L_matrices[0]); 
    # Recorro L_matrices para ver que todos los componentes tengan el mismo largo
    for i in range(len(L_matrices)):
        curr_M = L_matrices[i]; 
        # Veo que curr_M tenga largo len_m
        if len(curr_M)!=len_m:
            print('WARNING: Matriz numero ' + str(i+1) + ' tiene largo distinto que len_m.')
            len_m = min(len_m, len(curr_M)); 
        # Agrego el largo de curr_M[0] a L_row_len
        L_row_len.append(int(len(curr_M[0]))); 
    # Veo que L_titulos sea igual de largo (o mas) que L_row_len antes de agregar titulos
    if len(L_titulos) >= len(L_row_len):
        # Si L_titulos es mas grande que L_row_len, tiro warning
        if len(L_titulos) > len(L_row_len):
            print('WARNING: Largo de L_titulos es mayor a la cantidad de matrices dadas. Se descartan los titulos fuera del rango.')
        # Inicializo la lista de titulos
        L_head = []; 
        # Recorro L_row_len
        for l in range(len(L_row_len)):
            # Veo si L_titulos[l] es del largo en L_row_len[l]
            if len(L_titulos[l])==L_row_len[l]:
                # Si es del largo exacto, simplemente lo agrego
                L_head = L_head + L_titulos[l]; 
            elif len(L_titulos[l])>L_row_len[l]:
                # Si L_titulos[l] es mas largo, tiro warning y uso los primeros L_row_len[l] elementos
                L_head = L_head + L_titulos[l][:L_row_len[l]]; 
            elif len(L_titulos[l])<L_row_len[l]:
                # Si L_row_len[l] es mas largo, agrego numeros al final
                L_head = L_head + L_titulos[l] + [str(l)]*(L_row_len[l]-len(L_titulos[l])); 
        # Agrego L_head a M_out
        M_out.append(L_head); 
    # Con el largo de filas a recorrer definido, hago el recorrido
    for j in range(len_m):
        # Inicializo la lista correspondiente a la fila
        L_row = []; 
        # Recorro L_matrices para agarrar una fila de cada matriz
        for k in range(len(L_matrices)):
            curr_M = L_matrices[k]; 
            # Agrego la fila j de curr_M a L_row
            L_row = L_row + list(curr_M[j]); 
        # Agrego L_row a M_out
        M_out.append(L_row[:]); 
    return M_out


def list_to_color(x,y,z,mult=0.75):
    # Funcion para pasar tres listas de valores (del mismo largo) en tres listas valores entre 0 y 1, para usar como color
    # x, y, y z son listas de numeros que pueden incluir valores negativos

    # Uso list_to_color_transform() para devolver las listas con formato para color
    c_list = (list_to_color_transform(x,multiplier=mult),list_to_color_transform(y,multiplier=mult),list_to_color_transform(z,multiplier=mult)); 
    # Inicializo la lista que se devuelve
    c_out = []; 
    # Recorro c_list
    for i in range(len(c_list[0])):
        c_out.append((c_list[0][i], c_list[1][i], c_list[2][i])); 
    return c_out


def list_to_color_transform(L_x, multiplier=0.75):
    # Funcion para transformar una lista de valores positivos y negativos en una lista de valores entre 0 y 1
    # Si hay valores negativos, sumo el valor absoluto del minimo para que sean todos positivos
    if min(L_x)<0:
        L_abs = [i+abs(min(L_x)) for i in L_x]; 
    else:
        L_abs = L_x; 
    # Devuelvo la lista con todos los valores divididos por el maximo, para tener valores entre 0 y 1
    return [j*multiplier/max(L_abs) for j in L_abs]


def log_para_grafico(n, base=10):
    # Funcion que aplica escala logaritmica a un valor para graficar
    # Calcula el logaritmo del valor absoluto y aplica el signo de n
    # Devuelve 0 para valores entre 1 y -1

    # Veo si el valor esta entre 1 y -1
    if (n<1) and (n>-1):
        ret = 0; 
    else:
        ret = math.copysign(math.log(abs(n),base), n); 
    return ret


def pasar_a_lista(s, sep):
    # Funcion para pasar string a lista segun sep usando .split()
    # Solo se encarga del caso string vacio, que devuelve lista vacia en vez de una lista con string vacio
    if len(s):
        L_out = s.split(sep); 
    else:
        L_out = []; 
    return L_out


def pipeline_agregar_info_rnaseq(nom_su, nom_rnaseq, nom_out, path_su='', path_rnaseq='', path_out='.\\', col_listas=[], col_genes=4):
    # Funcion que abre el archivo generado por pipeline_generar_lista_sitios() y selecciona genes que esten up o downregulados en ensayos RNA-seq
    # Los datos de up o downregulacion de RNA-seq salen de script 9-ClasificacionPeaks3.py
    # Agrega info de up o downregulacion para cada gen registrado

    ### RNA-seq
    # Abro el archivo con info de genes up o downregulados
    M_rnaseq = abrir_csv(nom_rnaseq, path_arch=path_rnaseq, ignore_headers=True); 
    # Creo una lista de los Ensembl IDs en M_rnaseq
    L_ensembl = []; 
    # Recorro M_rnaseq
    for r in range(len(M_rnaseq)):
        curr_rnaseq = M_rnaseq[r]; 
        # Agrego la columna de Ensembl ID a L_ensembl
        L_ensembl.append(curr_rnaseq[1]); 
    ###
    ### Recorrer matriz de sitios de union
    # Abro el archivo con info de sitios de union
    M_su = abrir_csv_listas(nom_su, col_lista=col_listas, path_arch=path_su); 
    # Inicializo la matriz que se devuelve
    M_out = []; 
    # Recorro M_su
    for i in range(len(M_su)):
        curr_su = M_su[i]; 
        # Reviso los genes de curr_su
        curr_L_genes = curr_su[col_genes]; 
        # Booleano para ver si no se encuentra ningun gen
        sin_gen = True; 
        # Recorro la lista de genes
        for j in range(len(curr_L_genes)):
            curr_gene = curr_L_genes[j]; 
            # Reviso que curr_gene este en L_ensembl
            if curr_gene in L_ensembl:
                ### Chequeo
                if M_rnaseq[L_ensembl.index(curr_gene)][1]!=curr_gene:
                    print('WARNING: ' + str(curr_gene) + ' no encontrado en ' + str(M_rnaseq[L_ensembl.index(curr_gene)]))
                ###
                sin_gen = False; 
                # Agrego la informacion en M_rnaseq en el ID correspondiente a curr_gene en L_ensembl
                M_out.append(curr_su[:-1] + M_rnaseq[L_ensembl.index(curr_gene)][1:]); 
        # Si no se encuentra ningun gen, se hace append de curr_su[:-1]
        if sin_gen:
            M_out.append(curr_su[:-1]); 
    ###
    # Guardo M_out
    M_out = guardar_csv(M_out, nom_out, path_out=path_out, L_head=['chr_n', 'pos_ini', 'pos_end', 'seq', 'gene', 'log_fc']); 
    return M_out


def pipeline_generar_lista_sitios(nom_arch, genoma, nom_genoma, nom_out, path_csv='', path_fasta='', path_out='.\\', sep_csv=';', sep_out=';'):
    # Funcion para generar la lista de sitios
    # Falta procesar genes RNA-seq

    # Extraigo sitios en nom_arch con extraer_sitios()
    M_sitios = extraer_sitios(nom_arch, path_csv=path_csv, sep_csv=sep_csv); 
    # Proceso sitios extraidos con procesar_sitios()
    M_sitios = procesar_sitios(M_sitios, genoma, nom_genoma, path_fasta=path_fasta); 
    # Guardo M_sitios con guardar_csv()
    M_sitios = guardar_csv(M_sitios, nom_out, path_out=path_out, sep_out=sep_out, L_head=['chr_n','pos_ini','pos_end','seq','near_genes']); 
    return M_sitios


def pipeline_graficos_finales(nom_full, nom_mod='', path_full='', path_out='.\\'):
    # Funcion para hacer graficos 3d y 2d con los datos en las distintas tablas que fui armando

    # Columnas del archivo
    col_log_fc = 5; 
    col_seq = 6; 
    col_dist = 7; 
    col_updown = 8; 
    col_dist_inv = 9; 
    col_dist_log = 10; 
    col_dist_inv_log = 11; 
    col_pc_dist_log1 = 12; 
    col_pc_dist_log2 = 13; 
    col_pc_dist_log3 = 14; 
    col_pc_seq1 = 15; 
    col_pc_seq2 = 16; 
    col_pc_seq3 = 17; 
    col_red = 21; 
    col_green = 22; 
    col_blue = 23; 
    col_pc_dist_inv1 = 24; 
    col_pc_dist_inv2 = 25; 
    col_pc_dist_inv3 = 26; 
    col_pc_dist_norm1 = 27; 
    col_pc_dist_norm2 = 28; 
    col_pc_dist_norm3 = 29; 
    # Abro el csv con los datos full que necesito para esto
    M_full = abrir_csv(nom_full, path_arch=path_full, ignore_headers=True); 

    # Grafico de distancia vs updown (log_fc en eje z) con labels por secuencia
    #generar_grafico_2d(M_full, col_x=col_dist, col_y=col_log_fc, col_label=col_seq, L_col_color=[], plot_title='dist vs updown', 
    #                   nom_out='dist_vs_updown_'+nom_mod+'_unique_color', path_out=path_out, ext_out='.png', show_plot=False, L_axis=['dist_norm', 'log_fc']); 
    # Grafico de distancia vs updown (log_fc en eje z) con labels por secuencia y colores por PCA de secuencia
    #generar_grafico_2d(M_full, col_x=col_dist, col_y=col_log_fc, col_label=col_seq, L_col_color=[col_red,col_green,col_blue], plot_title='dist vs updown + color', 
    #                   nom_out='dist_vs_updown_'+nom_mod+'_pca_color', path_out=path_out, ext_out='.png', show_plot=False, L_axis=['dist_norm', 'log_fc']); 

    # Grafico de PCA_full para 2 PC con distancia en escala logaritmica
    generar_grafico_2d(M_full, col_x=col_pc_dist_log1, col_y=col_pc_dist_log2, col_label=col_seq, L_col_color=[col_red,col_green,col_blue], 
                       plot_title='pc dist log 1,2', nom_out='pc_dist_log_1_2_'+nom_mod+'_seq_color', path_out=path_out, ext_out='.png', 
                       show_plot=False, L_axis=['PC1 dist log', 'PC2 dist log']); 
    # Grafico de PCA_full para 2 PC con distancia en escala logaritmica vs distancia absoluta en eje z
    generar_grafico_3d(M_full, col_x=col_pc_dist_log1, col_y=col_pc_dist_log2, col_z=col_dist, col_label=col_seq, L_col_color=[col_red,col_green,col_blue], 
                       plot_title='pc dist log 1,2 vs dist', nom_out='pc_dist_log_1_2_vs_dist_'+nom_mod+'_seq_color', path_out=path_out, ext_out='.png', 
                       show_plot=False, L_axis=['PC1 dist log', 'PC2 dist log', 'dist_norm']); 

    # Grafico de PCA_full para 2 PC con distancia invertida
    generar_grafico_2d(M_full, col_x=col_pc_dist_inv1, col_y=col_pc_dist_inv2, col_label=col_seq, L_col_color=[col_red,col_green,col_blue], 
                       plot_title='pc dist inv 1,2', nom_out='pc_dist_inv_1_2_'+nom_mod+'_seq_color', path_out=path_out, ext_out='.png', 
                       show_plot=False, L_axis=['PC1 dist inv', 'PC2 dist inv']); 
    # Grafico de PCA_full para 2 PC con distancia invertida vs distancia absoluta en eje z
    generar_grafico_3d(M_full, col_x=col_pc_dist_inv1, col_y=col_pc_dist_inv2, col_z=col_dist, col_label=col_seq, L_col_color=[col_red,col_green,col_blue], 
                       plot_title='pc dist inv 1,2 vs dist', nom_out='pc_dist_inv_1_2_vs_dist_'+nom_mod+'_seq_color', path_out=path_out, ext_out='.png', 
                       show_plot=False, L_axis=['PC1 dist inv', 'PC2 dist inv', 'dist_norm']); 

    # Grafico de PCA_full para 2 PC con distancia absoluta
    generar_grafico_2d(M_full, col_x=col_pc_dist_norm1, col_y=col_pc_dist_norm2, col_label=col_seq, L_col_color=[col_red,col_green,col_blue], 
                       plot_title='pc dist norm 1,2', nom_out='pc_dist_norm_1_2_'+nom_mod+'_seq_color', path_out=path_out, ext_out='.png', 
                       show_plot=False, L_axis=['PC1 dist norm', 'PC2 dist norm']); 
    # Grafico de PCA_full para 2 PC con distancia absoluta vs distancia absoluta en eje z
    generar_grafico_3d(M_full, col_x=col_pc_dist_norm1, col_y=col_pc_dist_norm2, col_z=col_dist, col_label=col_seq, L_col_color=[col_red,col_green,col_blue], 
                       plot_title='pc dist norm 1,2 vs dist', nom_out='pc_dist_norm_1_2_vs_dist_'+nom_mod+'_seq_color', path_out=path_out, ext_out='.png', 
                       show_plot=True, L_axis=['PC1 dist norm', 'PC2 dist norm', 'dist_norm']); 
    
    '''
	X Hacer grafico de distancia vs updown con colores para secuencias
		X Primero con colores por secuencia unica sin importar nada
		X Despues con colores definidos por PCA de secuencias
	~ Usar distancia como eje z en grafico PCA que ya hice
		X Probar distancia normal e invertida
		X Ver en que tabla quedaron esos datos o que script los genera
	- Hacer grafico PCA con peso de variantes relativizado 
		* Para one-hot por nucleotido seria multiplicar variables dist y updown por sqrt(7) despues de normalizacion
		* Para one-hot por secuencia completa, no cambia nada (probar correr el PCA sin alterar los booleanos?)
    '''
    return M_full


def pipeline_pca(nom_arch, genoma, nom_out, path_out='.\\', path_arch='', sep_in=';', ext_in='.csv', n_pca=2, dist_t='i'):
    # Funcion que corre funciones necesarias para hacer analisis de componentes principales sobre datos en nom_arch

    # Abro el csv nom_arch en path_arch
    M_csv = abrir_csv(nom_arch, path_arch=path_arch, ext=ext_in, sep=sep_in, ignore_headers=True); 
    # Inicializo la matriz que se va a usar para PCA
    M_pca_in = []; 
    # Inicializo la matriz que devuelve la funcion
    M_out = []; 
    # Recorro M_csv
    for i in range(len(M_csv)):
        curr_row = M_csv[i]; 
        # Solo agarro curr_row si tiene info de gen cercano
        if len(curr_row) > 4:
            # Defino variables para PCA
            curr_fc = float(curr_row[5]); 
            curr_seq, seq_usada = procesar_seq(curr_row[3]); 
            curr_dist = dist_gene(int(curr_row[1]), int(curr_row[2]), curr_row[4], genoma); 
            # Defino como uso curr_dist
            curr_dist_usada = procesar_dist(curr_dist, transform=dist_t); 
            # Defino si el gen esta upregulado o downregulado
            if curr_fc > 0:
                updownreg = 1; 
            elif curr_fc < 0:
                updownreg = -1; 
            # Defino si agrego curr_dist_usada o no
            if str(curr_dist_usada) != '':
                # Agrego las variables PCA a M_pca_in
                M_pca_in.append([updownreg, curr_dist_usada] + curr_seq); 
                # Agrego las variables de PCA y de curr_row a M_out
                M_out.append(curr_row + [seq_usada, curr_dist, updownreg, curr_dist_usada] + curr_seq); 
            else:
                # Agrego las variables PCA a M_pca_in
                M_pca_in.append([updownreg] + curr_seq); 
                # Agrego las variables de PCA y de curr_row a M_out
                M_out.append(curr_row + [seq_usada, curr_dist, updownreg] + curr_seq); 
            
    # Uso generar_pca() con M_pca_in
    M_pca_out, pc_var, pc, _arr_in = generar_PCA(M_pca_in, n_pc=n_pca); 
    # Inicializo matriz con info para grafico
    M_pca_graph = []; 
    # Recorro M_pca_out y M_out
    for i in range(len(M_pca_out)): 
        curr_pca_out = M_pca_out[i]; 
        curr_out = M_out[i]; 
        # Agrego los valores relevantes a M_pca_graph
        M_pca_graph.append([curr_out[6]] + list(curr_pca_out)); 
    # Creo dataframe con M_pca_graph
    pca_df = pd.DataFrame(data=M_pca_graph, columns=['seq', 'PC1', 'PC2']); 
    sns.lmplot(x='PC1', y='PC2', data=pca_df, hue='seq', fit_reg=False, legend=True); 
    plt.title('PCA'); 
    plt.show(); 
    # Defino listas y matrices para guardar
    L_head_save = ['chr_n', 'pos_ini', 'pos_end', 'seq', 'gene', 'log_fc', 'seq_usada', 'dist', 'regulacion', 'dist_usada', 'seq_bin']; 
    M_save = [['PC_variance'], pc_var, pc, L_head_save] + M_out; 
    M_save = guardar_csv(M_save, nom_out, path_out=path_out); 
    return M_out, pc_var, pc


def procesar_dist(dist, transform='n'):
    # Funcion que determina la distancia que se usa para PCA de acuerdo a transform
    # transform='i' se usa 1/dist
    # transform='r' se usa simplificar_dist() para clasificar en rangos
    # transform='l' se usa log_para_grafico()
    # transform='li' o transform='il' se usa 1/log_para_grafico()
    # transform='n' o transform='' se usa la distancia normalmente
    # transform='na' o transform='0' no se usa la distancia
    # Cualquier otro valor tira warning y se usa la distancia normalmente

    # Pruebo los distintos transform determinados
    if transform.lower()=='n' or transform.lower()=='':
        # Si transform es 'n' o string vacio, uso dist normalmente
        ret = dist; 
    elif transform.lower()=='i':
        # Si transform es 'i', uso la inversa de dist
        ret = 1.0/dist; 
    elif transform.lower()=='r':
        # Si transform es 'r', uso la funcion simplificar_dist() para pasarlo a rangos
        ret = simplificar_dist(dist); 
    elif transform.lower()=='l':
        # Si transform es 'l', uso la funcion log_para_grafico() para pasar a escala logaritmica
        ret = log_para_grafico(dist, base=10); 
    elif transform.lower()=='li' or transform.lower()=='il':
        # Si transform es 'li' o 'il', uso 1/log_para_grafico()
        ret = 1/log_para_grafico(dist, base=10); 
    elif transform.lower()=='na' or transform.lower()=='0':
        # Si transform es 'na' o '0', devuelvo string vacio
        ret = ''; 
    else:
        # Si no es ninguno de los valores anteriores, tiro warning y uso dist normalmente
        print('WARNING: transform=' + transform + ' no se corresponde con las transformaciones aceptadas. Se usa dist normalmente.')
        ret = dist; 
    return ret


def procesar_seq_para_pca(M_csv, min_row_len=4, seq_col=3):
    # Funcion que procesa una matriz de sitios de union con genes cercanos y devuelve una matriz para dar a generar_pca()

    # Inicializo la matriz que se devuelve
    M_pca_in = []; 
    # Inicializo una matriz de registro de datos
    M_log = []; 
    # Recorro M_csv
    for i in range(len(M_csv)):
        curr_row = M_csv[i]; 
        # Solo agarro curr_row si tiene info de gen cercano
        if len(curr_row) > min_row_len:
            # Defino variables para PCA
            curr_seq, seq_usada = procesar_seq(curr_row[seq_col]); 
            # Agrego las variables PCA a M_pca_in
            M_pca_in.append(curr_seq[:]); 
            # Agrego las variables de PCA y de curr_row a M_out
            M_log.append(curr_row + [seq_usada] + curr_seq); 
    return M_pca_in, M_log


def procesar_seq(seq, len_target=7):
    # Funcion para transformar secuencia en lista de 0 y 1 para PCA
    # Incluye UN MONTON de cosas hardcodeadas
    # Pensado para len_target=7

    # Veo si seq tiene largo igual a len_target
    if len(seq) == len_target:
        # Veo si termina en GTG o empieza en CAC
        if seq[:3]=='CAC':
            seq_usada = invertir_seq(seq); 
        elif seq[-3:]=='GTG':
            seq_usada = seq; 
        else:
            print('WARNING: Secuencia ' + str(seq) + ' no se pudo parsear.')
    # Si no tiene largo igual a len_target, hago busqueda hardcodeada de largos 8 y 10
    elif len(seq) == 8:
        # busco GTG/GAG o CAC/CTC en las puntas
        if seq[1:4]=='CAC':
            seq_usada = invertir_seq(seq[1:]); 
        elif seq[-4:-1]=='GTG':
            seq_usada = seq[:-1]; 
        elif seq[1:4]=='CTC':
            seq_usada = invertir_seq(seq[1:]); 
        elif seq[-4:-1]=='GAG':
            seq_usada = seq[:-1]; 
        else:
            print('WARNING: Secuencia ' + str(seq) + ' no se pudo parsear.')
    elif len(seq) == 10:
        # busco GTG/GAG o CAC/CTC en las puntas
        if seq[1:4]=='CAC':
            seq_usada = invertir_seq(seq[1:-2]); 
        elif seq[-4:-1]=='GTG':
            seq_usada = seq[2:-1]; 
        elif seq[1:4]=='CTC':
            seq_usada = invertir_seq(seq[1:-2]); 
        elif seq[-4:-1]=='GAG':
            seq_usada = seq[2:-1]; 
        else:
            print('WARNING: Secuencia ' + str(seq) + ' no se pudo parsear.')
    # Inicializo la lista de 0/1 que se devuelve
    L_out = []; 
    # Defino diccionario para pasar de nucleotidos a 0 y 1
    dict_bin = {'A':[1,0,0,0], 'C':[0,1,0,0], 'T':[0,0,1,0], 'G':[0,0,0,1]}; 
    # Recorro seq_usada
    for n in range(len(seq_usada)):
        curr_c = seq_usada[n]; 
        # Paso cada caracter a lista de 0 y 1 con dict_bin
        L_out = L_out + dict_bin[curr_c]; 
    return L_out, seq_usada


def procesar_sitios(M_sitios, genoma, nom_genoma, dist_max=500000, path_fasta='', sub_sep=','):
    # Revisa sitios de union y agrega informacion relevante que voy a necesitar para PCA

    # Inicializo el elemento seq_data que consigue las secuencias
    sequence_data = seq_data(nom_genoma, genome_element=genoma, path_fasta=path_fasta); 
    # Inicializo la matriz que se devuelve
    M_out = []; 
    ###
    largo_M_sitios = len(M_sitios); 
    ###
    # Recorro M_sitios
    for i in range(largo_M_sitios):
        curr_sitio = M_sitios[i]; 
        # Defino contig, pos_ini y pos_end del sitio
        contig = curr_sitio[0]; 
        L_pos = curr_sitio[1].split('_'); 
        L_pos[0] = int(L_pos[0]); 
        L_pos[1] = int(L_pos[1]); 
        pos_ini = min(L_pos); 
        pos_end = max(L_pos); 
        # Defino secuencia
        if curr_sitio[2] != '':
            # Si ya hay una secuencia, me ahorro volver a buscarla
            seq = curr_sitio[2]; 
        else:
            # Consigo la secuencia del peak con seq_data._consulta_secuencia_fasta()
            seq = sequence_data._consulta_secuencia_fasta('chr'+str(contig), pos_ini, pos_end); 
        # Defino genes cercanos 
        L_genes_cerca_raw = genoma.genes_at_locus(contig, pos_ini-dist_max, pos_end+dist_max); 
        # Defino lista curada de genes cerca en formato string
        str_genes_cerca = ''; 
        # Recorro todos los genes encontrados
        for curr_gen in L_genes_cerca_raw:
            # Defino gene_start en base a strand
            if curr_gen.strand=='+':
                gen_start = curr_gen.start; 
                # Reviso que gen_start este entre pos_ini-dist_max y pos_end+dist_max
                registro_gen = (gen_start>(pos_ini-dist_max)) and (gen_start<(pos_end+dist_max)); 
            elif curr_gen.strand=='-':
                gen_start = curr_gen.end; 
                # Reviso que gen_start este entre pos_ini-dist_max y pos_end+dist_max
                registro_gen = (gen_start>(pos_ini-dist_max)) and (gen_start<(pos_end+dist_max)); 
            else:
                print('WARNING: Strand no se pudo procesar en gen ' + str(curr_gen))
                # Si strand da error, registro_gen queda como True, por las dudas
                registro_gen = True; 
            # Solo registro los genes con registro_gen=True y biotype protein_coding
            if registro_gen and (curr_gen.biotype == 'protein_coding'):
                # Uso sub_sep para separar los gene_id
                str_genes_cerca = str_genes_cerca + curr_gen.gene_id + sub_sep; 
        # Elimino la ultima ocurrencia de sub_sep
        str_genes_cerca = str_genes_cerca.rstrip(sub_sep); 
        # Agrego la informacion relevante a M_out en una lista
        M_out.append(['chr'+str(contig), int(pos_ini), int(pos_end), str(seq), str(str_genes_cerca)]); 
        ### Display
        if (i+1)%400 == 0:
            print('PROGRESO: ' + str(i+1) + '/' + str(largo_M_sitios))
        ###
    return M_out


def simplificar_dist(dist):
    # Funcion para transformar un valor numerico de dist en un valor simplificado
    # Funcion hardcodeada para funcionar con dist_max=500k
    # Asume que se pueden usar valores negativos

    # Pruebo de a rangos
    if dist < -500000:
        r = -5; 
    elif dist < -100000:
        r = -4; 
    elif dist < -50000:
        r = -3; 
    elif dist < -10000:
        r = -2; 
    elif dist < 0:
        r = -1; 
    elif dist == 0:
        r = 0; 
    elif dist <= 10000:
        r = 1; 
    elif dist <= 50000:
        r = 2; 
    elif dist <= 100000:
        r = 3; 
    elif dist <= 500000:
        r = 4; 
    else:
        r = 5; 
    return r


def std_len_head(L_head, head_len, ext=''):
    # Funcion para estandarizar el largo de una lista

    # Comparo el largo de L_head con head_len
    if len(L_head) < head_len: 
        # Si L_head es mas corto, agrego ext a la lista que se devuelve hasta que tenga el mismo largo
        L_out = L_head + [str(ext)]*(head_len-len(L_head)); 
    elif len(L_head) > head_len: 
        # Si L_head es mas largo, tiro warning y corto el largo de mas
        print('WARNING: L_head ' + str(L_head) + ' es mas largo que head_len. Se corta el largo de mas.')
        L_out = L_head[:head_len]; 
    else: 
        # Si tienen el mismo largo, solo devuelvo L_out
        L_out = L_head; 
    return L_out


def std_len_M_pca(M_pca):
    # Funcion que estandariza el largo de cada elemento de M_pca

    # Inicializo la matriz que se devuelve
    M_out = []; 
    # Inicializo el valor de largo maximo
    max_len = 0; 
    # Recorro M_pca para determinar el largo maximo de los elementos adentro
    for i in range(len(M_pca)):
        curr_len = len(M_pca[i]); 
        # Si el largo del elemento actual es mayor a max_len, lo anoto
        if curr_len>max_len:
            max_len = int(curr_len); 
    # Recorro M_pca una segunda vez para hacer que todos los elementos tengan el mismo largo
    for j in range(len(M_pca)):
        curr_len = len(M_pca[j]); 
        # Si el largo del elemento actual es menor a max_len, le agrego ceros al final
        if curr_len<max_len:
            len_dif = max_len-curr_len; 
            M_out.append(M_pca[j]+[0]*len_dif); 
        else:
            M_out.append(M_pca[j]); 
    return M_out


#################################### RESULTADOS ###################################

output_dump = []; 

if __name__=='__main__':
    output_dump.append(_main()); 

