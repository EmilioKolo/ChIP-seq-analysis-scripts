
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

# Generacion de graficos
import matplotlib.pyplot as plt

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
Funciones para reorganizar lo que hice en 11-Graficos_GO.py (ver ese archivo para generacion de archivos para generar los diccionarios de referencia)

Incluye generacion de graficos con datos descargados a mano de Biomart
'''

#################################### VARIABLES ####################################

# Direcciones de archivos y outputs para ib3 y casa
path_dropbox_casa = 'D:\\Users\\Admin\\Dropbox\\'; 
path_dropbox_ib3 = 'C:\\Users\\emili\\Dropbox\\'; 
path_fasta_casa = 'D:\\Archivos doctorado\\Genomas\\'; 
path_fasta_ib3 = 'X:\\Genomas\\'; 
path_output_dump_casa = 'D:\\Archivos doctorado\\Output_dump\\'; 
path_output_dump_ib3 = 'X:\\Output_dump\\'; 

# Nombres de archivos de input
nom_csv_human = '12-AndersonGO'; 
nom_csv_mouse = '12-DupaysGO'; 

# Nombres de archivos de peak usados para filtro
nom_peaks_human = '13-Anderson100k'; 
nom_peaks_mouse = '13-Dupays100k'; 

# Filtros usados em formato lista de tuplas (columna, valor esperado)
L_filtros_con_su = [(4, '1-Regulado directamente')]; 
L_filtros_sin_su = [(4, '2-Regulado indirectamente')]; 

col_genes_main = 17; 
n_bars_main=10; 


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
path_rnaseq_main = path_dropbox_main + 'Doctorado\\3-Genes transactivados rio abajo\\1-Interpretacion de datos ChIP-seq\\'; 
path_in_main = path_dropbox_main + 'Doctorado\\3-Genes transactivados rio abajo\\1-Interpretacion de datos ChIP-seq\\gene_ontology\\'; 
path_out_main = path_output_dump_main + 'GOterms\\'; 
path_hist_main = path_out_main + 'Graficos\\'
# Path que dependen de la especie
path_pwm_human = path_dropbox_main + 'Doctorado\\3-Genes transactivados rio abajo\\1-Interpretacion de datos ChIP-seq\\PWM_human\\'; 
path_pwm_mouse = path_dropbox_main + 'Doctorado\\3-Genes transactivados rio abajo\\1-Interpretacion de datos ChIP-seq\\PWM_mouse\\'; 



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


def abrir_dict(nom_arch, n=1, path_arch='', ext='.csv', sep=';', sep2=''):
    # Funcion para devolver un diccionario extraido de un .csv

    # Inicializo el diccionario que se devuelve
    dict_out = {}; 
    # Abro el csv como matriz
    M_arch = abrir_csv(nom_arch, dir_arch=path_arch, ext=ext, sep=sep, ignore_headers=False); 
    # Recorro M_arch
    for i in range(len(M_arch)):
        curr_L = M_arch[i]; 
        # Agrego curr_L[0] a dict_out, si es que no existe
        if not (curr_L[0] in dict_out.keys()):
            dict_out[curr_L[0]] = []; 
        else:
            print('ERROR: Mas de una fila asociada a ' + str(curr_L[0]) + ' en ' + str(nom_arch))
        # Si curr_L tiene uno o dos elementos, no recorro nada
        if len(curr_L)<=2:
            # Si curr_L es igual a 1, dejo dict_out como lista vacia
            if len(curr_L)>1:
                # Si sep2 tiene un valor, lo uso para separar curr_L[j]
                if sep2 != '':
                    curr_val = curr_L[1].split(sep=sep2); 
                else:
                    curr_val = [curr_L[1]]; 
                dict_out[curr_L[0]] = dict_out[curr_L[0]] + curr_val; 
        # Si curr_L tiene mas de dos elementos, lo recorro
        else:
            for j in range(1, len(curr_L)):
                # Si sep2 tiene un valor, lo uso para separar curr_L[j]
                if sep2 != '':
                    curr_val = curr_L[j].split(sep=sep2); 
                else:
                    curr_val = curr_L[j]; 
                # Si n es mayor que 1, inicializo sub-listas
                if n>1:
                    if j%n==0:
                        dict_out[curr_L[0]].append([]); 
                    # Agrego curr_val a dict_out[curr_L[0]][j//n]
                    dict_out[curr_L[0]][j//n].append(curr_val); 
                # Si n es igual o menor a 1, agrego los valores directamente a dict_out[curr_L[0]]
                else:
                    dict_out[curr_L[0]].append(curr_val); 
    return dict_out


def cargar_dict_go(nom_csv, path_csv=''):
    # Funcion que agarra los csv descargados manualmente de biomart y devuelve diccionarios de frecuencias de GO terms para histogramas

    # Abro el csv
    M_csv = abrir_csv(nom_csv, dir_arch=path_csv, sep='\t', ignore_headers=True, ext='.txt'); 
    # Inicializo el diccionario que se devuelve
    dict_freq = {}; 
    # Recorro M_csv
    for i in range(len(M_csv)): 
        curr_term = M_csv[i]; 
        # Reviso que el GO-term sea biological_process
        if len(curr_term) >= 4 and (curr_term[3]=='biological_process' and curr_term[2]!='biological_process'):
            # Reviso si el nombre del GO-term ya esta registrado en dict_freq
            if curr_term[2] in dict_freq.keys():
                # Si ya esta presente, sumo uno
                dict_freq[curr_term[2]] += 1; 
            else:
                # Si no esta presente, lo agrego, inicializado en 1
                dict_freq[curr_term[2]] = 1; 
    return dict_freq


def filtrar_dict_go(nom_dict_go, nom_dict_go_out, nom_arch_filtro, col_genes, filtro_usado=[], path_arch_filtro='', path_dict_go='', sep_genes=','):
    # Funcion que recibe el nombre de un archivo que recibiria cargar_dict_go() y el nombre de un archivo de peaks clasificados
    # Filtro usado contiene tuplas de formato (columna, valor esperado) para revisar en nom_arch_filtro
    # Devuelve un output como el de cargar_dict_go con los filtros aplicados
    # Devuelve diccionarios de frecuencias de GO terms para histogramas
    
    # Abro nom_arch_filtro para cargar los peaks en M_csv_filtro
    M_csv_filtro = abrir_csv(nom_arch_filtro, dir_arch=path_arch_filtro, ext='.csv', sep=';', ignore_headers=True); 
    # Inicializo lista de genes seleccionados
    L_gene_id_sel = []; 
    # Recorro M_csv_filtro
    for i in range(len(M_csv_filtro)):
        curr_row = M_csv_filtro[i]; 
        # Inicializo booleano para ver si me quedo con la fila
        keep_row = True; 
        # Recorro filtro_usado
        for curr_filtro in filtro_usado:
            # curr_row[curr_filtro[0]]==curr_filtro[1] tiene que valer para todos los curr_filtro
            if curr_row[curr_filtro[0]]!=curr_filtro[1]:
                # Si cualquiera de los filtros falla, descarto la fila
                keep_row = False; 
        # Veo si keep_row sigue siendo True
        if keep_row:
            # Reviso los genes en col_genes
            L_curr_genes = curr_row[col_genes].split(sep_genes); 
            # Recorro cada gen en L_curr_genes
            for curr_gen in L_curr_genes:
                # Veo si curr_gen ya esta en la lista de gene_id seleccionados
                if not (str(curr_gen) in L_gene_id_sel):
                    # Si no esta, lo agrego con append
                    L_gene_id_sel.append(str(curr_gen)); 
    # Abro nom_dict_go
    M_csv_dict_go = abrir_csv(nom_dict_go, dir_arch=path_dict_go, ext='.txt', sep='\t', ignore_headers=True); 
    # Inicializo la matriz que se usa para generar el diccionario de histogramas
    M_csv_dict_go_out = []; 
    # Recorro M_csv_dict_go
    for j in range(len(M_csv_dict_go)):
        curr_row_dict_go = M_csv_dict_go[j]; 
        # Veo que curr_row_dict_go[0] este en L_gene_id_sel
        if curr_row_dict_go[0] in L_gene_id_sel:
            # Si esta, agrego curr_row_dict_go entero a M_csv_dict_go_out
            M_csv_dict_go_out.append(M_csv_dict_go[j]); 
    ## De aca en adelante corro una copia de cargar_dict_go()
    # Inicializo el diccionario que se devuelve
    dict_freq = {}; 
    # Recorro M_csv
    for k in range(len(M_csv_dict_go_out)): 
        curr_term = M_csv_dict_go_out[k]; 
        # Reviso que el GO-term sea biological_process
        if len(curr_term) >= 4 and (curr_term[3]=='biological_process' and curr_term[2]!='biological_process'):
            # Reviso si el nombre del GO-term ya esta registrado en dict_freq
            if curr_term[2] in dict_freq.keys():
                # Si ya esta presente, sumo uno
                dict_freq[curr_term[2]] += 1; 
            else:
                # Si no esta presente, lo agrego, inicializado en 1
                dict_freq[curr_term[2]] = 1; 
    return dict_freq


def generar_M_hist(dict_gene_go, dict_go={}, id_term=2, mayor_a_menor=True):
    # Funcion que recibe un diccionario de GO asociados a genes y un diccionario de GO ID a nombre (opcional)
    # Devuelve una lista de listas con GO terms y cantidad de veces que aparecen en dict_gene_go, ordenada de mayor a menor
    # Si no se da dict_go, la lista de listas es con GO ID en vez de GO term
    # Si mayor_a_menor es False, se devuelve M_out de menor a mayor

    # Inicializo la matriz que se devuelve
    M_out = []; 
    # Inicializo un diccionario para hacer el conteo
    dict_cont = {}; 
    # Defino si se usa dict_go o no
    if dict_go == {}:
        usar_go_terms = False; 
    else:
        usar_go_terms = True; 
    # Recorro dict_gene_go.keys()
    dict_gene_go_keys = list(dict_gene_go.keys()); 
    for g in range(len(dict_gene_go_keys)):
        curr_gene = dict_gene_go_keys[g]; 
        curr_L_go = dict_gene_go[curr_gene]; 
        # Recorro cada uno de los GO terms asociados al gen
        for curr_go_id in curr_L_go:
            # Veo si curr_go_id esta en dict_cont.keys()
            if curr_go_id in dict_cont.keys():
                # Agrego uno al contador 
                dict_cont[curr_go_id] += 1; 
            else:
                # Inicializo la key en 1
                dict_cont[curr_go_id] = 1; 
    # Recorro dict_cont una vez cargados todos los valores
    dict_cont_keys = list(dict_cont.keys()); 
    for h in range(len(dict_cont_keys)):
        curr_go_id = dict_cont_keys[h]; 
        # Inicializo la lista que se agrega a M_out
        L_append = []; 
        # Veo si uso GO term o GO ID
        if usar_go_terms:
            L_append.append(str(dict_go[curr_go_id][id_term])); 
        else:
            L_append.append(str(curr_go_id)); 
        # Agrego el conteo
        L_append.append(int(dict_cont[curr_go_id])); 
        # Agrego L_append a M_out
        M_out.append(L_append[:]); 
    # Ordeno M_out
    M_out.sort(key=lambda x: x[1], reverse=mayor_a_menor); 
    return M_out


def guardar_dict(dict_out, nom_out, path_out='.\\', ext='.csv', sep=';', sep2=','):
    # Guarda los contenidos del diccionario dict_out en el archivo nom_out ubicado en path_out

    # Defino la direccion del archivo en base a path_out y nom_out
    dirarch = os.path.join(path_out, nom_out + ext); 
    # Creo el archivo
    with open(dirarch, 'w') as F_out:
        print('Archivo ' + nom_out + ext + ' creado.')
    # Lo vuelvo a abrir en modo append
    with open(dirarch, 'a') as F_out:
        # Recorro dict_out.keys()
        for curr_key in dict_out.keys():
            # dict_out[key] contiene una lista con strings
            curr_L = dict_out[curr_key]; 
            # Pruebo a ver si curr_L es una lista
            try:
                # Transformo curr_L en un string separado por sep2
                curr_str = ''; 
                for i in curr_L:
                    curr_str = curr_str + str(i) + sep2; 
                curr_str = curr_str.rstrip(sep2); 
            except TypeError:
                curr_str = str(curr_L); 
            # Agrego curr_key + sep + curr_str a F_out
            F_out.write(curr_key + sep + curr_str + '\n'); 
    return dict_out


def histograma_dict(dict_hist, nom_hist, guardar=True, path_hist='.\\', n_bars=10):
    # Funcion para generar un histograma en base a un diccionario con los valores del eje x en keys y sus frecuencias como valores del diccionario

    # Genero listas para el histograma
    keys_hist = list(dict_hist.keys()); 
    values_hist = dict_hist.values(); 
    # Ordeno las listas usadas por el histograma de mayor a menor (de acuerdo a la frecuencia asociada a cada key)
    values_hist, keys_hist = zip(*sorted(zip(values_hist, keys_hist), reverse=True)); 
    # Creo el histograma
    plt.figure(figsize=(6,4)); 
    n = plt.bar(keys_hist[:n_bars], values_hist[:n_bars], color='#040499', width=0.9); 
    plt.grid(axis='y', alpha=0.6); 
    plt.xticks(rotation=15, ha='right'); 
    plt.xlabel('GO term'); 
    plt.ylabel('Frecuencia'); 
    plt.title(nom_hist); 
    if guardar:
        plt.savefig(path_hist + 'Histograma_' + str(nom_hist) + '.png', dpi=300); 
    else:
        plt.show(); 
    plt.close(); 
    return n


def _main():
    # Funcion para probar funciones en ejecucion del archivo

    '''# Abro los diccionarios de GO terms
    dict_go_human = abrir_dict('GOid_info_human', path_arch=path_out_main, sep2=','); 
    dict_go_mouse = abrir_dict('GOid_info_mouse', path_arch=path_out_main, sep2=','); 

    # Abro los diccionarios de GO terms asociados a cada gen
    dict_gene_go_human = abrir_dict('HumanGenesGOdict', path_arch=path_in_main, sep2=','); 
    dict_gene_go_mouse = abrir_dict('MouseGenesGOdict', path_arch=path_in_main, sep2=','); 

    # Cuento las ocurrencias de GO terms asociados a todos los genes en formato para histograma
    M_hist_human = generar_M_hist(dict_gene_go_human, dict_go=dict_go_human, id_term=2); 
    M_hist_mouse = generar_M_hist(dict_gene_go_mouse, dict_go=dict_go_mouse, id_term=2); '''

    # Defino los titulos de los histogramas
    nom_hist_human = 'hist_human'; 
    nom_hist_mouse = 'hist_mouse'; 
    # Cargo los diccionarios desde los archivos
    dict_hist_human = cargar_dict_go(nom_csv_human, path_csv=path_in_main); 
    dict_hist_mouse = cargar_dict_go(nom_csv_mouse, path_csv=path_in_main); 
    # Filtro por genes cerca de sitios confirmados por RNAseq con o sin sitio de union de NKX2-5
    dict_hist_human_con_su = filtrar_dict_go(nom_csv_human, 'human_con_su', nom_peaks_human, col_genes_main, filtro_usado=L_filtros_con_su, path_arch_filtro=path_rnaseq_main, path_dict_go=path_in_main); 
    dict_hist_human_sin_su = filtrar_dict_go(nom_csv_human, 'human_sin_su', nom_peaks_human, col_genes_main, filtro_usado=L_filtros_sin_su, path_arch_filtro=path_rnaseq_main, path_dict_go=path_in_main); 
    dict_hist_mouse_con_su = filtrar_dict_go(nom_csv_mouse, 'mouse_con_su', nom_peaks_mouse, col_genes_main, filtro_usado=L_filtros_con_su, path_arch_filtro=path_rnaseq_main, path_dict_go=path_in_main); 
    dict_hist_mouse_sin_su = filtrar_dict_go(nom_csv_mouse, 'mouse_sin_su', nom_peaks_mouse, col_genes_main, filtro_usado=L_filtros_sin_su, path_arch_filtro=path_rnaseq_main, path_dict_go=path_in_main); 
    # Genero los histogramas con los diccionarios
    #out_human = histograma_dict(dict_hist_human, nom_hist_human, guardar=True, path_hist=path_hist_main, n_bars=n_bars_main); 
    #out_mouse = histograma_dict(dict_hist_mouse, nom_hist_mouse, guardar=True, path_hist=path_hist_main, n_bars=n_bars_main); 

    # Guardo los diccionarios
    dict_hist_human_con_su = guardar_dict(dict_hist_human_con_su, 'dict_hist_human_con_su', path_out=path_out_main, sep2=';'); 
    dict_hist_human_sin_su = guardar_dict(dict_hist_human_sin_su, 'dict_hist_human_sin_su', path_out=path_out_main, sep2=';'); 
    dict_hist_mouse_con_su = guardar_dict(dict_hist_mouse_con_su, 'dict_hist_mouse_con_su', path_out=path_out_main, sep2=';'); 
    dict_hist_mouse_sin_su = guardar_dict(dict_hist_mouse_sin_su, 'dict_hist_mouse_sin_su', path_out=path_out_main, sep2=';'); 
    # Creo los histogramas
    #out_human_con_su = histograma_dict(dict_hist_human_con_su, nom_hist_human+'_con_su', guardar=True, path_hist=path_hist_main, n_bars=n_bars_main); 
    #out_human_sin_su = histograma_dict(dict_hist_human_sin_su, nom_hist_human+'_sin_su', guardar=True, path_hist=path_hist_main, n_bars=n_bars_main); 
    #out_mouse_con_su = histograma_dict(dict_hist_mouse_con_su, nom_hist_mouse+'_con_su', guardar=True, path_hist=path_hist_main, n_bars=n_bars_main); 
    #out_mouse_sin_su = histograma_dict(dict_hist_mouse_sin_su, nom_hist_mouse+'_sin_su', guardar=True, path_hist=path_hist_main, n_bars=n_bars_main); 
    return ''


#################################### RESULTADOS ###################################

output_dump = []; 

if __name__=='__main__':
    output_dump.append(_main()); 

