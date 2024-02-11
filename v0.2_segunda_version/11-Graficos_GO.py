
ib3 = False; 
http_proxy_ib3 = 'http://proxy.fcen.uba.ar:8080'; 
https_proxy_ib3 = 'https://proxy.fcen.uba.ar:8080'; 
ftp_proxy_ib3 = 'ftp://proxy.fcen.uba.ar:8080'; 

# Generales
import copy
import os
import sys
import time
from random import shuffle
from pyensembl import EnsemblRelease
import numpy as np

# Analisis de secuencias y genomas
from Bio import Entrez, SeqIO, motifs, Seq
Entrez.email = 'ekolomenski@gmail.com'; 
Entrez.api_key = '9bccbd4949f226a5930868df36c211e8b408'; 
import biomart

# Generacion de graficos
import matplotlib.pyplot as plt

# Interaccion con REST API
import requests

# Importo analisischip
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
Scripts para trabajar con los GO id conseguidos para genes asociados a los resultados de sitios de union sacados de ChIP-seq y confirmados por RNA-seq
'''

#################################### VARIABLES ####################################

# Direcciones de archivos y outputs para ib3 y casa
path_dropbox_casa = 'D:\\Users\\Admin\\Dropbox\\'; 
path_dropbox_ib3 = 'C:\\Users\\emili\\Dropbox\\'; 
path_fasta_casa = 'D:\\Archivos doctorado\\Genomas\\'; 
path_fasta_ib3 = 'X:\\Genomas\\'; 
path_output_dump_casa = 'D:\\Archivos doctorado\\Output_dump\\'; 
path_output_dump_ib3 = 'X:\\Output_dump\\'; 

# Nombres de genomas usados
genome_name_mouse = 'mm9'; 
genome_name_human = 'hg19'; 

# Nombres de csv con peaks
#nom_csv_human = '13-Anderson100k'; 
#nom_csv_mouse = '13-Dupays100k'; 
nom_csv_human = '14-Anderson10k'; 
nom_csv_mouse = '14-Dupays10k'; 



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
path_out_graficos = path_out_main + 'Graficos\\'; 
# Path que dependen de la especie
path_pwm_human = path_dropbox_main + 'Doctorado\\3-Genes transactivados rio abajo\\1-Interpretacion de datos ChIP-seq\\PWM_human\\'; 
path_pwm_mouse = path_dropbox_main + 'Doctorado\\3-Genes transactivados rio abajo\\1-Interpretacion de datos ChIP-seq\\PWM_mouse\\'; 

col_class_main = 4; 
#col_genes_main = 17; 
#nom_mod_main = '_100k'; 
col_genes_main = 15; 
nom_mod_main = '_10k'; 


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


def abrir_csv_headers(nom_arch, path_arch='', ext='.csv', sep=';'):
    # Usa abrir_csv() y abrir_headers() para devolver una matriz extraida de un .csv y sus headers

    # Genero matriz csv sin headers
    M_csv = abrir_csv(nom_arch, dir_arch=path_arch, ext=ext, sep=sep, ignore_headers=True); 
    # Genero lista headers
    L_headers = abrir_headers(nom_arch, dir_arch=path_arch, ext=ext, sep=sep); 
    return M_csv, L_headers


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


def abrir_headers(nom_arch, dir_arch='', ext='.csv', sep=';'):
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


def dividir_lista(L, n):
    # Divide L en pedazos de largo n

    # Inicializo la lista que se devuelve
    L_out = []
    # Recorro L en saltos de a n
    for i in range(0, len(L), n):
        L_out.append(L[i:i+n])
    return L_out


def generar_data_used(server, dataset, attributes, decode='ascii'):
    # Funcion que recibe la variable server de biomart y devuelve data_used

    # Agarro el dataset dentro de server
    mart = server.datasets[dataset]; 
    # Busco attributes y los guardo en response
    response = mart.search({'attributes':attributes}); 
    # Saco la info de response y la guardo en data_used
    data_used = response.raw.data.decode(decode); 
    return data_used


def generar_dict_go(server, genome_name, verbose=True):
    # 

    # Atributos GO: go_id, go_linkage_type, goslim_goa_accession, goslim_goa_description
    attributes = ['go_id', 'go_linkage_type', 'goslim_goa_accession', 'goslim_goa_description']; 
    dataset = genome_to_dataset(genome_name); 
    if verbose:
        print('Iniciando descarga de diccionario para ' + str(dataset) + '.')
    data_used = generar_data_used(server, dataset, attributes, decode='ascii'); 
    if verbose:
        print('data_used obtenido. Iniciando procesamiento de los datos.')
    # Inicializo el diccionario que se devuelve
    dict_out = {}; 
    # Recorro data_used
    for curr_line in data_used.splitlines(): 
        # Uso el mismo orden que la lista de atributos (attributes_X)
        line = curr_line.split('\t'); 
        print(line)
        # Veo si line[0] (go_id) esta en dict_out.keys()
        if not (line[0] in dict_out.keys()):
            # Si no esta, agrego line[1], line[2] y line[3] a dict_out[line[0]]
            dict_out[line[0]] = line[1:]; 
        # Si line[0] esta en dict_out.keys(), veo que line[1:] este dentro de la lista en dict_out[line[0]]
        elif not (line[1:] in dividir_lista(dict_out[line[0]], len(line[1:]))):
            print('ERROR: ' + str(line[0]) + ' aparece mas de dos veces en dict_out con distintos valores.')
            dict_out[line[0]] = dict_out[line[0]] + line[1:]; 
    return dict_out


def generar_dict_go_api(L_ids_go, verbose_range=10):
    # Funcion para generar un diccionario con informacion sobre una lista de IDs de gene ontology

    # Inicializo el diccionario que se devuelve
    dict_out = {}; 
    ### Display
    len_ids = len(L_ids_go); 
    ###
    # Recorro L_ids_go
    for k in range(len(L_ids_go)):
        curr_id_go = L_ids_go[k]; 
        # Extraigo la info de curr_id_go
        L_info_id_go = traducir_go_id_dict(curr_id_go); 
        # Veo si curr_id_go esta entre las keys de dict_out
        if not (curr_id_go in dict_out.keys()):
            dict_out[curr_id_go] = []; 
        # Agrego la info extraida de curr_id_go
        dict_out[curr_id_go] = dict_out[curr_id_go] + L_info_id_go; 
        ### Display
        if (k+1)%verbose_range == 0:
            print('Progreso: ' + str(k+1) + '/' + str(len_ids))
        ###
    return dict_out


def generar_histograma(dict_hist, nombre_histograma, guardar=True, n_bars=10, path_hist='.\\'):
    # 

    # Genero listas para el histograma
    keys_hist = list(dict_hist.keys()); 
    print('Iniciando chequeo de GO terms.')
    len_keys = len(keys_hist); 
    for k in range(len_keys):
        keys_hist[k] = traducir_go_id(keys_hist[k]); 
        if (k+1)%10 == 0:
            print('Progreso: ' + str(k+1) + '/' + str(len_keys))
    values_hist = dict_hist.values(); 
    values_hist, keys_hist = zip(*sorted(zip(values_hist, keys_hist), reverse=True)); 
    # Creo el histograma
    plt.figure(figsize=(18,12)); 
    n = plt.bar(keys_hist[:n_bars], values_hist[:n_bars], color='#040499', width=0.9); 
    plt.grid(axis='y', alpha=0.6); 
    plt.xlabel('GO term'); 
    plt.ylabel('Frecuencia'); 
    plt.title(nombre_histograma); 
    if guardar:
        plt.savefig(path_hist + 'Histograma_' + str(nombre_histograma) + '.png', dpi=300); 
    else:
        plt.show(); 
    plt.close(); 
    return n


def genome_to_dataset(genome_name, modifier='_gene_ensembl'):
    # Defino dataset_name para usar con biomart 
    if genome_name.lower() == 'mm9' or genome_name.lower() == 'mm10' or 'mmusculus' in genome_name.lower():
        dataset_name = 'mmusculus' + modifier; 
    elif genome_name.lower() == 'hg19' or genome_name.lower() == 'hg38' or 'hsapiens' in genome_name.lower():
        dataset_name = 'hsapiens' + modifier; 
    else:
        print('ERROR: Nombre de genoma ' + str(genome_name) + ' desconocido. Usando genoma de raton.')
        dataset_name = 'mmusculus' + modifier; 
    return dataset_name


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
            # Transformo curr_L en un string separado por sep2
            curr_str = ''; 
            for i in curr_L:
                curr_str = curr_str + str(i) + sep2; 
            curr_str = curr_str.rstrip(sep2); 
            # Agrego curr_key + sep + curr_str a F_out
            F_out.write(curr_key + sep + curr_str + '\n'); 
    return dict_out


def guardar_lista_genes(nom_out, L_genes, sep=',', path_out='.\\', ext='.txt'): 
    # Guarda una lista de genes L_genes como archivo .txt con los elementos separados por sep

    # Defino la direccion del archivo en base a path_out y nom_out
    dirarch = os.path.join(path_out, nom_out + ext); 
    # Defino str_genes en base a L_genes y sep
    str_genes = ''; 
    # Recorro L_genes
    for curr_gen in L_genes:
        str_genes = str_genes + str(curr_gen) + sep; 
    # Elimino la ultima ocurrencia de sep
    str_genes = str_genes.rstrip(sep); 
    # Creo el archivo o lo sobreescribo con un archivo vacio si ya existe
    with open(dirarch, 'w') as F_out: 
        print('Archivo ' + nom_out + ext + ' creado.')
    # Lo vuelvo a abrir en modo append
    with open(dirarch, 'a') as F_out: 
        F_out.write(str_genes); 
    return L_genes


def pipeline_generar_dict_go_api(nom_in, nom_out, path_in='', path_out='.\\'):
    # Funcion que usa generar_dict_go_api() y guardar_dict() para generar un archivo con el diccionario correspondiente

    # Abro el diccionario con los GO-terms asociados a genes usando abrir_dict()
    dict_genes = abrir_dict(nom_in, path_arch=path_in, sep2=','); 
    # Genero L_ids_go con dict_genes
    L_ids_go = []; 
    # Recorro cada una de las keys de dict_genes
    for curr_key in dict_genes.keys(): 
        curr_L = dict_genes[curr_key]; 
        # Reviso cada uno de los elementos en curr_L
        for k in range(len(curr_L)):
            curr_val = curr_L[k]; 
            # Reviso si curr_val ya esta en L_ids_go y es distinto de ''
            if (not (curr_val in L_ids_go)) and curr_val!='':
                L_ids_go.append(curr_val); 
    # Armo el diccionario usando L_ids_go
    dict_go_api = generar_dict_go_api(L_ids_go); 
    # Guardo el diccionario usando guardar_dict()
    dict_go_api = guardar_dict(dict_go_api, nom_out, path_out=path_out); 
    return dict_go_api


def pipeline_histograma(nom_dict_genes, nom_hist, bars_hist=10, guardar_hist=True, path_dict_genes='', path_hist='.\\'):
    # 

    # Inicializo el histograma como diccionario
    hist_out = {}; 
    # Abro el diccionario con los GO-terms asociados a genes usando abrir_dict()
    dict_genes = abrir_dict(nom_dict_genes, path_arch=path_dict_genes, sep2=','); 
    # Recorro cada una de las keys de dict_genes
    for curr_key in dict_genes.keys(): 
        curr_L = dict_genes[curr_key]; 
        # Reviso cada uno de los elementos en curr_L
        for k in range(len(curr_L)):
            curr_val = curr_L[k]; 
            # Si curr_val no esta en hist_out, agrego la key
            if not (curr_val in hist_out.keys()):
                hist_out[curr_val] = 0; 
            # Agrego uno a hist_out[curr_val]
            hist_out[curr_val] += 1; 
    # Genero el histograma usando hist_out
    output_plt = generar_histograma(hist_out, nom_hist, guardar=guardar_hist, n_bars=bars_hist, path_hist=path_hist); 
    return hist_out


def pipeline_lista_genes(nom_peaks, col_genes, L_filtros=[], sep_genes=',', path_peaks='', nom_out='', path_out='.\\'):
    # Funcion que devuelve una lista de genes siguiendo un filtro
    # L_filtros tiene tuplas de dos elementos con posicion y valor esperado. e.g. [(1, 'X'), (4, '')]
        # con_su: [(4, '1-Regulado directamente')] // sin_su: [(4, '2-Regulado indirectamente')]

    # Inicializo la lista de genes que se devuelve
    L_genes = []; 
    # Abro la tabla de peaks y los guardo en una matriz 
    M_peaks = abrir_csv(nom_peaks, dir_arch=path_peaks, ignore_headers=True); 
    # Recorro M_peaks
    for i in range(len(M_peaks)):
        curr_peak = M_peaks[i]; 
        # Booleano para definir si el pico se usa o no
        bool_filtro = True; 
        # Recorro L_filtros y veo que se cumplan todos
        for curr_filtro in L_filtros:
            # Reviso que el filtro tenga por lo menos 2 elementos
            if len(curr_filtro)>=2:
                # Veo que curr_peak[curr_filtro[0]] sea igual a curr_filtro[1]
                if curr_peak[curr_filtro[0]]!=curr_filtro[1]:
                    bool_filtro = False; 
            else:
                print('WARNING: Filtro ' + str(curr_filtro) + ' tiene menos de 2 elementos.')
        # Veo que bool_filtro sea True
        if bool_filtro:
            # Selecciono la lista de genes de la columna col_genes
            curr_L_genes = curr_peak[col_genes].split(sep_genes); 
            # Recorro curr_L_genes
            for curr_gen in curr_L_genes:
                # Veo que curr_gen no este en L_genes y que curr_gen no sea un string vacio
                if str(curr_gen)!='' and (not (str(curr_gen) in L_genes)):
                    L_genes.append(str(curr_gen)); 
    # Si nom_out tiene un string, guardo la lista
    if len(nom_out):
        L_genes = guardar_lista_genes(nom_out, L_genes, sep='\n', path_out=path_out); 
    # Si nom_out es un string vacio, solo printeo L_genes
    else:
        print(L_genes)
    return L_genes


def responseBody_to_list(responseBody, L_go_id):
    # Transforma responseBody con varios elementos en lista de nombres correspondientes a los IDs de L_go_id

    # Inicializo la lista que se devuelve
    L_out = []; 
    # Inicializo un diccionario para registrar go_ids con nombres
    dict_check = {}; 
    # Recorro responseBody['results']
    for i in range(len(responseBody['results'])):
        curr_res = responseBody['results'][i]; 
        # Reviso que curr_res['id'] este en dict_check
        if not (curr_res['id'] in dict_check.keys()):
            # Si no esta, agrego la key
            dict_check[curr_res['id']] = []; 
        # Agrego curr_res['name'] a ese elemento 
        dict_check[curr_res['id']].append(curr_res['name']); 
    # Inicializo contadores para L_go_id
    cont_mult_name = 0; 
    cont_no_name = 0; 
    cont_ok = 0; 
    # Una vez recorrido responseBody['results'], recorro L_go_id para ver que todos los ids esten en dict_check
    for curr_go_id in L_go_id:
        if curr_go_id in dict_check.keys():
            if len(dict_check[curr_go_id]) > 1:
                # Mas de un nombre asociado a curr_go_id
                cont_mult_name += 1; 
            else:
                # Un nombre asociado a curr_go_id
                cont_ok += 1; 
            L_out = L_out + dict_check[curr_go_id]; 
        else:
            # curr_go_id sin nombre asociado
            cont_no_name += 1; 
    print('IDs con todo bien: ' + str(cont_ok))
    if cont_mult_name:
        print('IDs con mas de un nombre asociado: ' + str(cont_mult_name))
    if cont_no_name:
        print('IDs sin ningun nombre asociado: ' + str(cont_no_name))
    return L_out


def server_biomart(biomart_server='http://uswest.ensembl.org/biomart'):
    # Obtiene la variable server de biomart
    if biomart_server=='':
        if ib3:
            server = biomart.BiomartServer(http_proxy=http_proxy_ib3, https_proxy=https_proxy_ib3); 
        else:
            server = biomart.BiomartServer(); 
    else:
        if ib3:
            server = biomart.BiomartServer(biomart_server, http_proxy=http_proxy_ib3, https_proxy=https_proxy_ib3); 
        else:
            server = biomart.BiomartServer(biomart_server); 
    return server


def traducir_go_id(go_id):
    # Devuelve el nombre del gene ontology ID

    # Segundos a esperar si no anda requests.get()
    wait_time = 60; 
    wait_time_long = 600; 
    # Agarro solo el numero de go_id
    id_num = go_id[3:]; 
    # Defino la URL de la REST API
    requestURL = 'https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/GO%3A' + id_num; 
    # Uso GET para obtener el elemento json
    try:
        r = requests.get(requestURL, headers={ "Accept" : "application/json"}); 
    except:
        print('Conexion perdida. Esperando ' + str(wait_time) + ' segundos.')
        time.sleep(wait_time); 
        try:
            r = requests.get(requestURL, headers={ "Accept" : "application/json"}); 
        except:
            print('Segundo intento fallido. Esperando ' + str(wait_time_long) + ' segundos.')
            time.sleep(wait_time_long); 
            r = requests.get(requestURL, headers={ "Accept" : "application/json"}); 
    # Reviso que el pedido sea ok
    if not r.ok:
        print('ERROR: go_id ' + str(go_id) + ' no ok.')
        r.raise_for_status(); 
        ret = ''; 
    else:
        responseBody = r.json(); 
        if responseBody['numberOfHits'] > 1:
            print('WARNING: Mas de un hit.')
            for r in responseBody['results']:
                print(r)
        ret=responseBody['results'][0]['name']; 
    return ret


def traducir_go_id_dict(go_id):
    # Devuelve informacion del gene ontology ID

    # Segundos a esperar si no anda requests.get()
    wait_time = 60; 
    wait_time_long = 600; 
    # Agarro solo el numero de go_id
    id_num = go_id[3:]; 
    # Defino la URL de la REST API
    requestURL = 'https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/GO%3A' + id_num; 
    # Uso GET para obtener el elemento json
    try:
        r = requests.get(requestURL, headers={ "Accept" : "application/json"}); 
    except:
        print('Conexion perdida. Esperando ' + str(wait_time) + ' segundos.')
        time.sleep(wait_time); 
        try:
            r = requests.get(requestURL, headers={ "Accept" : "application/json"}); 
        except:
            print('Segundo intento fallido. Esperando ' + str(wait_time_long) + ' segundos.')
            time.sleep(wait_time_long); 
            r = requests.get(requestURL, headers={ "Accept" : "application/json"}); 
    # Defino la variable que se devuelve
    ret = []; 
    # Reviso que el pedido sea ok
    if not r.ok:
        print('ERROR: go_id ' + str(go_id) + ' no ok.')
        r.raise_for_status();  
    else:
        responseBody = r.json(); 
        # Reviso que responseBody tenga un solo hit (tiro warning sino)
        if responseBody['numberOfHits'] > 1:
            print('WARNING: Mas de un hit.')
            for r in responseBody['results']:
                ret.append(r['id']); 
                ret.append(r['isObsolete']); 
                ret.append(r['name']); 
        else:
            ret.append(responseBody['results'][0]['id']); 
            ret.append(responseBody['results'][0]['isObsolete']); 
            ret.append(responseBody['results'][0]['name']);  
    return ret


def traducir_go_id_list(L_go_id):
    # Devuelve una lista de nombres para una lista de gene ontology IDs

    # Inicializo la lista que se devuelve
    L_out = []; 
    # Defino el requestURL
    requestURL = 'https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/'; 
    # Inicializo una lista para no repetir go_ids
    L_go_unique = []; 
    for go_id in L_go_id:
        if not (go_id in L_go_unique):
            L_go_unique.append(go_id); 
            requestURL = requestURL + 'GO%3A' + str(go_id[3:]) + '%2C'; 
    # Elimino el ultimo %2C de requestURL
    requestURL = requestURL[:-3]; 
    # Busco requestURL
    r = requests.get(requestURL, headers={ "Accept" : "application/json"}); 
    # Reviso que el pedido sea ok
    if not r.ok:
        print('ERROR')
        r.raise_for_status(); 
    else:
        responseBody = r.json(); 
        # Reviso que numberOfHits sea igual al largo de L_go_id
        if responseBody['numberOfHits'] != len(L_go_id):
            print('WARNING: Numero de hits no coincide con numero de gene ontology IDs pedidos. Numero de hits: ' + 
            str(responseBody['numberOfHits']) + '; numero de pedidos: ' + str(len(L_go_id)))
        # Uso responseBody_to_list() para obtener L_out
        L_out = responseBody_to_list(responseBody, L_go_id); 
    return L_out


def traducir_go_id_list_recursivo(L_go_id, len_cutoff=20):
    # Usa recursion para implementar traducir_go_id_list() de a pedazos

    # Si L_go_id es de largo menor que len_cutoff, uso traducir_go_id_list()
    if len(L_go_id) <= len_cutoff:
        ret = traducir_go_id_list(L_go_id); 
    # Si L_go_id es de largo mayor que len_cutoff, 
    else:
        ret = traducir_go_id_list(L_go_id[:len_cutoff]) + traducir_go_id_list_recursivo(L_go_id[len_cutoff:], len_cutoff=len_cutoff); 
    return ret


def _main():
    # Funcion para probar funciones en ejecucion del archivo

    # Pruebo generar un diccionario para humano y raton
    #dict_go_human = pipeline_generar_dict_go_api('HumanGenesGOdict', 'GOid_info_human', path_in=path_in_main, path_out=path_out_main); 
    #dict_go_mouse = pipeline_generar_dict_go_api('MouseGenesGOdict', 'GOid_info_mouse', path_in=path_in_main, path_out=path_out_main); 

    # Funciones para generar listas de genes
    print('### Humano')
    print('# Total')
    L_genes_human_total = pipeline_lista_genes(nom_csv_human, col_genes=col_genes_main, L_filtros=[], path_peaks=path_rnaseq_main, 
                                               nom_out='L_genes_human_total'+nom_mod_main, path_out=path_out_main); 
    print('# Con sitio de union')
    L_genes_human_con_su = pipeline_lista_genes(nom_csv_human, col_genes=col_genes_main, L_filtros=[(col_class_main, '1-Regulado directamente')], 
                                                path_peaks=path_rnaseq_main, nom_out='L_genes_human_con_su'+nom_mod_main, path_out=path_out_main); 
    print('# Sin sitio de union')
    L_genes_human_sin_su = pipeline_lista_genes(nom_csv_human, col_genes=col_genes_main, L_filtros=[(col_class_main, '2-Regulado indirectamente')], 
                                                path_peaks=path_rnaseq_main, nom_out='L_genes_human_sin_su'+nom_mod_main, path_out=path_out_main); 
    print('### Raton')
    print('# Total')
    L_genes_mouse_total = pipeline_lista_genes(nom_csv_mouse, col_genes=col_genes_main, L_filtros=[], path_peaks=path_rnaseq_main, 
                                               nom_out='L_genes_mouse_total'+nom_mod_main, path_out=path_out_main); 
    print('# Con sitio de union')
    L_genes_mouse_con_su = pipeline_lista_genes(nom_csv_mouse, col_genes=col_genes_main, L_filtros=[(col_class_main, '1-Regulado directamente')], 
                                                path_peaks=path_rnaseq_main, nom_out='L_genes_mouse_con_su'+nom_mod_main, path_out=path_out_main); 
    print('# Sin sitio de union')
    L_genes_mouse_sin_su = pipeline_lista_genes(nom_csv_mouse, col_genes=col_genes_main, L_filtros=[(col_class_main, '2-Regulado indirectamente')], 
                                                path_peaks=path_rnaseq_main, nom_out='L_genes_mouse_sin_su'+nom_mod_main, path_out=path_out_main); 
    '''
    # Funciones para histogramas
    hist_human = pipeline_histograma('HumanGenesGOdict', 'HumanGenesGO', bars_hist=20, guardar_hist=True, path_dict_genes=path_in_main, path_hist=path_out_graficos); 
    hist_mouse = pipeline_histograma('MouseGenesGOdict', 'MouseGenesGO', bars_hist=20, guardar_hist=True, path_dict_genes=path_in_main, path_hist=path_out_graficos); 
    '''
    '''# Pruebo hacer los diccionarios para humano y raton
    print('Inicializando servidor biomart.')
    server = server_biomart(); 
    print('Server inicializado. Obteniendo diccionarios de GO terms.')

    print()
    dict_human = generar_dict_go(server, genome_name_human); 
    print('dict_human creado. Guardando.')
    dict_human = guardar_dict(dict_human, 'dict_go_id_human', path_out=path_out_main, sep2=';'); 

    print()
    dict_mouse = generar_dict_go(server, genome_name_mouse); 
    print('dict_mouse creado. Guardando.')
    dict_mouse = guardar_dict(dict_mouse, 'dict_go_id_mouse', path_out=path_out_main, sep2=';');''' 

    return ''


#################################### RESULTADOS ###################################

output_dump = []; 

if __name__=='__main__':
    output_dump.append(_main()); 

