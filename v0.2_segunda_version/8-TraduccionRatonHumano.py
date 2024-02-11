
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
Funciones para pasar identificadores de genes de raton a genes de humano

pipeline_mouse_to_human() para contar coincidencias entre genes de humano y raton

# Databases:
>Homo sapiens
hsapiens_encode, hsapiens_external_feature, hsapiens_gene_ensembl, hsapiens_genomic_sequence, hsapiens_karyotype_end, hsapiens_karyotype_start
hsapiens_marker_end, hsapiens_marker_start, hsapiens_mirna_target_feature, hsapiens_regulatory_feature, hsapiens_snp, hsapiens_snp_som
hsapiens_structvar, hsapiens_structvar_som
>Mus musculus
mmusculus_external_feature, mmusculus_gene_ensembl, mmusculus_genomic_sequence, mmusculus_karyotype_end, mmusculus_karyotype_start
mmusculus_mirna_target_feature, mmusculus_regulatory_feature, mmusculus_snp, mmusculus_structvar
'''

#################################### VARIABLES ####################################

# Direcciones de archivos y outputs
path_dropbox_casa = 'D:\\Users\\Admin\\Dropbox\\'; 
path_dropbox_ib3 = 'C:\\Users\\emili\\Dropbox\\'; 
path_fasta_casa = 'D:\\Archivos doctorado\\Genomas\\'; 
path_fasta_ib3 = 'X:\\Genomas\\'; 
path_output_dump_casa = 'D:\\Archivos doctorado\\Output_dump\\'; 
path_output_dump_ib3 = 'X:\\Output_dump\\'; 

# Ver atributos en key_attributes.txt en carpeta 3-Genes transactivados rio abajo\1-Interpretacion de datos ChIP-seq
attributes_general = ['mgi_symbol', 'ensembl_gene_id']; 
attributes_hsapiens_homolog = ['hsapiens_homolog_ensembl_gene', 'hsapiens_homolog_associated_gene_name', 'hsapiens_homolog_orthology_confidence', 'hsapiens_homolog_orthology_type']; 
attributes_extra = ['external_gene_name', 'hgnc_symbol', 'pdb']; # Hay predictores y otras cosas interesantes
attributes_mmusculus_homolog = ['mmusculus_paralog_ensembl_gene', 'mmusculus_paralog_associated_gene_name', 'mmusculus_paralog_orthology_type']; 
attributes = ['ensembl_gene_id', 'hsapiens_homolog_ensembl_gene', 'hsapiens_homolog_associated_gene_name']; 

nombre_rnaseq_raton = '6-Dupays2015_RNAseq'; 
nombre_rnaseq_humano = '9-Anderson2018_RNAseq'; 

nom_entrez_ensembl_human = '8-ensembl_to_entrez_human'; 
nom_entrez_ensembl_mouse = '8-ensembl_to_entrez_mouse'; 
nom_entrez_ensembl_mouse = 'dict_entrez_ensembl_mouse'; 


# Variables main()

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
path_pwm_main = path_dropbox_main + 'Doctorado\\3-Genes transactivados rio abajo\\1-Interpretacion de datos ChIP-seq\\PWM\\'; 
path_out_main = path_output_dump_main + 'TraduccionRatonHumano\\'; 
path_out_entrez_ensembl = path_output_dump_main + 'TraduccionEntrezEnsembl\\'; 

nom_rnaseq_main = nombre_rnaseq_humano; 


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
    with open(filepath, 'r', encoding="latin-1") as F:
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


def buscar_gene_id(gene_symbol, L_aliases, L_ensembl_id, genoma, verbose=0):
    # Busca a que gen/genes hacen referencia gene_symbol, L_aliases y L_ensembl_id
    # verbose<=0 sin prints, verbose>=1 prints de error, verbose>=2 prints de warning, verbose>=3 print del output

    # Inicializo la lista de genes que se devuelve
    L_gene_id_out = []; 
    # Inicializo contador y lista para almacenar resultados de ciclo while
    L_gene_id_symbol = []; 
    L_gene_id_aliases = []; 
    cont = -1; 
    # Armo un ciclo while para revisar gene_symbol y L_aliases
    while cont<len(L_aliases):
        # Primero pruebo con gene_symbol
        if cont == -1:
            try: 
                L_gene_id = genoma.gene_ids_of_gene_name(gene_symbol); 
                L_gene_id_symbol = L_gene_id[:]; 
            except:
                pass
        # Despues de probar con gene_symbol, pruebo con L_aliases
        else:
            try: 
                L_gene_id = genoma.gene_ids_of_gene_name(L_aliases[cont]); 
                L_gene_id_aliases = L_gene_id_aliases + L_gene_id[:]; 
            except:
                pass
        # Sumo 1 para que avance el ciclo
        cont = cont + 1; 
    # Reviso L_ensembl_id para asegurarme que sean IDs validos
    L_ensembl_id_confirmados = []; 
    for curr_ensembl_id in L_ensembl_id:
        # Pruebo gene_by_id con cada ensembl_id
        try:
            genoma.gene_by_id(curr_ensembl_id); 
            # Si no tira error, lo agrego a confirmados
            L_ensembl_id_confirmados.append(str(curr_ensembl_id)); 
        except:
            pass
    # Reviso si hay elementos en L_ensembl_id_confirmados
    if len(L_ensembl_id_confirmados) > 0:
        L_gene_id_aliases_total = L_gene_id_symbol + L_gene_id_aliases; 
        # Defino lista de ensembl_id que no se encontraron con aliases
        L_ensembl_id_ausente = []; 
        # Recorro L_ensembl_id
        for ensembl_id in L_ensembl_id_confirmados:
            # Si ensembl_id no esta entre los ids encontrados con aliases
            if not (ensembl_id in L_gene_id_aliases_total):
                L_ensembl_id_ausente.append(str(ensembl_id)); 
            # Si esta, lo agrego a L_gene_id_out
            else:
                L_gene_id_out.append(str(ensembl_id)); 
        # Defino lista de alias_id que no se encontraron con ensembl_id
        L_alias_id_ausente = []; 
        for alias_id in L_gene_id_aliases_total:
            # Si alias_id no esta en L_ensembl_id
            if not (alias_id in L_ensembl_id_confirmados):
                L_alias_id_ausente.append(str(alias_id)); 
        # Si L_gene_id_out es lista vacia, reviso L_ensembl_id_ausente y L_alias_id_ausente
        if len(L_gene_id_out) == 0:
            L_gene_id_out = L_ensembl_id_ausente + L_alias_id_ausente; 
            if verbose>=2:
                print('### WARNING en procesamiento de ' + str(gene_symbol) + ', no se confirmo ningun ensembl_id. Usando ' + str(L_gene_id_out))
    elif (len(L_gene_id_symbol) > 0) or (len(L_gene_id_aliases) > 0):
        if verbose>=2:
            print('### WARNING en procesamiento de ' + str(gene_symbol) + ', usando L_gene_id_symbol y L_gene_id_aliases. ' + str(L_gene_id_symbol) + ' ' + str(L_gene_id_aliases))
        L_gene_id_out = L_gene_id_symbol[:] + L_gene_id_aliases[:]; 
    else:
        if verbose>=1:
            print('###### ERROR en procesamiento de ' + str(gene_symbol) + ', no hay ensembl_id valido asociado.')
    if verbose>=3:
        print(L_gene_id_out)
    return L_gene_id_out


def cargar_dict(nom_dict, path_dict='.\\', ext='.csv', sep=';', sep2=','):
    # Carga los contenidos de un diccionario guardado en .csv por guardar_dict()

    # Inicializo el diccionario que se devuelve
    dict_out = {}; 
    # Defino la direccion del archivo en base a path_dict y nom_dict
    dirarch = os.path.join(path_dict, nom_dict + ext); 
    # Abro el archivo
    with open(dirarch, 'r', encoding='latin-1') as F_dict:
        print('Archivo ' + nom_dict + ext + ' abierto.')
        # Recorro F_dict con readlines()
        for curr_line in F_dict.readlines():
            curr_L = curr_line.rstrip().split(sep=sep); 
            # Reviso que curr_L solo tenga dos elementos
            if len(curr_L) != 2:
                print('ERROR en linea ' + str(curr_L) + ', no tiene dos elementos.')
            else:
                if not(curr_L[0] in dict_out.keys()):
                    dict_out[curr_L[0]] = []; 
                for id1 in curr_L[1].split(sep=sep2):
                    dict_out[curr_L[0]].append(str(id1)); 
                #dict_out[curr_L[0]] = curr_L[1].split(sep=sep2); 
    return dict_out


def cargar_dict_base(nom_dict, path_dict='.\\', ext='.csv', sep=';'):
    # Carga los contenidos de un diccionario guardado en .csv por guardar_dict()

    # Inicializo el diccionario que se devuelve
    dict_out = {}; 
    # Defino la direccion del archivo en base a path_dict y nom_dict
    dirarch = os.path.join(path_dict, nom_dict + ext); 
    # Abro el archivo
    with open(dirarch, 'r') as F_dict:
        print('Archivo ' + nom_dict + ext + ' abierto.')
        # Recorro F_dict con readlines()
        for curr_line in F_dict.readlines():
            curr_L = curr_line.rstrip().split(sep=sep); 
            # Reviso que curr_L solo tenga dos elementos
            if len(curr_L) != 2:
                print('ERROR en linea ' + str(curr_L) + ', no tiene dos elementos.')
            else:
                dict_out[curr_L[0]] = curr_L[1]; 
    return dict_out


def comparar_2dict(dict1, dict2):
    # Compara 2 diccionarios y printea cualquier cosa que no tengan en comun

    print()
    print('##### Empezando comparacion de dict1 y dict2')
    print()
    # Inicializo la variable que se devuelve
    dict_out = {'key1_missing':[], 'key1_in_dict2':[], 'key_igual':[], 'key_diferente':[], 'key2_in_dict1':[], 'key2_missing':[]}; 
    # Recorro dict1.keys()
    for key1 in dict1.keys():
        # Veo si key1 esta en dict2.keys()
        if not(key1 in dict2.keys()):
            print()
            print('# Diferencia en keys: ' + str(key1) + ' ausente en dict2.keys()')
            print()
            dict_out['key1_missing'].append(str(key1)); 
        else:
            dict_out['key1_in_dict2'].append(str(key1)); 
            if dict1[key1] == dict2[key1]: 
                dict_out['key_igual'].append(str(key1)); 
            else:
                print("> Diferencia en dict1['" + str(key1) + "']=" + str(dict1[key1]) + " con dict2['" + str(key1) + "']=" + str(dict2[key1]))
                dict_out['key_diferente'].append(str(key1)); 
    # Recorro dict2.keys()
    for key2 in dict2.keys():
        # Veo si key2 esta en dict_out['key1_in_dict2']
        if key2 in dict_out['key1_in_dict2']:
            if key2 in dict1.keys():
                # Funcionamiento normal
                dict_out['key2_in_dict1'].append(str(key2)); 
            else:
                # Error 1
                dict_out['key2_missing'].append(str(key2)); 
                print('### ERROR 1: ' + str(key2) + ' se encuentra buscando key1 en dict2 pero no se encuentra buscando key2 en dict1.')
        else:
            if key2 in dict1.keys():
                # Error 2
                dict_out['key2_in_dict1'].append(str(key2)); 
                print('### ERROR 2: ' + str(key2) + ' no se encuentra buscando key1 en dict2 pero se encuentra buscando key2 en dict1.')
                # Hay que comparar dict2[key2] con dict1[key2]
                if dict2[key2] == dict1[key2]:
                    dict_out['key_igual'].append(str(key2)); 
                else:
                    print("> Diferencia en dict2['" + str(key2) + "']=" + str(dict2[key2]) + " con dict1['" + str(key2) + "']=" + str(dict1[key2]))
                    dict_out['key_diferente'].append(str(key2)); 
            else:
                # Funcionamiento normal
                dict_out['key2_missing'].append(str(key2)); 
                print()
                print('# Diferencia en keys: ' + str(key2) + ' ausente en dict1.keys()')
                print()
    print()
    print('##### Comparacion de dict1 y dict2 finalizada')
    print()
    if 'y' in input('Desea mostrar dict_out de la comparacion? ("y" hace print)').lower():
        print(dict_out)
    return dict_out


def crear_dict_entrez_ensembl(nom_dict_entrez_ensembl, col_entrez=3, col_ensembl=0, col_nombre=1, col_alias=2, path_dict='', ext='.txt'):
    # Abre el archivo nom_dict_entrez_ensembl y lo transforma en diccionario
    # Devuelve dos diccionarios, uno con todos los valores de Entrez y sus ensembl ids asociados, otro con los nombres y alias de genes y sus ensembl ids asociados

    # Inicializo los diccionarios que se devuelven
    dict_entrez_ensembl = {}; 
    dict_nom_genes = {}; 
    # Abro el archivo nom_dict_entrez_ensembl
    M_dict_entrez_ensembl = abrir_csv(nom_dict_entrez_ensembl, dir_arch=path_dict, sep=';', ignore_headers=True, ext=ext); 
    ### Display
    len_M_dict = len(M_dict_entrez_ensembl); 
    ###
    # Recorro M_dict_entrez_ensembl
    for i in range(len_M_dict): 
        ### Display
        if i%500==0:
            print('Progreso: ' + str(i) + '/' + str(len_M_dict))
        ###
        L_curr = M_dict_entrez_ensembl[i]; 
        # Defino los valores dados en columnas entrez, ensembl, nombre y alias
        curr_entrez_id = L_curr[col_entrez]; 
        curr_ensembl_id = L_curr[col_ensembl]; 
        curr_nom_gene = L_curr[col_nombre]; 
        curr_alias_gene = L_curr[col_alias]; 
        # Cargo entrez_id y ensembl_id a dict_entrez_ensembl
        dict_entrez_ensembl = dict_append(dict_entrez_ensembl, curr_entrez_id, curr_ensembl_id); 
        # Cargo nom_gene y ensembl_id a dict_nom_genes
        dict_nom_genes = dict_append(dict_nom_genes, curr_nom_gene, curr_ensembl_id); 
        # Cargo alias y ensembl_id a dict_nom_genes
        dict_nom_genes = dict_append(dict_nom_genes, curr_alias_gene, curr_ensembl_id); 
    return dict_entrez_ensembl, dict_nom_genes


def crear_dict_orthologs(genome1, L_attributes, nombre_out, id_key=0, id_str=1, ids_extra=[], path_out='.\\'):
    # Crea un diccionario de ortologos y lo guarda

    print('>Creando dict_2genomes.')
    # Uso generar_dict_2genomes() y guardar_dict() 
    dict_2genomes = generar_dict_2genomes(genome1, L_attributes, id_key=id_key, id_str=id_str, ids_extra=ids_extra, biomart_server='http://uswest.ensembl.org/biomart'); 
    print('>dict_2genomes creado. Guardando.')
    dict_2genomes = guardar_dict(dict_2genomes, nombre_out, path_out=path_out); 
    print('>dict_2genomes guardado como "' + nombre_out + '.csv" en "' + str(path_out) + '".')
    return dict_2genomes


def dict_append(dict_in, curr_key, curr_elem):
    # Funcion que ve si dict_in tiene curr_key 
    # Si no la tiene, agrega una lista vacia a esa key y le carga curr_elem
    # Si la tiene, ve si ya esta curr_elem cargado
    # Si no lo tiene, lo agrega

    # Creo una copia de dict_in
    dict_out = copy.deepcopy(dict_in); 
    # Veo si curr_key ya esta en dict_out.keys
    if curr_key in dict_out.keys():
        # Veo si curr_elem ya esta anotado en dict_out[curr_key]
        if not (str(curr_elem) in dict_out[curr_key]):
            # Si no esta anotado, lo agrego a la lista
            dict_out[curr_key].append(str(curr_elem)); 
    elif curr_key!='':
        # Si no esta anotado, agrego curr_key y le hago append de curr_elem
        dict_out[curr_key] = []; 
        dict_out[curr_key].append(str(curr_elem)); 
    return dict_out


def extraer_genes_peaks(nom_arch, path_arch='', sep=';', ext='.csv', col_genes=0):
    # Funcion que extrae la lista de genes en columna col_genes de archivo nom_arch

    # Abro el archivo nom_arch en path_arch
    M_arch = abrir_csv(nom_arch, dir_arch=path_arch, sep=sep, ext=ext, ignore_headers=False); 
    # Inicializo la lista de genes que se devuelven
    L_genes = []; 
    # Recorro M_arch
    for i in range(len(M_arch)):
        curr_row = M_arch[i]; 
        # Reviso si el gen ya esta anotado en L_genes
        if not (curr_row[col_genes] in L_genes):
            # Si no esta anotado, lo agrego
            L_genes.append(curr_row[col_genes]); 
        else:
            # Si ya esta anotado, tiro warning (no deberia pasar)
            print('WARNING: Gen ' + str(curr_row[col_genes]) + ' repetido.')
    return L_genes


def _extraer_genes_peaks_OLD(nom_arch, path_arch='', sep=';', sep2=',', col_genes=17, col_class=4):
    # Funcion que extrae una lista de genes de un archivo de peaks clasificados

    ## Datos hardcodeados
    # Texto seleccionado en col_class
    class_select = ['1-Regulado directamente', '2-Regulado indirectamente']; 
    # Abro el archivo nom_arch en path_arch
    M_arch = abrir_csv(nom_arch, dir_arch=path_arch, sep=sep, ignore_headers=True); 
    # Inicializo la lista de genes que se devuelven
    L_genes = []; 
    # Recorro M_arch
    for i in range(len(M_arch)):
        curr_row = M_arch[i]; 
        # Reviso que la fila tenga clase correcta
        if curr_row[col_class] in class_select:
            # Extraigo los genes en col_genes
            curr_genes = curr_row[col_genes].split(sep2); 
            # Recorro curr_genes
            for j in range(len(curr_genes)):
                # Veo si curr_genes[j] esta en L_genes
                if not (str(curr_genes[j]) in L_genes):
                    # Si no esta, lo agrego a L_genes
                    L_genes.append(str(curr_genes[j])); 
    return L_genes


def generar_dict_2genomes(genome1, L_attributes, id_key=0, id_str=1, ids_extra=[], biomart_server='http://uswest.ensembl.org/biomart'):
    # Genera un diccionario que pasa de IDs de ensembl para un genoma a id de ensembl para otro genoma
    # Pensado (y testeado) para pasar de Mus musculus a Homo sapiens

    # Conexion al servidor de biomart
    server = generar_dict_2genomes_server_biomart(biomart_server); 
    # Defino los nombres de los genomas
    dataset_biomart = genome_to_dataset(genome1, modifier='_gene_ensembl'); 
    # Consigo data_used con generar_dict_2genomes_data_used()
    data_used = generar_dict_2genomes_data_used(server, dataset_biomart, L_attributes); 
    # Inicializo el diccionario que se devuelve
    genome1_to_genome2 = {}; 
    # Recorro data_general
    for curr_line in data_used.splitlines(): 
        # Uso el mismo orden que la lista de atributos (attributes_X)
        line = curr_line.split('\t'); 
        # Uso generar_dict_2genomes_line() para procesar line
        genome1_to_genome2 = generar_dict_2genomes_line(genome1_to_genome2, line, id_key, id_str, ids_extra=ids_extra); 
    return genome1_to_genome2


def generar_dict_2genomes_data_used(server, dataset, attributes, decode='ascii'):
    # Funcion que recibe la variable server de biomart y devuelve data_used

    # Agarro el dataset dentro de server
    mart = server.datasets[dataset]; 
    # Busco attributes y los guardo en response
    response = mart.search({'attributes':attributes}); 
    # Saco la info de response y la guardo en data_used
    data_used = response.raw.data.decode(decode); 
    return data_used


def generar_dict_2genomes_line(dict_in, line, id_key, id_str, ids_extra=[], sep_in='|'):
    # Funcion para procesar la linea extraida de data_usado y agregarla al diccionario
    # line es una lista de attributes
    # id_key es la posicion en line que corresponde a la key de dict_in agregada
    # id_str es la posicion en line que corresponde al valor agregado a dict_in[id_key]

    # Inicializo el diccionario que se devuelve
    dict_out = copy.deepcopy(dict_in); 
    # Reviso que los valores importantes no esten vacios y los guardo en dict
    if len(line[id_key]) and len(line[id_str]):
        # Defino line_str
        line_str = str(line[id_str]); 
        # Si hay ids en ids_extra, los agrego (a lo bruto)
        for id_extra in ids_extra:
            line_str = line_str + sep_in + str(line[id_extra]); 
        # Reviso si el valor de line[id_key] esta en keys del diccionario
        if not (str(line[id_key]) in dict_out.keys()):
            # Si no esta presente, agrego una lista vacia al diccionario primero
            dict_out[str(line[id_key])] = []; 
        # Agrego el elemento a la lista
        dict_out[str(line[id_key])].append(line_str); 
    return dict_out


def generar_dict_2genomes_server_biomart(biomart_server):
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


def generar_dict_entrez_ensembl(genome_name, biomart_server='www.ensembl.org/biomart'):
    # Genera un diccionario que pasa de IDs de entrez a IDs de ensembl
    # Ver mas info en https://autobencoder.com/2021-10-03-gene-conversion/

    # Primero defino el nombre del genoma
    if genome_name == 'mm9' or genome_name == 'mm10':
        dataset_name = 'mmusculus_gene_ensembl'; 
    elif genome_name == 'hg19' or genome_name == 'hg38':
        dataset_name = 'hsapiens_gene_ensembl'; 
    else:
        print('ERROR: Nombre de genoma ' + str(genome_name) + ' desconocido. Usando genoma de raton.')
        dataset_name = 'mmusculus_gene_ensembl'; 
    # Conexion al servidor de biomart
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
    mart = server.datasets[dataset_name]; 
    # Atributos a buscar 
    # Simbolo del gen: 'mgi_symbol'
    # Atributos ensembl: 'ensembl_transcript_id', 'ensembl_gene_id', 'ensembl_peptide_id'
    # Atributos entrez: 'entrezgene_accession', 'entrezgene_description', 'entrezgene_id', 'entrezgene_trans_name'
    attributes_general = ['mgi_symbol', 'ensembl_gene_id', 'entrezgene_accession', 'entrezgene_id']; 
    # Por si es necesario buscar mas en profundidad
    attributes_entrez = ['entrezgene_accession', 'entrezgene_description', 'entrezgene_id', 'entrezgene_trans_name']; 
    attributes_ensembl = ['ensembl_transcript_id', 'ensembl_gene_id', 'ensembl_peptide_id']; 
    # Busco attributes y los guardo en response
    response_general = mart.search({'attributes': attributes_general}); 
    # Saco la info de response y la guardo en data
    data_general = response_general.raw.data.decode('ascii'); 
    # Inicializo el diccionario que se devuelve
    entrez_to_ensembl = {}; 
    # Recorro data_general
    for line in data_general.splitlines(): 
        line = line.split('\t'); 
        # Uso el mismo orden que attributes_general
        gene_symbol = line[0]; 
        ensembl_gene_id = line[1]; 
        entrez_accession = line[2]; 
        entrez_gene_id = line[3]; 
        # Reviso que entrez_accession o entrez_gene_id no sean valores vacios y los cargo al dict
        if len(entrez_gene_id):
            entrez_to_ensembl[entrez_gene_id] = ensembl_gene_id; 
        if len(entrez_accession):
            entrez_to_ensembl[entrez_accession] = ensembl_gene_id; 
    return entrez_to_ensembl


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
    with open(dirarch, 'w', encoding='latin-1') as F_out:
        print('Archivo ' + nom_out + ext + ' creado.')
    # Lo vuelvo a abrir en modo append
    with open(dirarch, 'a', encoding='latin-1') as F_out:
        # Recorro dict_out.keys()
        for curr_key in dict_out.keys():
            # dict_out[key] contiene una lista con strings
            curr_L = dict_out[curr_key]; 
            # Transformo curr_L en un string separado por sep2
            curr_str = ''; 
            for i in curr_L:
                curr_str = curr_str + i + sep2; 
            curr_str = curr_str.rstrip(sep2); 
            # Agrego curr_key + sep + curr_str a F_out
            F_out.write(curr_key + sep + curr_str + '\n'); 
    return dict_out


def guardar_dict_simple(dict_out, nom_out, path_out='.\\', ext='.csv', sep=';'):
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
            curr_str = dict_out[curr_key]; 
            # Agrego curr_key + sep + curr_str a F_out
            F_out.write(curr_key + sep + curr_str + '\n'); 
    return dict_out


def guardar_txt(L_out, nom_out, path_out=''):
    # Guarda una lista en un archivo .txt usando str() para cada elemento

    # Defino dir_arch en base a nombre y path a usar
    if path_out == '':
        dir_arch = nom_out + '.txt'; 
    else:
        dir_arch = os.path.join(path_out, nom_out + '.txt'); 
    # Abro el archivo en dir_arch como write para borrar otro si ya existe
    with open(dir_arch, 'w') as F_out:
        print('Archivo ' + nom_out + '.txt creado.')
    # Vuelvo a abrir el archivo en modo append para guardar L_out
    with open(dir_arch, 'a') as F_out:
        # Recorro L_out
        for curr_line in L_out:
            # Guardo curr_line como str() y agrego \n al final
            F_out.write(str(curr_line) + '\n'); 
    return L_out


def pipeline_human_entrez_to_ensembl(nom_rnaseq, genoma_humano, nom_dict_entrez_ensembl='', crear_diccionario=False, path_dict_entrez_ensembl='.\\', path_rnaseq='', path_out_dict='.\\'):
    # Funcion para seleccionar los id de Entrez de resultados de RNA-seq en Anderson 2018 y pasarlos a id de Ensembl

    # Cargo la tabla de resultados de RNA-seq
    M_rnaseq, L_headers_rnaseq = abrir_csv_headers(nom_rnaseq, dir_arch=path_rnaseq); 
    # Nombres de diccionarios entrez_ensembl y nom_genes
    nom_entrez_ensembl = 'dict_entrez_ensembl_human_l1'; 
    nom_nom_genes = 'dict_nom_genes_ensembl_human_l1'; 
    # Si crear_diccionario es True, lo creo en base a nom_dict_entrez_ensembl
    if crear_diccionario:
        # Cargo el diccionario Entrez-Ensembl del archivo
        dict_entrez_ensembl, dict_nom_genes = crear_dict_entrez_ensembl(nom_dict_entrez_ensembl, path_dict=path_dict_entrez_ensembl); 
        print('Carga de diccionarios y matrices finalizada.')
        dict_entrez_ensembl = guardar_dict(dict_entrez_ensembl, nom_entrez_ensembl, path_out=path_out_dict); 
        dict_nom_genes = guardar_dict(dict_nom_genes, nom_nom_genes, path_out=path_out_dict); 
        print('Diccionarios guardados en archivos.')
    # Si es false, lo cargo de haberlo creado antes
    else:
        print('Cargando diccionarios ya generados')
        dict_entrez_ensembl = cargar_dict(nom_entrez_ensembl, path_dict=path_out_dict); 
        dict_nom_genes = cargar_dict(nom_nom_genes, path_dict=path_out_dict); 
        print('Diccionarios cargados.')
    # Inicializo la lista que se devuelve
    L_out = []; 
    # Inicializo contadores de genes bien traducidos, traducidos por nombre y no traducidos
    cont_trans_ok = 0; 
    cont_trans_name = 0; 
    cont_no_trans = 0; 
    # Recorro los genes de M_rnaseq
    for j in range(len(M_rnaseq)):
        curr_gene_L = M_rnaseq[j]; 
        # Booleano para chequear si hay que seguir buscando
        ensembl_id_agregado = False; 
        # Reviso si curr_gene_L[1] es 'NA' o un id de gen
        if curr_gene_L[1]!='NA':
            curr_entrez_id = curr_gene_L[1]; 
            # Veo si curr_entrez_id esta en dict_entrez_ensembl.keys()
            if str(curr_entrez_id) in dict_entrez_ensembl.keys():
                # Reviso los elemento de dict_entrez_ensembl[curr_entrez_id]
                for curr_ensembl_id in dict_entrez_ensembl[curr_entrez_id]:
                    # Veo que curr_ensembl_id no este en L_out
                    if not(str(curr_ensembl_id) in L_out):
                        # Veo que curr_ensembl_id se encuentre en genoma_humano
                        try:
                            # Si corre este script, curr_ensembl_id es un ID usable
                            genoma_humano.gene_by_id(curr_ensembl_id); 
                            # Agrego curr_ensembl_id a L_out y hago ensembl_id_agregado=True
                            L_out.append(str(curr_ensembl_id)); 
                            ensembl_id_agregado = True; 
                        except:
                            # Si falla genoma_humano.gene_by_id(curr_ensembl_id), asumo que curr_ensembl_id no es usable
                            print('# WARNING: No se encontro Ensembl ID ' + str(curr_ensembl_id) + ' en genoma dado.')
                    else:
                        # Si curr_ensembl_id ya esta en L_out, pongo ensembl_id_agregado como True
                        ensembl_id_agregado = True; 
            else:
                # Si curr_entrez_id no esta en dict_entrez_ensembl.keys()
                print('# WARNING: No se encontro Entrez ID ' + str(curr_entrez_id) + ' en dict_entrez_ensembl.')
        else:
            # Si curr_gene_L[1] es 'NA', tiro warning
            print('WARNING: curr_gene_L[1] es ' + str(curr_gene_L[1]))
        # Si no se encontro ningun Ensembl ID con dict_entrez_ensembl (o curr_gene_L[1] es 'NA'), pruebo dict_nom_genes
        if not ensembl_id_agregado:
            # Defino gene_name
            curr_gene_name = curr_gene_L[0]; 
            # Veo si curr_gene_name esta en dict_nom_genes.keys()
            if str(curr_gene_name) in dict_nom_genes.keys():
                # Reviso los elemento de dict_nom_genes[curr_gene_name]
                for curr_ensembl_id in dict_nom_genes[curr_gene_name]:
                    # Veo que curr_ensembl_id no este en L_out
                    if not(str(curr_ensembl_id) in L_out):
                        # Veo que curr_ensembl_id se encuentre en genoma_humano
                        try:
                            # Si corre este script, curr_ensembl_id es un ID usable
                            genoma_humano.gene_by_id(curr_ensembl_id); 
                            # Agrego curr_ensembl_id a L_out y hago ensembl_id_agregado=True
                            L_out.append(str(curr_ensembl_id)); 
                            ensembl_id_agregado = True; 
                        except:
                            # Si falla genoma_humano.gene_by_id(curr_ensembl_id), asumo que curr_ensembl_id no es usable
                            print('# WARNING: No se encontro Ensembl ID ' + str(curr_ensembl_id) + ' en genoma dado.')
                    else:
                        # Si curr_ensembl_id ya esta en L_out, pongo ensembl_id_agregado como True
                        ensembl_id_agregado = True; 
                if ensembl_id_agregado:
                    cont_trans_name += 1; 
                else:
                    cont_no_trans += 1; 
            else:
                # Si curr_gene_name no se encuentra en dict_nom_genes.keys(), ya no queda nada que hacer
                print('### WARNING: No se encontro gene name ' + str(curr_gene_name) + ' en dict_nom_genes.')
                cont_no_trans += 1; 
        else:
            # Si ya esta ensembl_id_agregado, sumo uno al contador de transcripcion ok
            cont_trans_ok += 1; 
    ### Display
    print('Cantidad de genes traducidos bien: ' + str(cont_trans_ok))
    print('Cantidad de genes traducidos por nombre: ' + str(cont_trans_name))
    print('Cantidad de genes no traducidos: ' + str(cont_no_trans))
    ###
    return L_out


def pipeline_mouse_entrez_to_ensembl(nom_rnaseq, genoma_raton, nom_dict_entrez_ensembl='', crear_diccionario=False, path_dict_entrez_ensembl='.\\', path_rnaseq='', path_out_dict='.\\'):
    # Funcion para seleccionar los id de Entrez de resultados de RNA-seq en Anderson 2018 y pasarlos a id de Ensembl

    # Cargo la tabla de resultados de RNA-seq
    M_rnaseq, L_headers_rnaseq = abrir_csv_headers(nom_rnaseq, dir_arch=path_rnaseq); 
    # Nombres de diccionarios entrez_ensembl y nom_genes
    nom_entrez_ensembl = 'dict_entrez_ensembl_mouse_l1'; 
    nom_nom_genes = 'dict_nom_genes_ensembl_mouse_l1'; 
    # Si crear_diccionario es True, lo creo en base a nom_dict_entrez_ensembl
    if crear_diccionario:
        # Cargo el diccionario Entrez-Ensembl del archivo
        dict_entrez_ensembl, dict_nom_genes = crear_dict_entrez_ensembl(nom_dict_entrez_ensembl, path_dict=path_dict_entrez_ensembl); 
        print('Carga de diccionarios y matrices finalizada.')
        dict_entrez_ensembl = guardar_dict(dict_entrez_ensembl, nom_entrez_ensembl, path_out=path_out_dict); 
        dict_nom_genes = guardar_dict(dict_nom_genes, nom_nom_genes, path_out=path_out_dict); 
        print('Diccionarios guardados en archivos.')
    # Si es false, lo cargo de haberlo creado antes
    else:
        print('Cargando diccionarios ya generados')
        dict_entrez_ensembl = cargar_dict(nom_entrez_ensembl, path_dict=path_out_dict); 
        dict_nom_genes = cargar_dict(nom_nom_genes, path_dict=path_out_dict); 
        print('Diccionarios cargados.')
    # Inicializo la lista que se devuelve
    L_out = []; 
    # Inicializo contadores de genes bien traducidos, traducidos por nombre y no traducidos
    cont_trans_ok = 0; 
    cont_trans_name = 0; 
    cont_no_trans = 0; 
    # Recorro los genes de M_rnaseq
    for j in range(len(M_rnaseq)):
        curr_gene_L = M_rnaseq[j]; 
        # Booleano para chequear si hay que seguir buscando
        ensembl_id_agregado = False; 
        # Reviso si curr_gene_L[1] es 'NA' o un id de gen
        if curr_gene_L[1]!='NA':
            curr_entrez_id = curr_gene_L[1]; 
            # Veo si curr_entrez_id esta en dict_entrez_ensembl.keys()
            if str(curr_entrez_id) in dict_entrez_ensembl.keys():
                # Reviso los elemento de dict_entrez_ensembl[curr_entrez_id]
                for curr_ensembl_id in dict_entrez_ensembl[curr_entrez_id]:
                    # Veo que curr_ensembl_id no este en L_out
                    if not(str(curr_ensembl_id) in L_out):
                        # Veo que curr_ensembl_id se encuentre en genoma_humano
                        try:
                            # Si corre este script, curr_ensembl_id es un ID usable
                            genoma_raton.gene_by_id(curr_ensembl_id); 
                            # Agrego curr_ensembl_id a L_out y hago ensembl_id_agregado=True
                            L_out.append(str(curr_ensembl_id)); 
                            ensembl_id_agregado = True; 
                        except:
                            # Si falla genoma_humano.gene_by_id(curr_ensembl_id), asumo que curr_ensembl_id no es usable
                            print('# WARNING: No se encontro Ensembl ID ' + str(curr_ensembl_id) + ' en genoma dado.')
                    else:
                        # Si curr_ensembl_id ya esta en L_out, pongo ensembl_id_agregado como True
                        ensembl_id_agregado = True; 
            else:
                # Si curr_entrez_id no esta en dict_entrez_ensembl.keys()
                print('# WARNING: No se encontro Entrez ID ' + str(curr_entrez_id) + ' en dict_entrez_ensembl.')
        else:
            # Si curr_gene_L[1] es 'NA', tiro warning
            print('WARNING: curr_gene_L[1] es ' + str(curr_gene_L[1]))
        # Si no se encontro ningun Ensembl ID con dict_entrez_ensembl (o curr_gene_L[1] es 'NA'), pruebo dict_nom_genes
        if not ensembl_id_agregado:
            # Defino gene_name
            curr_gene_name = curr_gene_L[0]; 
            # Veo si curr_gene_name esta en dict_nom_genes.keys()
            if str(curr_gene_name) in dict_nom_genes.keys():
                # Reviso los elemento de dict_nom_genes[curr_gene_name]
                for curr_ensembl_id in dict_nom_genes[curr_gene_name]:
                    # Veo que curr_ensembl_id no este en L_out
                    if not(str(curr_ensembl_id) in L_out):
                        # Veo que curr_ensembl_id se encuentre en genoma_humano
                        try:
                            # Si corre este script, curr_ensembl_id es un ID usable
                            genoma_raton.gene_by_id(curr_ensembl_id); 
                            # Agrego curr_ensembl_id a L_out y hago ensembl_id_agregado=True
                            L_out.append(str(curr_ensembl_id)); 
                            ensembl_id_agregado = True; 
                        except:
                            # Si falla genoma_humano.gene_by_id(curr_ensembl_id), asumo que curr_ensembl_id no es usable
                            print('# WARNING: No se encontro Ensembl ID ' + str(curr_ensembl_id) + ' en genoma dado.')
                    else:
                        # Si curr_ensembl_id ya esta en L_out, pongo ensembl_id_agregado como True
                        ensembl_id_agregado = True; 
                if ensembl_id_agregado:
                    cont_trans_name += 1; 
                else:
                    cont_no_trans += 1; 
            else:
                # Si curr_gene_name no se encuentra en dict_nom_genes.keys(), ya no queda nada que hacer
                print('### WARNING: No se encontro gene name ' + str(curr_gene_name) + ' en dict_nom_genes.')
                cont_no_trans += 1; 
        else:
            # Si ya esta ensembl_id_agregado, sumo uno al contador de transcripcion ok
            cont_trans_ok += 1; 
    ### Display
    print('Cantidad de genes traducidos bien: ' + str(cont_trans_ok))
    print('Cantidad de genes traducidos por nombre: ' + str(cont_trans_name))
    print('Cantidad de genes no traducidos: ' + str(cont_no_trans))
    ###
    return L_out


def pipeline_mouse_to_human(nom_peaks_humano, nom_peaks_raton, nom_dict_raton_humano, path_peaks='', path_dict='', ext_dict='.txt'):
    # Funcion que compara una lista de genes de raton con una lista de genes de humano

    # Cargo las listas de genes de nom_peaks en path_peaks
    L_genes_humano = extraer_genes_peaks(nom_peaks_humano, path_arch=path_peaks); 
    L_genes_raton = extraer_genes_peaks(nom_peaks_raton, path_arch=path_peaks); 
    # Inicializo los valores que se cuentan
    cont_emparejados = 0; 
    L_emparejados = []; 
    # Cargo el diccionario humano-raton desde una matriz
    ### Creado con 8-dict_biomart_mmusculus_to_hsapiens.txt en mente
    M_dict_humano_raton = abrir_csv(nom_dict_raton_humano, dir_arch=path_dict, ext=ext_dict, sep='\t', ignore_headers=True); 
    # Recorro M_dict_humano_raton
    for i in range(len(M_dict_humano_raton)):
        curr_trad = M_dict_humano_raton[i]; 
        # Raton en posicion 0, humano en posicion 3
        if (curr_trad[0] in L_genes_raton) and (curr_trad[3] in L_genes_humano):
            L_emparejados.append([str(curr_trad[0]), str(curr_trad[3])])
            cont_emparejados += 1; 
    ### Display
    print('>Genes en lista de genes de humano: ' + str(len(L_genes_humano)))
    print('>Genes en lista de genes de raton: ' + str(len(L_genes_raton)))
    print('>Genes en comun entre las listas de humano y raton: ' + str(cont_emparejados))
    ###
    return L_emparejados, cont_emparejados


def seleccionar_genes_rnaseq(nom_rnaseq, genoma, dict_entrez_ensembl, path_rnaseq=''):
    # Selecciona la lista de genes que aparezcan en un archivo de resultados de RNA-seq
    # El archivo de resultados tiene que tener headers y ser un .csv

    # Abro el archivo RNA-seq
    M_rnaseq, L_headers_rnaseq = abrir_csv_headers(nom_rnaseq, dir_arch=path_rnaseq); 
    # Agarro los ids de las columnas que me interesan
    gene_id_col = L_headers_rnaseq.index('Gene ID'); 
    gene_symbol_col = L_headers_rnaseq.index('Gene Symbol'); 
    aliases_col = L_headers_rnaseq.index('Aliases'); 
    entrez_id_col = L_headers_rnaseq.index('Entrez ID'); 
    ensembl_id_col = L_headers_rnaseq.index('Ensembl ID'); 
    # Inicializo la lista de genes que se devuelven
    L_genes = []; 
    # Recorro M_rnaseq
    for i in range(len(M_rnaseq)):
        curr_L = M_rnaseq[i]; 
        # Contadores para estatus de busqueda de id
        curr_buscando = True; 
        curr_status = 0; 
        # Extraigo la info de las columnas que me interesan
        curr_gene_id = curr_L[gene_id_col]; 
        curr_gene_symbol = curr_L[gene_symbol_col]; 
        curr_aliases = curr_L[aliases_col].split('|'); 
        curr_entrez_id = curr_L[entrez_id_col]; 
        curr_ensembl_id = curr_L[ensembl_id_col].split('|'); 
        # Busco curr_entrez_id en dict_entrez_ensembl
        if curr_entrez_id in dict_entrez_ensembl.keys():
            id_ensembl = dict_entrez_ensembl[curr_entrez_id]; 
            # Pruebo si el ensembl_id encontrado esta en el genoma
            try:
                genoma.gene_by_id(id_ensembl); 
                L_genes.append(str(id_ensembl)); 
                curr_buscando = False; 
            except:
                curr_status = 1; 
        # Si no se encontro id_ensembl con dict_entrez_ensembl
        if curr_buscando:
            # Uso una funcion para buscar el id del gen que referencian gene_symbol, aliases y/o ensembl_id
            L_gene_symbol = buscar_gene_id(curr_gene_symbol, curr_aliases, curr_ensembl_id, genoma); 
            # Si L_gene_symbol tiene elementos, los sumo a L_genes
            if len(L_gene_symbol) > 0:
                # Si hay mas de un id, tiro warning
                if len(L_gene_symbol) > 1:
                    print()
                    print('WARNING: mas de un gene symbol asociado a ' + str(curr_gene_symbol))
                    print()
                L_genes = L_genes + L_gene_symbol[:]; 
                curr_buscando = False; 
        # Hago una ultima revision 
        if curr_buscando:
            # Busco primero por gene_symbol y despues por aliases
            try:
                L_ids = genoma.gene_ids_of_gene_name(curr_gene_symbol); 
                print('WARNING. Otros IDs encontrados: ' + str(L_ids))
            except:
                pass
    return L_genes


def traducir_raton_humano_ensembl(id_gen, dict_raton_humano, genoma=''):
    # Pasa id_gen en raton a id_gen en humano, usando una serie de diccionarios y haciendo algunos chequeos
    # Genoma tiene que ser el genoma al cual se traducen los ids (probado para hg19)

    # Inicializo la lista de ids que se devuelve y el codigo de error en 0
    L_ids_out = []; 
    cod_error = 0; 
    # Primero veo que id_gen este en dict_raton_humano.keys()
    if id_gen in dict_raton_humano.keys():
        # El elemento dict_raton_humano[id_gen] deberia ser una lista
        L_id_gen_dict = dict_raton_humano[id_gen]; 
        # Recorro cada elemento de L_id_gen_dict
        for id_gen_dict in L_id_gen_dict:
            # Veo si genoma tiene un elemento
            if str(genoma) != '':
                # Si tiene un elemento, uso try: genoma.gene_by_id()
                try:
                    genoma.gene_by_id(id_gen_dict); 
                    # Si funciona, agrego id_gen_dict a L_ids out
                    L_ids_out.append(str(id_gen_dict)); 
                except:
                    print('WARNING: id ' + str(id_gen_dict) + ' no encontrado en genoma.')
                    cod_error = 2; 
            # Si no hay genoma para chequear, devuelvo todos los elementos
            else:
                L_ids_out.append(str(id_gen_dict)); 
    # Si id_gen no esta en dict_raton_humano.keys() tiro error
    else:
        print('ERROR: id ' + str(id_gen) + ' no encontrado en diccionario.')
        cod_error = 1; 
    # cod_error=1 significa que el id no se encontro en dict_raton_humano
    # cod_error=2 significa que el id no se encontro en genoma
    return L_ids_out, cod_error


def _main():
    # Funcion para probar funciones en ejecucion del archivo

    '''# Ver atributos en key_attributes.txt en carpeta 3-Genes transactivados rio abajo\1-Interpretacion de datos ChIP-seq
    #attributes_general = ['mgi_symbol', 'ensembl_gene_id']; 
    #attributes_hsapiens_homolog = ['hsapiens_homolog_ensembl_gene', 'hsapiens_homolog_associated_gene_name', 'hsapiens_homolog_orthology_confidence', 'hsapiens_homolog_orthology_type']; 
    #attributes_extra = ['external_gene_name', 'hgnc_symbol', 'pdb']; # Hay predictores y otras cosas interesantes
    #attributes_mmusculus_homolog = ['mmusculus_paralog_ensembl_gene', 'mmusculus_paralog_associated_gene_name', 'mmusculus_paralog_orthology_type']; '''
    ### Lineas que generan el diccionario mmusculus_to_hsapiens
    #attributes_mmusculus_hsapiens = ['ensembl_gene_id', 'hsapiens_homolog_ensembl_gene']; 
    #mmusculus_to_hsapiens = crear_dict_orthologs('mm9', attributes_mmusculus_hsapiens, 'mmusculus_to_hsapiens_base', path_out=path_out_main); 
    #mmusculus_to_hsapiens_load = cargar_dict('mmusculus_to_hsapiens_base', path_dict=path_out_main); 
    #dict_out = comparar_2dict(mmusculus_to_hsapiens, mmusculus_to_hsapiens_load); 

    ### Lineas que usan mmusculus_to_hsapiens
    #mmusculus_to_hsapiens = cargar_dict('mmusculus_to_hsapiens_base', path_dict=path_out_main); 
    #L_rnaseq_human = _test_traducir_genes_rnaseq(nom_rnaseq_main, mm9, hg19, nom_dict_raton='mmusculus_to_hsapiens_base', nom_dict_entrez='dict_entrez_ensembl', 
    #                                             path_dict_raton_humano=path_out_main, path_rnaseq=path_rnaseq_main); 
    #L_rnaseq_human = guardar_txt(L_rnaseq_human, 'genes_rnaseq_dupays_traduccion_hg19', path_out=path_out_main); 

    ### Lineas para cargar datos de RNA-seq de Anderson 2018
    #L_rnaseq_human = pipeline_human_entrez_to_ensembl(nombre_rnaseq_humano, hg19, nom_dict_entrez_ensembl=nom_entrez_ensembl_human, crear_diccionario=False, 
    #                                                  path_dict_entrez_ensembl=path_rnaseq_main, path_rnaseq=path_rnaseq_main, path_out_dict=path_out_entrez_ensembl); 
    #L_rnaseq_human = guardar_txt(L_rnaseq_human, nom_out='9-lista_genes_rnaseq_anderson', path_out=path_rnaseq_main); 

    ### Lineas para crear diccionarios entrez to ensembl
    #dict_entrez_ensembl_mouse = generar_dict_entrez_ensembl('mm9', biomart_server='www.ensembl.org/biomart'); 
    #dict_entrez_ensembl_mouse = guardar_dict_simple(dict_entrez_ensembl_mouse, nom_entrez_ensembl_mouse, path_out=path_out_entrez_ensembl, ext='.csv'); 
    #dict_entrez_ensembl_mouse = crear_dict_entrez_ensembl(nom_entrez_ensembl_mouse, col_entrez=3, col_ensembl=0, col_nombre=1, col_alias=2, path_dict=path_rnaseq_main); 
    ### Lineas para cargar datos de RNA-seq de Dupays 2018
    #L_rnaseq_mouse = pipeline_human_entrez_to_ensembl(nombre_rnaseq_raton, mm9, nom_dict_entrez_ensembl=nom_entrez_ensembl_mouse, crear_diccionario=False, 
    #                                                  path_dict_entrez_ensembl=path_rnaseq_main, path_rnaseq=path_rnaseq_main, path_out_dict=path_out_entrez_ensembl); 
    #L_rnaseq_mouse = guardar_txt(L_rnaseq_mouse, nom_out='9-lista_genes_rnaseq_anderson', path_out=path_rnaseq_main); 

    # Cuento coincidencias entre humano y raton con pipeline_mouse_to_human
    nom_human = '15-anderson_100kpb_updown_genes_procesado'; 
    nom_mouse = '15-dupays_100kpb_updown_genes_procesado'; 
    nom_dict = '8-dict_biomart_mmusculus_to_hsapiens'; 
    L_out, cont_out = pipeline_mouse_to_human(nom_human, nom_mouse, nom_dict, path_peaks=path_rnaseq_main, path_dict=path_rnaseq_main, ext_dict='.txt'); 
    return ''


def _test_traducir_genes_rnaseq(nom_rnaseq, genoma_raton, genoma_humano, nom_dict_raton='mmusculus_to_hsapiens_base', nom_dict_entrez='dict_entrez_ensembl', 
                                path_dict_raton_humano='.\\', path_rnaseq=''):
    # Corre seleccionar_genes_rnaseq() para obtener una lista de genes y los busca con traducir_raton_humano_ensembl()

    # Inicializo lista que se devuelve
    L_out = []; 
    # Cargo el diccionario entrez-ensembl
    dict_entrez_ensembl = cargar_dict_base(nom_dict_entrez, path_dict=path_rnaseq_main); 
    # Genero la lista de genes en nom_rnaseq
    L_genes_rnaseq = seleccionar_genes_rnaseq(nom_rnaseq, genoma_raton, dict_entrez_ensembl, path_rnaseq=path_rnaseq); 
    # Cargo el diccionario raton-humano
    dict_raton_humano = cargar_dict(nom_dict_raton, path_dict=path_dict_raton_humano); 
    # Inicializo contadores
    cont_genes_encontrados = 0; 
    cont_genes_no_encontrados = 0; 
    cont_genes_no_encontrados_en_dict = 0; 
    cont_genes_no_encontrados_en_genoma = 0; 
    cont_error_grave = 0; 
    # Recorro cada uno de los genes
    for gene_rnaseq in L_genes_rnaseq:
        # Paso el gen por traducir_raton_humano_ensembl()
        L_traducido, cod_error = traducir_raton_humano_ensembl(gene_rnaseq, dict_raton_humano, genoma=genoma_humano); 
        # cod_error=1 significa que el id no se encontro en dict_raton_humano
        # cod_error=2 significa que el id no se encontro en genoma
        # Veo si L_traducido contiene algun elemento
        if len(L_traducido):
            cont_genes_encontrados += 1; 
            # Agrego L_traducido a L_out
            L_out = L_out + L_traducido[:]; 
        # Si no tiene ningun elemento
        else:
            cont_genes_no_encontrados += 1; 
            if cod_error==1:
                cont_genes_no_encontrados_en_dict += 1; 
            elif cod_error==2:
                cont_genes_no_encontrados_en_genoma += 1; 
            else:
                print('ERROR: cod_error=' + str(cod_error) + ' con len(L_traducido)==0')
                cont_error_grave += 1; 
    ### Display final
    print('>>>>>>>>>>>> RESULTADOS FINALES')
    print('Total de genes RNA-seq: ' + str(len(L_genes_rnaseq)))
    print('Genes de RNA-seq encontrados en dict_raton_humano: ' + str(cont_genes_encontrados))
    print('Genes de RNA-seq no traducidos: ' + str(cont_genes_no_encontrados))
    print('Genes de RNA-seq no encontrados en dict_raton_humano: ' + str(cont_genes_no_encontrados_en_dict))
    print('Genes de RNA-seq no encontrados en el genoma dado: ' + str(cont_genes_no_encontrados_en_genoma))
    print('Errores graves: ' + str(cont_error_grave))
    ###
    return L_out


#################################### RESULTADOS ###################################

output_dump = []; 

if __name__=='__main__':
    output_dump.append(_main()); 

