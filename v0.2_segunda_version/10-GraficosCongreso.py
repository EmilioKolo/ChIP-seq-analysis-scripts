
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
Congreso A2B2C 2022
Funciones para generar graficos y tablas con numeros que voy a presentar

_print_valores_venn() para mostrar conteos de distintas combinaciones de criterios de busqueda sobre la tabla

subset_M() permite seleccionar algunos elementos de una matriz por criterios medianamente complejos

contar_genes_cercanos_matriz() cuenta la cantidad de genes en una lsita de peaks cortada usando subset_M()

_venn_v2() genera conteos usados para los diagramas de Venn y de torta a presentar el 09-11

dist_2_sitios() cuenta la distancia entre dos sitios de union en formato de las tablas PeaksClasificados: [ini]_[end]

generar_tabla_procesada() agrega las distancias entre sitios de union y crea tablas nuevas con informacion mas importante mejor ordenada
    * Sigue teniendo problemas con el display de las distancias

extraer_genes_rnaseq() extrae la lista de genes (no repetidos) en una columna dada (pensado para genes rnaseq)

pipeline_go_terms() genera un diccionario de Ensembl IDs a GO terms y hace un conteo para histograma

### FALTA:
- Revisar GO terms de listas de genes
	- Pasar GO terms de ID a valor
	- Ver sitios con mas de un gen RNA-seq cerca
		* Si son parte de cluster, misma via, etc.
###
'''


#################################### VARIABLES ####################################

# Direcciones de archivos y outputs para ib3 y casa
path_dropbox_casa = 'D:\\Users\\Admin\\Dropbox\\'; 
path_dropbox_ib3 = 'C:\\Users\\emili\\Dropbox\\'; 
path_fasta_casa = 'D:\\Archivos doctorado\\Genomas\\'; 
path_fasta_ib3 = 'X:\\Genomas\\'; 
path_output_dump_casa = 'D:\\Archivos doctorado\\Output_dump\\'; 
path_output_dump_ib3 = 'X:\\Output_dump\\'; 

# Nombres de archivos
L_nombres_in = ['Anderson1M', 'Anderson100k', 'Anderson50k', 'Dupays1M', 'Dupays100k', 'Dupays50k']; 
# Nombres usados (principalmente 0 y 1)
L_nombres_main = ['Anderson100k', 'Dupays100k', 'Anderson1M', 'Dupays1M']; 
L_especies_main = ['humano', 'raton', 'humano', 'raton']; 


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
path_in_main = path_dropbox_main + 'Doctorado\\2022 Congreso A2B2C\\Tablas crudas\\'; 
path_out_tablas_main = path_output_dump_main + 'ArchivosCongreso\\Tablas procesadas\\'; 
path_out_main = path_output_dump_main + 'ArchivosCongreso\\Graficos\\'; 
path_out_go = path_output_dump_main + 'GOterms\\'; 
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


def abrir_csv_headers(nom_arch, path_arch='', ext='.csv', sep=';'):
    # Usa abrir_csv() y abrir_headers() para devolver una matriz extraida de un .csv y sus headers

    # Genero matriz csv sin headers
    M_csv = abrir_csv(nom_arch, dir_arch=path_arch, ext=ext, sep=sep, ignore_headers=True); 
    # Genero lista headers
    L_headers = abrir_headers(nom_arch, dir_arch=path_arch, ext=ext, sep=sep); 
    return M_csv, L_headers


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


def buscar_elemento_en_matriz(M_in, busq, col):
    # Busca un elemento busq en columna col de M_in
    # Devuelve cantidad de instancias de busq en M_in

    # Inicializo el valor que se devuelve
    cont_busq = 0; 
    # Recorro M_in
    for i in range(len(M_in)):
        curr_L = M_in[i]; 
        # Veo si curr_L[col] es igual a busq
        if curr_L[col] == busq:
            cont_busq += 1; 
    return cont_busq


def buscar_L_elementos_en_M(M_in, L_busq, L_col):
    # Busca una lista de elementos L_busq en una lista de columnas correspondientes L_col de M_in
    # Devuelve la cantidad de instancias donde todo corresponde

    # Inicializo el valor que se devuelve
    cont_busq = 0; 
    # Reviso que L_col y L_busq tengan el mismo largo
    if len(L_busq) != len(L_col):
        print('ERROR: L_busq y L_col no tienen el mismo largo. Se devuelve 0.')
    # Si tienen el mismo largo, sigo
    else:
        # Recorro M_in
        for i in range(len(M_in)):
            curr_L = M_in[i]; 
            # Inicializo el booleano que confirma que todo esta bien
            bool_match = True; 
            # Veo para cada elemento en L_busq, que se encuentre en su elemento correspondiente de L_col
            for b in range(len(L_busq)):
                busq = L_busq[b]; 
                col = L_col[b]; 
                # Veo que curr_L[col] sea igual a busq
                if curr_L[col] != busq:
                    # Si alguno de los elementos es distinto, bool_match da False
                    bool_match = False; 
            # Despues de pasar por todos los elementos de L_busq y L_col, veo si bool_match sigue siendo True
            if bool_match:
                cont_busq += 1; 
    return cont_busq


def contar_genes_cercanos_matriz(M_in, L_head, titulo_col='genes_rnaseq', sub_sep=',', verbose=False, devolver_L=False):
    # Funcion que cuenta la cantidad de genes unicos en columna titulo_col de L_head dentro de matriz M_in

    # Inicializo la lista de genes que se cuenta o se devuelve
    L_genes = []; 
    # Recorro M_in
    for i in range(len(M_in)):
        curr_peak = M_in[i]; 
        # Busco genes en el peak con devolver_genes_cercanos_peak()
        curr_L_genes = devolver_genes_cercanos_peak(curr_peak, L_head, titulo_col=titulo_col, sub_sep=sub_sep, verbose=verbose); 
        # Recorro curr_L_genes
        for curr_gen in curr_L_genes:
            # Veo si cada gen esta en L_genes
            if not (str(curr_gen) in L_genes):
                # Si no esta, lo agrego
                L_genes.append(str(curr_gen)); 
    # Defino el valor que se devuelve
    if devolver_L:
        ret = L_genes; 
    else:
        ret = len(L_genes)
    return ret


def devolver_genes_cercanos_peak(L_peak, L_head, titulo_col='genes_rnaseq', sub_sep=',', verbose=False):
    # Funcion que recibe un peak y devuelve la lista de genes (sin repetidos) en la columna titulo_col

    # Inicializo la lista que se devuelve
    L_out = []; 
    # Defino columna y string con genes
    index_col = L_head.index(titulo_col); 
    str_genes = L_peak[index_col]; 
    # Defino la lista de genes en base a str_genes
    L_genes = str_genes.split(sep=sub_sep); 
    # Recorro L_genes para ver que ninguno este repetido
    for gen in L_genes:
        # Reviso que el gen no sea string vacio
        if str(gen) != '':
            if not(str(gen) in L_out):
                L_out.append(str(gen)); 
            elif verbose:
                print('Gen ' + str(gen) + ' repetido en peak ' + str(L_peak[:4]))
    return L_out


def dist_2_sitios(sitio1, sitio2, sep='_'):
    # Funcion que recibe dos sitios en formato [ini]sep[end] y devuelve su distancia
    # Si se solapan, devuelve 0

    # Defino posiciones iniciales y finales de sitio1 y sitio2
    L_sitio1 = sitio1.split(sep=sep); 
    L_sitio1[0] = int(L_sitio1[0]); 
    L_sitio1[1] = int(L_sitio1[1]); 
    L_sitio2 = sitio2.split(sep=sep); 
    L_sitio2[0] = int(L_sitio2[0]); 
    L_sitio2[1] = int(L_sitio2[1]); 
    # Defino sitio min y sitio max
    L_sitio_min = min(L_sitio1, L_sitio2); 
    L_sitio_max = max(L_sitio1, L_sitio2); 
    # Si se solapan, sitio_max[0] es menor a sitio_min[1]
    if L_sitio_max[0] <= L_sitio_min[1]:
        dist = 0; 
    # Si no se solapan, la distancia es L_sitio_max[0] - L_sitio_min[1]
    else: 
        dist = L_sitio_max[0] - L_sitio_min[1]; 
    return dist


def extraer_genes_rnaseq(M_csv, col_rnaseq, sub_sep=','):
    # Extrae todos los genes en la columna col_rnaseq de M_csv

    # Inicializo la lista que se devuelve
    L_genes = []; 
    # Recorro M_csv
    for i in range(len(M_csv)):
        curr_row = M_csv[i]; 
        # Defino el string donde se encuentra la lista de genes RNA-seq
        str_rnaseq = curr_row[col_rnaseq]; 
        # Paso el string a lista
        L_rnaseq = str_rnaseq.split(sep=sub_sep); 
        # Reviso que L_rnaseq tenga elementos y no sea una lista con un elemento vacio
        if len(L_rnaseq) and L_rnaseq != ['']:
            # Recorro L_rnaseq
            for gene_id in L_rnaseq:
                # Veo que gene_id no este en L_genes
                if not (str(gene_id) in L_genes):
                    L_genes.append(str(gene_id)); 
    return L_genes


def generar_data_used(server, dataset, attributes, decode='ascii'):
    # Funcion que recibe la variable server de biomart y devuelve data_used

    # Agarro el dataset dentro de server
    mart = server.datasets[dataset]; 
    # Busco attributes y los guardo en response
    response = mart.search({'attributes':attributes}); 
    # Saco la info de response y la guardo en data_used
    data_used = response.raw.data.decode(decode); 
    return data_used


def generar_tabla_procesada(nom_in, nom_out, nom_nkx25, L_nom_tf, path_in='', path_out='.\\', val_true='X', sub_sep=',', num_sep='+', sep_sitio='_'):
    # Agrega informacion relevante en orden relevante a tablas crudas

    # Abro el archivo nom_in y guardo todo en M_in y L_head
    M_in, L_head = abrir_csv_headers(nom_in, path_arch=path_in); 
    # Defino la matriz M_out con sus titulos L_head_out que se van a guardar en archivo nom_out y devolver en return
    M_out = []; 
    L_head_out = L_head[:]; 
    # Recorro M_in
    for i in range(len(M_in)):
        curr_peak = M_in[i]; 
        # Defino la lista de valores que se devuelven, inicializada con curr_peak
        L_out = curr_peak[:]; 

        ## Datos de sitios de union de NKX2-5
        # Defino la columna de sitio de union para NKX2-5
        col_su_nkx25 = L_head.index(nom_nkx25 + '_su_confirmado'); 
        # Reviso si hay sitio de union en el peak o no
        hay_nkx25 = curr_peak[col_su_nkx25]==val_true; 
        # Posicion de hay_nkx25
        pos_hay_nkx25 = 5; 
        # Reorganizo la confirmacion de sitio de union de NKX2-5 para ir al principio
        if hay_nkx25:
            L_out = L_out[:pos_hay_nkx25] + [str(val_true)] + L_out[pos_hay_nkx25:]; 
        else:
            L_out = L_out[:pos_hay_nkx25] + [''] + L_out[pos_hay_nkx25:]; 
        # Solo hago esto la primera vuelta
        if i == 0: 
            L_head_out = L_head_out[:pos_hay_nkx25] + [nom_nkx25 + '_su_conf'] + L_head_out[pos_hay_nkx25:]; 

        ## Datos de sitios de union de otros factores de transcripcion
        # Defino la lista de columnas para factores de transcripcion en L_nom_tf
        L_col_tf = nom_a_col_tf(L_head, L_nom_tf); 
        # Agrego valor para registrar si hay otros fatores de transcripcion
        hay_tf = ver_tf_peak(curr_peak, L_col_tf, val_true=val_true); 
        # Posicion de hay_tf
        pos_hay_tf = pos_hay_nkx25+1; 
        # Agrego 'X' si hay tf o '' si no hay tf
        if hay_tf:
            L_out = L_out[:pos_hay_tf] + [str(val_true)] + L_out[pos_hay_tf:]; 
        else:
            L_out = L_out[:pos_hay_tf] + [''] + L_out[pos_hay_tf:]; 
        # Solo hago esto la primera vuelta
        if i == 0: 
            L_head_out = L_head_out[:pos_hay_tf] + ['otro_su_pssm'] + L_head_out[pos_hay_tf:]; 

        ## Distancias entre sitios de union de NKX2-5 y de otros factores de transcripcion
        # Inicializo la lista de valores que se agregan a L_out
        L_dist_su = []; 
        L_head_dist = []; 
        # Agrego un elemento por factor de transcripcion en L_nom_tf a las listas de valores
        for t in range(len(L_nom_tf)):
            # Agrego una lista a L_dist_su 
            L_dist_su.append(['','']); 
            # Tambien aprovecho a agregar elementos a L_head_dist
            L_head_dist.append(L_nom_tf[t]+'_dist'); 
            L_head_dist.append(L_nom_tf[t]+'_min_dist'); 
        # Si hay sitios de union para NKX2-5 y otro factor de transcripcion, calculo distancias
        if hay_nkx25 and hay_tf:
            # Primero defino las posiciones de los sitios de union de NKX2-5
            col_su_lista = L_head.index('pos_su_lista'); 
            col_su_pssm_nkx25 = L_head.index('pos_pssm_' + nom_nkx25); 
            # Strings correspondientes a esas posiciones
            str_su_lista = curr_peak[col_su_lista]; 
            str_su_pssm_nkx25 = curr_peak[col_su_pssm_nkx25]; 
            # Listas correspondientes a las posiciones
            if str_su_lista != '':
                L_su_lista = str_su_lista.split(sep=sub_sep); 
            else:
                L_su_lista = []; 
            if str_su_pssm_nkx25 != '':
                L_su_pssm_nkx25 = str_su_pssm_nkx25.split(sep=sub_sep); 
            else:
                L_su_pssm_nkx25 = []; 
            # Junto ambas listas para simplificar el analisis
            L_su_total = L_su_lista + L_su_pssm_nkx25; 
            # Recorro cada sitio de union
            for curr_su_nkx25 in L_su_total:
                # Recorro cada uno de los factores de transcripcion
                for t in range(len(L_nom_tf)):
                    curr_tf = L_nom_tf[t]; 
                    # Defino la columna de posiciones de sitios de union para ese tf
                    col_pos_tf = L_head.index('pos_pssm_'+curr_tf); 
                    # String correspondiente a las posiciones de los sitios de union para el tf
                    str_pos_su_tf = curr_peak[col_pos_tf]; 
                    # Lista correspondiente a las posiciones de los sitios de union para el tf
                    if str_pos_su_tf == '':
                        L_pos_su_tf = []; 
                    else:
                        L_pos_su_tf = str_pos_su_tf.split(sep=sub_sep); 
                    # Inicializo los valores que se agregan a L_out
                    str_dist_su = ''; 
                    min_dist = ''; 
                    # Recorro cada uno de los sitios de union en L_pos_su_tf
                    for pos_su_tf in L_pos_su_tf:
                        # Calculo la distancia entre curr_su_nkx25 y pos_su_tf
                        dist = dist_2_sitios(curr_su_nkx25, pos_su_tf, sep=sep_sitio); 
                        # Agrego la distancia a str_dist_su
                        str_dist_su = str_dist_su + str(dist) + num_sep; 
                        # Si dist es menor que min_dist, actualizo min_dist
                        if min_dist=='' or dist < min_dist:
                            min_dist = int(dist); 
                    # Agrego str_dist_su y min_dist a L_dist_su
                    L_dist_su[t][0] = L_dist_su[t][0] + str_dist_su; 
                    if L_dist_su[t][1] == '':
                        L_dist_su[t][1] = str(min_dist); 
                    else:
                        L_dist_su[t][1] = str(min(int(L_dist_su[t][1]), min_dist)); 
        # Agrego las columnas a L_out y L_head
        for u in range(len(L_nom_tf)):
            L_out.append(L_dist_su[u][0].rstrip(num_sep)); 
            L_out.append(L_dist_su[u][1].rstrip(num_sep)); 
        # Solo hago esto la primera vuelta
        if i == 0: 
            L_head_out = L_head_out + L_head_dist; 
        # Agrego L_out a M_out
        M_out.append(L_out[:]); 
    # Guardo M_out y L_head out usando guardar_tabla_final()
    M_out = guardar_tabla_final(nom_out, M_out, path_out=path_out, L_head=L_head_out); 
    return M_out


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
                curr_str = curr_str + i + sep2; 
            curr_str = curr_str.rstrip(sep2); 
            # Agrego curr_key + sep + curr_str a F_out
            F_out.write(curr_key + sep + curr_str + '\n'); 
    return dict_out


def guardar_tabla_final(nom_out, M_out, path_out='.\\', ext='.csv', sep=';', L_head=[]):
    # Guarda los contenidos de la matriz M_out en el archivo nom_out ubicado en path_out

    # Defino la direccion del archivo en base a path_out y nom_out
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
                str_head = str_head + str(h) + sep; 
            # Elimino la ultima ocurrencia de sep
            str_head = str_head.rstrip(sep); 
            # Agrego end of line y guardo en F_out
            F_out.write(str_head + '\n'); 
        # Recorro M_out
        for L_out in M_out:
            curr_str = ''; 
            # Recorro L_out
            for i in L_out:
                curr_str = curr_str + str(i) + sep; 
            # Elimino la ultima ocurrencia de sep
            curr_str = curr_str.rstrip(sep); 
            # Agrego end of line y guardo en F_out
            F_out.write(curr_str + '\n'); 
    return M_out


def nom_a_col_tf(L_head, L_nom_tf, mod_despues='_pssm'):
    # Recibe una lista de nombres de tf y una lista de titulos en headers
    # Busca la posicion en L_head de cada uno de los tf (agrega modificador)

    # Inicializo la lista que se devuelve
    L_out = []; 
    # Recorro L_nom_tf
    for curr_nom_tf in L_nom_tf:
        curr_tf_col = L_head.index(curr_nom_tf + mod_despues)
        L_out.append(int(curr_tf_col))
    return L_out


def pipeline_go_terms(nom_csv, attributes, genome_name, col_rnaseq='genes_rnaseq', dir_csv=''):
    # 

    # Abro el archivo nom_csv
    M_csv, L_headers = abrir_csv_headers(nom_csv, path_arch=dir_csv); 
    # Extraigo todos los genes RNA-seq
    L_genes_rnaseq = extraer_genes_rnaseq(M_csv, L_headers.index(col_rnaseq)); 
    # Inicializo server de biomart
    print('Inicializando servidor biomart.')
    server = server_biomart(); 
    dataset = genome_to_dataset(genome_name); 
    data_used = generar_data_used(server, dataset, attributes, decode='ascii'); 
    print('Server inicializado. Iniciando recorrido de data_used.')
    # Inicializo un diccionario de GO terms
    dict_go_terms = {}; 
    dict_go_gene = {}; 
    # Recorro data_used
    for curr_line in data_used.splitlines(): 
        # Uso el mismo orden que la lista de atributos (attributes_X)
        line = curr_line.split('\t'); 
        print(line)
        # Veo que ensembl_gene_id este en line y que line[1] sea un GO term
        if (line[0] in L_genes_rnaseq) and line[1]!='':
            # Veo que str(line[1]) no este en dict_go_term
            if not (str(line[1]) in dict_go_terms.keys()):
                dict_go_terms[str(line[1])] = 1; 
            else:
                dict_go_terms[str(line[1])] += 1; 
            # Veo si str(line[0]) no esta en dict_go_gene
            if not (str(line[0]) in dict_go_gene.keys()):
                # Agrego str(line[0]) como key con lista vacia
                dict_go_gene[str(line[0])] = []; 
            # Reviso que str(line[1]) no este en dict_go_gene[str(line[0])]
            if not (str(line[1]) in dict_go_gene[str(line[0])]):
                # Agrego el go term a dict_go_gene[str(line[0])]
                dict_go_gene[str(line[0])].append(str(line[1])); 
    return dict_go_terms, dict_go_gene


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


def subset_M(M_in, L_busq, L_col, corte_or=True):
    # Funcion que devuelve un subset de M_in que cumpla con criterios de L_busq y L_col
    # Si corte_or es True, se asegura que alguno de los elementos L_busq se encuentre en la columna correspondiente de L_col para cada elemento de M_in
    # Si corte_or es False, se asegura que todos los elementos L_busq se encuentren en la columna correspondiente de L_col para cada elemento de M_in

    # Inicializo la matriz que se devuelve
    M_out = []; 
    # Reviso que L_busq tenga el mismo largo que L_col
    if len(L_busq) != len(L_col):
        print('ERROR: L_busq y L_col de largos distintos. Se devuelve matriz vacia.')
    # Si el largo de L_busq y L_col es 0, devuelvo M_in
    elif len(L_busq) == 0:
        M_out = M_in; 
    # En el resto de los casos, sigo con la funcion
    else:
        # Recorro cada elemento de M_in
        for i in range(len(M_in)):
            curr_L = M_in[i]; 
            # Inicializo el booleano que determina si devuelvo curr_L (depende de corte_or)
            devolver_L = not corte_or; 
            # Recorro cada elemento de L_busq y L_col
            for b in range(len(L_busq)):
                busq = L_busq[b]; 
                col = L_col[b]; 
                # Veo que curr_L[col] sea igual a busq
                if (curr_L[col] == busq) == corte_or:
                    # Si corte_or es True y curr_L[col] == busq, hay que setear devolver_L en corte_or
                    # Si corte_or es False y curr_L[col] != busq, hay que setear devolver_L en corte_or
                    devolver_L = corte_or; 
            # Despues de pasar por todos los elementos de L_busq y L_col, veo si devolver_L es True
            if devolver_L:
                M_out.append(curr_L[:]); 
    return M_out


def ver_tf_peak(L_peak, L_col, val_true='X'):
    # Funcion que busca en posiciones dentro de L_col para ver si L_peak[col]==val_true
    # Devuelve True si cualquiera de los valores es igual a val_true

    # Inicializo el booleano que se devuelve
    tf_peak = False; 
    # Recorro L_col
    for col in L_col:
        # Veo si el valor de L_peak en la columan col es igual a val_true
        if L_peak[col]==val_true:
            tf_peak = True; 
    return tf_peak


def _main():
    # Funcion para probar funciones en ejecucion del archivo

    # Defino los archivos usados
    dir_csv = path_out_tablas_main; 
    # Atributos GO: go_id, go_linkage_type, goslim_goa_accession, goslim_goa_description
    attributes = ['ensembl_gene_id', 'go_id', 'goslim_goa_accession', 'goslim_goa_description']; 
    nom_csv_mouse = L_nombres_main[1]; 
    nom_csv_human = L_nombres_main[0]; 
    genome_name_mouse = 'mm9'; 
    genome_name_human = 'hg19'; 
    # Genero los diccionarios para raton y humano
    print('>Iniciando procesamiento en raton.')
    dict_go_mouse, dict_go_gene_mouse = pipeline_go_terms(nom_csv_mouse, attributes, genome_name_mouse, dir_csv=dir_csv); 
    dict_go_gene_mouse = guardar_dict(dict_go_gene_mouse, 'MouseGOdict', path_out=path_out_go); 
    print()
    print('>Iniciando procesamiento en humano.')
    dict_go_human, dict_go_gene_human = pipeline_go_terms(nom_csv_human, attributes, genome_name_human, dir_csv=dir_csv); 
    dict_go_gene_human = guardar_dict(dict_go_gene_human, 'HumanGOdict', path_out=path_out_go); 
    '''# Defino nombre del archivo
    id_main = 1; 
    nom_arch_main = L_nombres_main[id_main]; 
    especie = L_especies_main[id_main]; 
    # Variables asociadas a especie
    if especie.lower() in ['humano']:
        nom_nkx25 = 'NKX25'; 
        nom_tbx20 = 'TBX20'; 
        nom_meis1 = 'MEIS1'; 
        nom_tgif1 = 'TGIF1'; 
        nom_hand1 = 'HAND1'; 
        nom_maf = 'MAF'; 
        nom_gata1 = 'GATA1'; 
        nom_gata6 = 'GATA6'; 
        nom_gata4 = 'GATA4'; 
    elif especie.lower() in ['raton']:
        nom_nkx25 = 'Nkx25'; 
        nom_tbx20 = 'Tbx20'; 
        nom_meis1 = 'Meis1'; 
        nom_tgif1 = 'Tgif1'; 
        nom_hand1 = 'Hand1'; 
        nom_maf = 'Maf'; 
        nom_gata1 = 'Gata1'; 
        nom_gata6 = 'Gata6'; 
        nom_gata4 = 'Gata4'; 
    else:
        print('ERROR: Especie ' + str(especie) + ' no reconocida.')
        return
    # Defino la lista de nombres de factores de transcripcion a leer
    L_nom_tf = [nom_tbx20, nom_meis1, nom_tgif1, nom_hand1, nom_maf, nom_gata1, nom_gata6, nom_gata4]; 
    # Corro generar_tabla_procesada
    M_csv = generar_tabla_procesada(nom_arch_main, nom_arch_main + '_un_min', nom_nkx25, L_nom_tf, path_in=path_in_main, path_out=path_out_tablas_main); '''

    return dict_go_human


def _abrir_pssm(nom_arch, path_arch='', solo_pssm=True, pseudocounts=0.5):
    # Extrae la matriz de counts de nom_arch y la pasa a PSSM

    # Extraigo el motif con extraer_motif_pcm()
    m = _extraer_motif_pcm(nom_arch, path_arch=path_arch); 
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


def _extraer_motif_pcm(nom_arch, path_arch=''):
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


def _generar_sitios_posibles(seq, target_len):
    # Genera todos los sitios posibles de largo target_len que incluyan seq, de manera recursiva

    # Si target_len es menor al largo de seq, tiro error y devuelvo seq
    if target_len < len(seq):
        print('ERROR: Se dio target_len menor a len(seq), se devuelve una lista con seq')
        # Lista que se devuelve
        L_out = [seq]; 
    # Si target_len es igual al largo de seq, devuelvo seq
    elif target_len == len(seq):
        # Lista que se devuelve
        L_out = [seq]; 
    # Si target_len es mayor que el largo de seq, agrego una lista a cada lado de todas las secuencias posibles
    elif target_len > len(seq):
        # Lista que se devuelve
        L_out = []; 
        # Diferencia entre target_len y len(seq)
        len_diff = target_len - len(seq); 
        # Genero todas las secuencias posibles de largo len_diff
        alphabet = ['A', 'C', 'T', 'G']; 
        L_curr = ['']; 
        L_possible = []; 
        for i in range(len_diff): 
            # Recorro cada elemento de L_curr
            for j in L_curr:
                # Agrego cada letra de las posibles a cada elemento de L_curr
                for k in alphabet:
                    L_possible.append(j + k); 
            # Una vez terminado el ciclo, copio L_possible sobre L_curr y reinicio L_possible
            L_curr = L_possible[:]; 
            L_possible = []; 
        # Agrego cada elemento de L_curr a cada lado de seq en cada combinacion posible
        for left in L_curr:
            for right in L_curr:
                L_out.append(left + seq + right); 
    return L_out


def _print_valores_venn(nom_arch, path_in, nom_tf='NKX25', hacer_prints=True):
    # Funcion para printear valores relevantes para diagrama de venn de archivo nom_arch

    # Abro el archivo 
    M_csv, L_headers = abrir_csv_headers(nom_arch, path_arch=path_in); 
    # Si hacer_prints es true, uso _print_un_valor_venn() para hacer el conteo
    if hacer_prints:
        print()
        # Totales
        text_tot = 'totales'; 
        _print_un_valor_venn(M_csv, L_headers, text_tot, [], []); 
        n_tot = contar_genes_cercanos_matriz(M_csv, L_headers, titulo_col='genes_rnaseq', devolver_L=False); 
        print('Numero de genes cercanos a peaks: ' + str(n_tot))
        print()
        # Con sitios de union
        text_su = 'con sitios de union para nkx25'; 
        _print_un_valor_venn(M_csv, L_headers, text_su, [nom_tf+'_su_confirmado'], ['X']); 
        col_su = L_headers.index(nom_tf+'_su_confirmado'); 
        L_col_su = [col_su]; 
        M_su = subset_M(M_csv, ['X'], L_col_su, corte_or=False); 
        n_su = contar_genes_cercanos_matriz(M_su, L_headers, titulo_col='genes_rnaseq', devolver_L=False); 
        print('Numero de genes cercanos a peaks: ' + str(n_su))
        print()
        # Con genes cerca
        text_genes_cerca = 'con genes cerca'; 
        _print_un_valor_venn(M_csv, L_headers, text_genes_cerca, ['genes_cerca'], ['X']); 
        col_gc = L_headers.index('genes_cerca'); 
        L_col_gc = [col_gc]; 
        M_gc = subset_M(M_csv, ['X'], L_col_gc, corte_or=False); 
        n_gc = contar_genes_cercanos_matriz(M_gc, L_headers, titulo_col='genes_rnaseq', devolver_L=False); 
        print('Numero de genes cercanos a peaks: ' + str(n_gc))
        print()
        # Con sitios de union y genes cerca
        text_su_genes_cerca = 'con sitios de union para nkx25 y genes cerca'; 
        _print_un_valor_venn(M_csv, L_headers, text_su_genes_cerca, ['genes_cerca',nom_tf+'_su_confirmado'], ['X','X']); 
        L_col_su_gc = [col_gc, col_su]; 
        M_su_gc = subset_M(M_csv, ['X','X'], L_col_su_gc, corte_or=False); 
        n_su_gc = contar_genes_cercanos_matriz(M_su_gc, L_headers, titulo_col='genes_rnaseq', devolver_L=False); 
        print('Numero de genes cercanos a peaks: ' + str(n_su_gc))
        print()
        # Con genes confirmados cerca
        text_genes_confirmados = 'con genes confirmados cerca'; 
        _print_un_valor_venn(M_csv, L_headers, text_genes_confirmados, ['genes_rnaseq_cerca'], ['X']); 
        col_conf = L_headers.index('genes_rnaseq_cerca'); 
        L_col_conf = [col_conf]; 
        M_conf = subset_M(M_csv, ['X'], L_col_conf, corte_or=False); 
        n_conf = contar_genes_cercanos_matriz(M_conf, L_headers, titulo_col='genes_rnaseq', devolver_L=False); 
        print('Numero de genes cercanos a peaks: ' + str(n_conf))
        print()
        # Con sitios de union y genes confirmados cerca
        text_su_genes_confirmados = 'con sitios de union para nkx25 y genes confirmados cerca'; 
        _print_un_valor_venn(M_csv, L_headers, text_su_genes_confirmados, ['genes_rnaseq_cerca',nom_tf+'_su_confirmado'], ['X','X']); 
        L_col_su_conf = [col_conf, col_su]; 
        M_su_conf = subset_M(M_csv, ['X','X'], L_col_su_conf, corte_or=False); 
        n_su_conf = contar_genes_cercanos_matriz(M_su_conf, L_headers, titulo_col='genes_rnaseq', devolver_L=False); 
        print('Numero de genes cercanos a peaks: ' + str(n_su_conf))
        print()
    return M_csv, L_headers


def _print_un_valor_venn(M_in, L_head, print_title, L_col_name, L_col_values):
    # Funcion para contar la cantidad de ocurrencia de L_col_values en L_col_name de M_in (buscando con L_head)
    # Busca que todas las columnas (L_col_name) correspondan con todos los valores (L_col_values)
    # Hace un print al final usando print_title como nombre del peak

    # Inicializo la lista de columnas buscadas
    L_col_search = []; 
    # Cargo los index de acuerdo a L_head y L_col_name
    for col_name in L_col_name:
        index_col = L_head.index(col_name); 
        # Agrego index_col a L_col_search
        L_col_search.append(int(index_col)); 
    # Veo que L_col_name y L_col_values tengan el mismo largo
    if len(L_col_name) != len(L_col_values):
        print('ERROR: L_col_name y L_col_value tienen que tener el mismo largo.')
        n_peaks = 0; 
    # Si tienen el mismo valor, hago el conteo
    else:
        # Cuento con buscar_L_elementos_en_M()
        n_peaks = buscar_L_elementos_en_M(M_in, L_col_values, L_col_search); 
        # Hago el print
        print('Peaks ' + str(print_title) + ': ' + str(n_peaks))
    return n_peaks


def _prueba_scores_L_sitios():
    # Funcion para probar los scores en PSSM de los sitios confirmados experimentalmente para NKX2-5
    # Abandonada hasta despues de congreso (confirmacion poco importante)

    dict_output = {}; 
    # Sitios confirmados para NKX2-5
    L_sitios = ['GCAAGTG', 'GGAAGTG', 'GAAAGTG', 'ATAAGTG', 'GTAAGTG', 'CTAAGTG', 'TCAAGTG', 'TGAAGTG', 'TAAAGTG', 'TTAAGTG']; 
    # Nombre de los archivos con las matrices de pesos para NKX2-5 en humano y raton
    nom_pcm_mouse = 'NKX25_MOUSE.H11MO.0.A.pcm'; 
    nom_pcm_human = 'NKX25_HUMAN.H11MO.0.B.pcm'; 
    # Obtengo las matrices pssm para humano y raton
    pssm_nkx25_mouse = _abrir_pssm(nom_pcm_mouse, path_arch=path_pwm_mouse, solo_pssm=True); 
    pssm_nkx25_human = _abrir_pssm(nom_pcm_human, path_arch=path_pwm_human, solo_pssm=True); 
    # Corro para humano y raton
    for especie in ['humano', 'raton']:
        # Defino la matriz pssm usada
        if especie=='humano':
            curr_pssm_nkx25 = pssm_nkx25_human; 
        elif especie=='raton':
            curr_pssm_nkx25 = pssm_nkx25_mouse; 
        else:
            print('ERROR: especie "' + str(especie) + '" no reconocida.')
        # Defino el largo para la matriz
        len_pssm = len(curr_pssm_nkx25.consensus); 
        # Defino el score maximo para la matriz y el threshold en base a eso
        max_score = curr_pssm_nkx25.max; 
        score_limit = max_score * 0.9; 
        print('Especie: ' + especie)
        print('Cutoff score: ' + str(score_limit))
        # Recorro L_sitios
        for curr_sitio in L_sitios:
            dict_output[curr_sitio+especie] = {}; 
            # Si curr_sitio es mas corto que len_pssm, puede haber problemas
            L_sitios_posibles = _generar_sitios_posibles(curr_sitio, len_pssm); 
            # Corro pssm.calculate() para todos los elementos en L_sitios_posibles
            for sitio_posible in L_sitios_posibles:
                max_forward = max(curr_pssm_nkx25.calculate(sitio_posible)); 
                seq_reverse = _reverse_seq(sitio_posible); 
                max_reverse = max(curr_pssm_nkx25.calculate(seq_reverse)); 
                score_max_sitio = max(max_forward, max_reverse); 
                # Calcular el maximo entre forward y reverse
                dict_output[curr_sitio+especie][sitio_posible] = [float(score_max_sitio), score_max_sitio < score_limit]; 
            ### FALTA
            # Calcular el score de curr_sitio para curr_pssm_nkx25
            ###
    return dict_output


def _reverse_seq(seq):
    # 
    dict_reverse = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}; 
    seq_out = ''; 
    for n in seq:
        seq_out = dict_reverse[n] + seq_out; 
    return seq_out


def _test_print_dict_scores():
    dict_scores = _prueba_scores_L_sitios(); 
    for key in dict_scores.keys():
        cont_menor = 0; 
        cont_mayor = 0; 
        max_score_sitio = 0; 
        min_score_sitio = 10; 
        for keykey in dict_scores[key].keys(): 
            curr_score = dict_scores[key][keykey][0]; 
            menor_que_limite = dict_scores[key][keykey][1]; 
            # Veo si curr_score es mayor/menor que el maximo/minimo por sitio
            if curr_score > max_score_sitio:
                max_score_sitio = float(curr_score); 
            if curr_score < min_score_sitio:
                min_score_sitio = float(curr_score); 
            # Cuento cuantos scores son mayores o menores que el limite
            if menor_que_limite:
                cont_menor += 1; 
            else:
                cont_mayor += 1; 
        print('### Sitio actual: ' + str(key))
        print('Score maximo: ' + str(max_score_sitio))
        print('Score minimo: ' + str(min_score_sitio))
        print('Scores menores que el limite: ' + str(cont_menor))
        print('Scores mayores que el limite: ' + str(cont_mayor))
        print()
    return dict_scores


def _venn_v2(nom_usado, especie):
    # Defino variable de testeo
    nom_test = nom_usado; 
    if especie == 'humano':
        tf_test = 'NKX25'; 
        nom_tbx20 = 'TBX20'; 
        nom_meis1 = 'MEIS1'; 
        nom_tgif1 = 'TGIF1'; 
        nom_hand1 = 'HAND1'; 
        nom_maf = 'MAF'; 
        nom_gata1 = 'GATA1'; 
        nom_gata6 = 'GATA6'; 
        nom_gata4 = 'GATA4'; 
    elif especie == 'raton':
        tf_test = 'Nkx25'; 
        nom_tbx20 = 'Tbx20'; 
        nom_meis1 = 'Meis1'; 
        nom_tgif1 = 'Tgif1'; 
        nom_hand1 = 'Hand1'; 
        nom_maf = 'Maf'; 
        nom_gata1 = 'Gata1'; 
        nom_gata6 = 'Gata6'; 
        nom_gata4 = 'Gata4'; 
    else:
        print('ERROR: Especie ' + str(especie) + ' no reconocida.')
        return
    # Uso _print_valores_venn() para nom_test
    M_csv, L_headers = _print_valores_venn(nom_test, path_in_main, nom_tf=tf_test, hacer_prints=False); 

    # Defino M_conf
    col_class = L_headers.index('peak_class')
    L_busq_conf = ['1-Regulado directamente', '2-Regulado indirectamente']; 
    L_col_conf = [col_class, col_class]; 
    M_conf = subset_M(M_csv, L_busq_conf, L_col_conf, corte_or=True); 
    # Total confirmados por RNA-seq
    text_conf = 'confirmados por RNA-seq'; 
    _print_un_valor_venn(M_conf, L_headers, text_conf, [], []); 
    n_conf = contar_genes_cercanos_matriz(M_conf, L_headers, titulo_col='genes_rnaseq', devolver_L=False); 
    print('Numero de genes cercanos a peaks: ' + str(n_conf))
    print()
    # Sin sitio de union
    text_sin_su = 'sin sitio de union nkx25'; 
    _print_un_valor_venn(M_conf, L_headers, text_sin_su, [tf_test+'_su_confirmado'], ['']); 
    col_su = L_headers.index(tf_test+'_su_confirmado'); 
    L_col_su = [col_su]; 
    L_busq_sin_su = ['']; 
    M_sin_su = subset_M(M_conf, L_busq_sin_su, L_col_su, corte_or=False); 
    n_sin_su = contar_genes_cercanos_matriz(M_sin_su, L_headers, titulo_col='genes_rnaseq', devolver_L=False); 
    print('Numero de genes cercanos a peaks: ' + str(n_sin_su))
    print()
    # Con sitio de union
    text_con_su = 'con sitio de union nkx25'; 
    _print_un_valor_venn(M_conf, L_headers, text_con_su, [tf_test+'_su_confirmado'], ['X']); 
    L_busq_con_su = ['X']; 
    M_con_su = subset_M(M_conf, L_busq_con_su, L_col_su, corte_or=False); 
    n_con_su = contar_genes_cercanos_matriz(M_con_su, L_headers, titulo_col='genes_rnaseq', devolver_L=False); 
    print('Numero de genes cercanos a peaks: ' + str(n_con_su))
    print()
    # Defino M_conf_con_su_otros_tf y M_conf_sin_su_otros_tf
    col_tbx20 = L_headers.index(nom_tbx20 + '_pssm'); 
    col_meis1 = L_headers.index(nom_meis1 + '_pssm'); 
    col_tgif1 = L_headers.index(nom_tgif1 + '_pssm'); 
    col_hand1 = L_headers.index(nom_hand1 + '_pssm'); 
    col_maf = L_headers.index(nom_maf + '_pssm'); 
    col_gata1 = L_headers.index(nom_gata1 + '_pssm'); 
    col_gata6 = L_headers.index(nom_gata6 + '_pssm'); 
    col_gata4 = L_headers.index(nom_gata4 + '_pssm'); 
    L_col_su_otros_tf = [col_tbx20, col_meis1, col_tgif1, col_hand1, col_maf, col_gata1, col_gata6, col_gata4]; 
    L_busq_con_su_otros_tf = ['X', 'X', 'X', 'X', 'X', 'X', 'X', 'X']; 
    M_conf_con_su_otros_tf = subset_M(M_conf, L_busq_con_su_otros_tf, L_col_su_otros_tf, corte_or=True); 
    L_busq_sin_su_otros_tf = ['', '', '', '', '', '', '', '']; 
    M_conf_sin_su_otros_tf = subset_M(M_conf, L_busq_sin_su_otros_tf, L_col_su_otros_tf, corte_or=False); 
    # Total con SU otros TF
    text_conf_con_su_otros_tf = 'con SU otros TF'; 
    _print_un_valor_venn(M_conf_con_su_otros_tf, L_headers, text_conf_con_su_otros_tf, [], []); 
    n_conf_con_su_otros_tf = contar_genes_cercanos_matriz(M_conf_con_su_otros_tf, L_headers, titulo_col='genes_rnaseq', devolver_L=False); 
    print('Numero de genes cercanos a peaks: ' + str(n_conf_con_su_otros_tf))
    print()
    # Total sin SU otros TF
    text_conf_sin_su_otros_tf = 'sin SU otros TF'; 
    _print_un_valor_venn(M_conf_sin_su_otros_tf, L_headers, text_conf_sin_su_otros_tf, [], []); 
    n_conf_sin_su_otros_tf = contar_genes_cercanos_matriz(M_conf_sin_su_otros_tf, L_headers, titulo_col='genes_rnaseq', devolver_L=False); 
    print('Numero de genes cercanos a peaks: ' + str(n_conf_sin_su_otros_tf))
    print()
    # Sin sitio de union nkx25 sin SU otros TF
    text_sin_su_sin_su_otros_tf = 'sin sitio de union nkx25 y sin sitio de union otros tf'; 
    _print_un_valor_venn(M_conf_sin_su_otros_tf, L_headers, text_sin_su_sin_su_otros_tf, [tf_test+'_su_confirmado'], ['']); 
    col_su = L_headers.index(tf_test+'_su_confirmado'); 
    L_col_su = [col_su]; 
    L_busq_sin_su = ['']; 
    M_sin_su_sin_su_otros_tf = subset_M(M_conf_sin_su_otros_tf, L_busq_sin_su, L_col_su, corte_or=False); 
    n_sin_su_sin_su_otros_tf = contar_genes_cercanos_matriz(M_sin_su_sin_su_otros_tf, L_headers, titulo_col='genes_rnaseq', devolver_L=False); 
    print('Numero de genes cercanos a peaks: ' + str(n_sin_su_sin_su_otros_tf))
    print()
    # Con sitio de union nkx25 sin SU otros TF
    text_con_su_sin_su_otros_tf = 'con sitio de union nkx25 y sin sitio de union otros tf'; 
    _print_un_valor_venn(M_conf_sin_su_otros_tf, L_headers, text_con_su_sin_su_otros_tf, [tf_test+'_su_confirmado'], ['X']); 
    L_busq_con_su = ['X']; 
    M_con_su_sin_su_otros_tf = subset_M(M_conf_sin_su_otros_tf, L_busq_con_su, L_col_su, corte_or=False); 
    n_con_su_sin_su_otros_tf = contar_genes_cercanos_matriz(M_con_su_sin_su_otros_tf, L_headers, titulo_col='genes_rnaseq', devolver_L=False); 
    print('Numero de genes cercanos a peaks: ' + str(n_con_su_sin_su_otros_tf))
    print()
    # Sin sitio de union nkx25 con SU otros TF
    text_sin_su_con_su_otros_tf = 'sin sitio de union nkx25 y con sitio de union otros tf'; 
    _print_un_valor_venn(M_conf_con_su_otros_tf, L_headers, text_sin_su_con_su_otros_tf, [tf_test+'_su_confirmado'], ['']); 
    col_su = L_headers.index(tf_test+'_su_confirmado'); 
    L_col_su = [col_su]; 
    L_busq_sin_su = ['']; 
    M_sin_su_con_su_otros_tf = subset_M(M_conf_con_su_otros_tf, L_busq_sin_su, L_col_su, corte_or=False); 
    n_sin_su_con_su_otros_tf = contar_genes_cercanos_matriz(M_sin_su_con_su_otros_tf, L_headers, titulo_col='genes_rnaseq', devolver_L=False); 
    print('Numero de genes cercanos a peaks: ' + str(n_sin_su_con_su_otros_tf))
    print()
    # Con sitio de union nkx25 con SU otros TF
    text_con_su_con_su_otros_tf = 'con sitio de union nkx25 y con sitio de union otros tf'; 
    _print_un_valor_venn(M_conf_con_su_otros_tf, L_headers, text_con_su_con_su_otros_tf, [tf_test+'_su_confirmado'], ['X']); 
    L_busq_con_su = ['X']; 
    M_con_su_con_su_otros_tf = subset_M(M_conf_con_su_otros_tf, L_busq_con_su, L_col_su, corte_or=False); 
    n_con_su_con_su_otros_tf = contar_genes_cercanos_matriz(M_con_su_con_su_otros_tf, L_headers, titulo_col='genes_rnaseq', devolver_L=False); 
    print('Numero de genes cercanos a peaks: ' + str(n_con_su_con_su_otros_tf))
    return M_csv


#################################### RESULTADOS ###################################


output_dump = []; 

if __name__=='__main__':
    output_dump.append(_main()); 


