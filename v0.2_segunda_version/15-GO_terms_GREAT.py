
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
Funciones para seleccionar peaks de ChIP-seq cerca de genes up o downregulados en ensayos de RNA-seq
Busca genes segun criterio GREAT (gen mas cercano a cada direccion del peak)

Incluye pipeline para procesar listas de genes GREAT para agregar nombres de los genes
Se puede agregar otra info a la tabla procesada

### Pseudocodigo
# Lista de genes y si estan up- o down-regulados: funcion abrir_dict_entrez_ensembl()
# Abrir archivo .bed (definir nombre)
# Agarrar cada uno de los intervalos
    # Ver genes cercanos: funcion genes_cercanos()
    # Si alguno de los genes aparecen en lista de up- o down-regulados, quedarse con el intervalo
'''

#################################### VARIABLES ####################################

# Direcciones de archivos y outputs para ib3 y casa
path_dropbox_casa = 'D:\\Users\\Admin\\Dropbox\\'; 
path_dropbox_ib3 = 'C:\\Users\\emili\\Dropbox\\'; 
path_fasta_casa = 'D:\\Archivos doctorado\\Genomas\\'; 
path_fasta_ib3 = 'X:\\Genomas\\'; 
path_output_dump_casa = 'D:\\Archivos doctorado\\Output_dump\\'; 
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
path_rnaseq_main = path_dropbox_main + 'Doctorado\\3-Genes transactivados rio abajo\\1-Interpretacion de datos ChIP-seq\\'; 
path_in_main = path_rnaseq_main; 
path_in_great = path_rnaseq_main + 'GREAT\\'; 
path_in_L_genes = path_rnaseq_main + 'gene_ontology\\listas_genes_GREAT\\'; 
path_out_main = path_output_dump_main + 'tablas_GREAT\\bed\\'; 
# Path que dependen de la especie
path_pwm_human = path_dropbox_main + 'Doctorado\\3-Genes transactivados rio abajo\\1-Interpretacion de datos ChIP-seq\\PWM_human\\'; 
path_pwm_mouse = path_dropbox_main + 'Doctorado\\3-Genes transactivados rio abajo\\1-Interpretacion de datos ChIP-seq\\PWM_mouse\\'; 

# Valores para Anderson y Dupays
bed_anderson = 'anderson_original'; 
bed_dupays = 'dupays_original'; 
nom_anderson = '9-Anderson2018_RNAseq'; 
nom_dupays = '9-Dupays2015_RNAseq'; 
nom_dict_entrez_ensembl_anderson = '9-dict_entrez_ensembl_human'; 
nom_dict_entrez_ensembl_dupays = '9-dict_entrez_ensembl_mouse'; 
col_log_fc_anderson = 3; 
col_log_fc_dupays = 2; 
col_entrez_id_anderson = 1; 
col_entrez_id_dupays = 0; 

# Listas de valores para recorrer
L_dist = [1000*1000, 100*1000, 10*1000]; 
L_nom_dist = ['1mpb', '100kpb', '10kpb']; 
L_L_updown = [['up', 'down'], ['up'], ['down']]; 

# Variables para listas de genes GREAT
L_nombres_genes_great = os.listdir(path_in_L_genes); 
L_genomas_genes_great = [hg19]*9 + [mm9]*9; 


#################################### FUNCIONES ####################################


def abrir_csv(nom_arch, dir_arch='', ext='.csv', sep=';', ignore_headers=True):
    # Abre un archivo en formato de matriz con filas por lineas y columnas separadas por sep
    # Devuelve una matriz con una lista de filas, cada una con una lista de columnas
    # ignore_headers ignora la primer fila si es True, asume que no hay titulos si es False

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


def abrir_dict_entrez_ensembl(nom_dict, path_dict='.\\', ext='.csv', sep=';', sub_sep=','):
    # Funcion para abrir un diccionario en formato pensado para traduccion entrez a ensembl
    # Cada key (entrez_id) puede tener mas de un elemento asociado (ensembl_id)

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
            # Reviso si curr_L[0] ya esta entre las keys del diccionario
            elif curr_L[0] in dict_out.keys():
                # Si ya esta, agrego curr_L[1] a los elementos en dict[key]
                print('ID ' + str(curr_L[0]) + ' repetido.')
                dict_out[curr_L[0]] = dict_out[curr_L[0]] + curr_L[1].split(sub_sep); 
            else:
                # Si la key no esta, la agrego con curr_L[1]
                dict_out[curr_L[0]] = curr_L[1].split(sub_sep); 
    return dict_out


def abrir_genes_rnaseq(nom_arch, col_entrez_id, col_log_fc, genoma, nom_dict_entrez_ensembl, path_arch='', path_dict='.\\', col_ensembl_id=-1, sub_sep='|'):
    # Funcion para abrir resultados de RNA-seq y conseguir matriz con genes y si estan up o down regulados
    # nom_arch es el nombre del archivo con los datos de RNA-seq
    # log_fc negativo significa que NKX2-5 upregula, log_fc positivo significa que NKX2-5 downregula
    # genoma es necesario para revisar el pasaje de entrez id a ensembl id

    ### Display
    cont_traducciones_intentadas = 0; 
    cont_traducciones_exitosas = 0; 
    ###
    # Inicializo la matriz que se devuelve
    M_out = []; 
    # Abro nom_arch en path_arch
    M_rnaseq = abrir_csv(nom_arch, dir_arch=path_arch, ext='.csv', sep=';', ignore_headers=True); 
    # Abro el diccionario entrez_ensembl en path_dict
    dict_entrez_ensembl = abrir_dict_entrez_ensembl(nom_dict_entrez_ensembl, path_dict=path_dict); 
    # Recorro M_rnaseq
    for i in range(len(M_rnaseq)):
        curr_row = M_rnaseq[i]; 
        # Selecciono informacion de expresion diferencial
        log_fc = float(curr_row[col_log_fc]); 
        # Defino si el gen esta upregulado o downregulado en base a log_fc
        if log_fc < 0:
            regulation = 'up'; 
        elif log_fc > 0:
            regulation = 'down'; 
        else:
            print('WARNING: Valores de log fold change no pueden valer 0.')
            regulation = 'ERROR'; 
        # Selecciono informacion de identificadores
        entrez_id = str(curr_row[col_entrez_id]); 
        # Busco ensembl IDs directamente en la tabla solo si col_ensembl_id es mayor a 0
        if col_ensembl_id>=0:
            # Agarro una lista usando split, puede generar elementos vacios en la lista
            L_ensembl_ids_raw = curr_row[col_ensembl_id].split(sub_sep); 
            # Genero la lista final
            L_ensembl_ids = []; 
            # Reviso cada elemento de L_ensembl_ids_raw
            for j in range(len(L_ensembl_ids_raw)):
                # Solo registro los elementos que no sean strings vacios
                if L_ensembl_ids_raw[j]!='':
                    L_ensembl_ids.append(L_ensembl_ids_raw[j]); 
        else:
            L_ensembl_ids = []; 
        # Si L_ensembl_ids tiene largo 0, pruebo con entrez_id
        if len(L_ensembl_ids)==0:
            # Consigo la lista de Ensembl IDs asociados a entrez_id usando traducir_entrez_to_ensembl()
            L_ensembl_ids = traducir_entrez_to_ensembl(entrez_id, genoma, dict_entrez_ensembl); 
            ### Display
            cont_traducciones_intentadas+=1; 
            if len(L_ensembl_ids)>0:
                cont_traducciones_exitosas+=1; 
            ###
        # Recorro L_ensembl_ids y los cargo en M_out con tupla junto a regulation
        for k in range(len(L_ensembl_ids)):
            curr_ensembl_id = L_ensembl_ids[k]; 
            # Agrego la tupla a M_out
            M_out.append((str(curr_ensembl_id), str(regulation))); 
    ### Display
    print('Cantidad de traducciones intentadas de entrez a ensembl: ' + str(cont_traducciones_intentadas))
    print('Cantidad de traducciones exitosas de entrez a ensembl: ' + str(cont_traducciones_exitosas))
    ###
    return M_out


def definir_pos0(gene_start, gene_end, direction):
    # Recibe el start, end y strand de un elemento gen de EnsemblRelease y devuelve el +1

    # Veo si direction es + o -
    if direction == '+':
        ret = gene_start; 
    elif direction == '-':
        ret = gene_end; 
    # Si no es ninguno de los dos, devuelvo gene_start por defecto
    else:
        print('WARNING: Strand ' + str(direction) + ' no reconocido. Se devuelve gene_start.')
        ret = gene_start; 
    return ret


def eliminar_repetidos(L_rep):
    # Funcion que elimina elementos repetidos en una lista
    # Funciona para matrices

    # Inicializo la lista que se devuelve
    L_out = []; 
    # Recorro L_rep
    for i in range(len(L_rep)):
        curr_elem = L_rep[i]; 
        # Veo si curr_elem esta en L_out
        if not (curr_elem in L_out):
            # Si no esta, lo agrego
            L_out.append(copy.deepcopy(curr_elem)); 
    return L_out


def genes_cercanos(intervalo, genoma, dist_max):
    # Obtiene los genes mas cerca de un dado intervalo, a cada lado
    # Calcula segun distancia al +1, si un +1 cae dentro del intervalo se agrega a la lista y se devuelven mas de 2 genes
    # Intervalo en formato (chr_n, pos_ini, pos_end)
    # Genoma tiene que ser elemento EnsemblRelease()

    # Inicializo la lista que se devuelve
    L_out = []; 
    # Defino contig en base a chr_n en intervalo
    contig = intervalo[0][3:]; 
    # Busco genes cercanos con genoma.genes_at_locus(contig, ini, end); 
    L_genes_upstream = genoma.genes_at_locus(contig, int(intervalo[1])-dist_max, int(intervalo[1])); 
    L_genes_in_range = genoma.genes_at_locus(contig, int(intervalo[1]), int(intervalo[2])); 
    L_genes_downstream = genoma.genes_at_locus(contig, int(intervalo[2]), int(intervalo[2])+dist_max); 
    ## Busco genes dentro del intervalo
    # Recorro L_genes_in_range
    for i in range(len(L_genes_in_range)):
        curr_gen = L_genes_in_range[i]; 
        # Defino el +1 de curr_gen con definir_pos0()
        pos0_curr_gen = definir_pos0(curr_gen.start, curr_gen.end, curr_gen.strand); 
        # Veo si pos0_curr_gen cae dentro del intervalo
        if pos0_curr_gen <= int(intervalo[2]) and pos0_curr_gen >= int(intervalo[1]):
            # Si cae dentro del intervalo, lo agrego a L_out directamente
            L_out.append([str(curr_gen.gene_id), str(curr_gen.gene_name)]); 
    ## Busco genes upstream
    # Inicializo las variables de gen mas cercano y distancia minima
    min_dist = dist_max+1; 
    gen_mas_cerca_up = []; 
    # Recorro L_genes_upstream
    for j in range(len(L_genes_upstream)):
        curr_gen = L_genes_upstream[j]; 
        # Defino el +1 de curr_gen con definir_pos0()
        pos0_curr_gen = definir_pos0(curr_gen.start, curr_gen.end, curr_gen.strand); 
        # Veo si pos0_curr_gen cae dentro del intervalo
        if pos0_curr_gen <= int(intervalo[1]) and pos0_curr_gen >= int(intervalo[1])-dist_max:
            # Si cae dentro del intervalo, calculo la distancia del gen a intervalo[1]
            curr_dist = int(intervalo[1]) - pos0_curr_gen; 
            # Veo si curr_dist es menor que min_dist
            if curr_dist < min_dist:
                # Si curr_dist es menor que min_dist, actualizo min_dist y gen_mas_cerca_up
                min_dist = int(curr_dist); 
                gen_mas_cerca_up = [str(curr_gen.gene_id), str(curr_gen.gene_name)]; 
            elif curr_dist == min_dist:
                print('WARNING: Dos genes tienen la misma min_dist=' + str(curr_dist) + ', no se carga el gen actual.')
                print('Gen actual: ' + str(curr_gen))
    # Una vez recorrido L_genes_upstream, veo si se cargo algun gen como mas cercano
    if len(gen_mas_cerca_up):
        # Si se cargo, lo agrego a L_out
        L_out.append(gen_mas_cerca_up); 
    ## Busco genes downstream
    # Inicializo las variables de gen mas cercano y distancia minima
    min_dist = dist_max+1; 
    gen_mas_cerca_down = []; 
    # Recorro L_genes_downstream
    for k in range(len(L_genes_downstream)):
        curr_gen = L_genes_downstream[k]; 
        # Defino el +1 de curr_gen con definir_pos0()
        pos0_curr_gen = definir_pos0(curr_gen.start, curr_gen.end, curr_gen.strand); 
        # Veo si pos0_curr_gen cae dentro del intervalo
        if pos0_curr_gen <= int(intervalo[2])+dist_max and pos0_curr_gen >= int(intervalo[2]):
            # Si cae dentro del intervalo, calculo la distancia del gen a intervalo[2]
            curr_dist = pos0_curr_gen - int(intervalo[2]); 
            # Veo si curr_dist es menor que min_dist
            if curr_dist < min_dist:
                # Si curr_dist es menor que min_dist, actualizo min_dist y gen_mas_cerca_down
                min_dist = int(curr_dist); 
                gen_mas_cerca_down = [str(curr_gen.gene_id), str(curr_gen.gene_name)]; 
            elif curr_dist == min_dist:
                print('WARNING: Dos genes tienen la misma min_dist=' + str(curr_dist) + ', no se carga el gen actual.')
                print('Gen actual: ' + str(curr_gen))
    # Una vez recorrido L_genes_downstream, veo si se cargo algun gen como mas cercano
    if len(gen_mas_cerca_down):
        # Si se cargo, lo agrego a L_out
        L_out.append(gen_mas_cerca_down); 
    return L_out


def guardar_bed(M_bed, nom_bed, path_bed='.\\', ext='.bed', sep='\t'):
    # Funcion para guardar una matriz en un archivo formato .bed

    # Defino la direccion del archivo en base a dir_bed y nom_bed
    dirarch = os.path.join(path_bed, nom_bed + ext); 
    # Creo el archivo
    with open(dirarch, 'w') as F_out:
        print('Archivo ' + nom_bed + ext + ' creado.')
    # Lo vuelvo a abrir en modo append
    with open(dirarch, 'a') as F_out:
        # Recorro M_bed
        for L_bed in M_bed:
            curr_str = ''; 
            # Recorro L_bed
            for i in L_bed:
                curr_str = curr_str + str(i) + sep; 
            # Elimino la ultima ocurrencia de sep
            curr_str = curr_str.rstrip(sep); 
            # Agrego end of line y guardo en F_out
            F_out.write(curr_str + '\n'); 
    return M_bed


def guardar_csv(M_csv, nom_out, path_out='.\\', ext='.csv', sep=';'):
    # Funcion para guardar una matriz en formato .csv

    # Defino la direccion del archivo en base a path_out y nom_our
    dirarch = os.path.join(path_out, nom_out + ext); 
    # Creo el archivo
    with open(dirarch, 'w') as F_out:
        print('Archivo ' + nom_out + ext + ' creado.')
    # Lo vuelvo a abrir en modo append
    with open(dirarch, 'a') as F_out:
        # Recorro M_csv
        for L_csv in M_csv:
            curr_str = ''; 
            # Recorro L_csv
            for i in L_csv:
                curr_str = curr_str + str(i) + sep; 
            # Elimino la ultima ocurrencia de sep
            curr_str = curr_str.rstrip(sep); 
            # Agrego end of line y guardo en F_out
            F_out.write(curr_str + '\n'); 
    return M_csv


def pipeline_genes_great(L_nom_arch, L_genomas, path_arch='', path_out='.\\'):
    # Funcion para abrir archivos con ID de Ensembl de genes up/down regulados
    # Los vuelve a guardar en path_out con info agregada

    # Recorro L_nom_arch y L_genomas
    for i in range(len(L_nom_arch)):
        # Selecciono el nombre de archivo y la extension
        curr_arch = L_nom_arch[i][:-4]; 
        curr_ext = L_nom_arch[i][-4:]; 
        curr_genoma = L_genomas[i]; 
        # Extraigo la matriz del archivo
        M_arch = abrir_csv(curr_arch, dir_arch=path_arch, ext=curr_ext, ignore_headers=False); 
        print(M_arch)
        # Elimino elementos repetidos
        M_arch = eliminar_repetidos(M_arch); 
        print(M_arch)
        # Proceso la matriz con procesar_arch_genes()
        M_procesada = procesar_arch_genes(M_arch, curr_genoma); 
        # Guardo la matriz con guardar_csv()
        M_procesada = guardar_csv(M_procesada, curr_arch+'_procesado', path_out=path_out); 
    return ''


def pipeline_peaks_genes_confirmados(nom_bed, nom_rnaseq, nom_dict_entrez_ensembl, nom_out, col_entrez_id, col_log_fc, genoma, dist_max, path_bed='', path_rnaseq='', 
                                     path_dict='.\\', path_out='', sep='\t', ext='.bed', L_updown=['up', 'down']):
    # Pipeline que selecciona picos de un archivo .bed de acuerdo a si alguno de los 2 genes mas cercanos esta upregulado o downregulado en RNAseq
    # L_updown define el tipo de genes que se seleccionan segun RNAseq, si solo upregulados, solo downregulados o ambos

    # Abro el archivo .bed
    M_bed = abrir_csv(nom_bed, dir_arch=path_bed, ext=ext, sep=sep, ignore_headers=False); 
    # Inicializo la matriz que se devuelve
    M_out = []; 
    # Consigo la lista de genes confirmados por RNAseq
    L_genes_rnaseq = abrir_genes_rnaseq(nom_rnaseq, col_entrez_id, col_log_fc, genoma, nom_dict_entrez_ensembl, path_arch=path_rnaseq, path_dict=path_dict); 
    # Inicializo una lista de genes RNAseq cerca de peaks
    L_genes_rnaseq_cerca = []; 
    # Recorro M_bed
    for i in range(len(M_bed)):
        curr_range = M_bed[i]; 
        # Consigo la lista de genes cercanos con funcion genes_cercanos
        L_genes_cerca = genes_cercanos((curr_range[0], curr_range[1], curr_range[2]), genoma, dist_max); 
        L_genes_cerca_id = []; 
        # Recorro L_genes_cerca para pasar solo los IDs a L_genes_cerca_id
        for gen_cerca in L_genes_cerca:
            L_genes_cerca_id.append(str(gen_cerca[0])); 
        # Recorro L_genes_rnaseq y veo si alguno de los genes en L_genes cerca aparecen
        for j in range(len(L_genes_rnaseq)):
            curr_gen_rnaseq = L_genes_rnaseq[j][0]; 
            curr_regulation = L_genes_rnaseq[j][1]; 
            # Veo que curr_regulation este en L_updown y curr_gen_rnaseq este en L_genes_cerca_id
            if (curr_gen_rnaseq in L_genes_cerca_id) and (curr_regulation in L_updown):
                # Si eso pasa, se guarda el rango actual
                M_out.append(curr_range[:]); 
                # Guardo el gen de L_genes_rnaseq encontrado en L_genes_rnaseq_cerca
                L_genes_rnaseq_cerca.append(L_genes_rnaseq[j]); 
                ### Aca solia haber un break
    # Guardo M_out en archivo .bed
    M_out = guardar_bed(M_out, nom_out, path_bed=path_out); 
    # Guardo L_genes_rnaseq_cerca con guardar_bed
    L_genes_rnaseq_cerca = guardar_bed(L_genes_rnaseq_cerca, nom_out+'_genes', path_bed=path_out, ext='.csv', sep=';'); 
    return M_out


def procesar_arch_genes(M_arch, genoma):
    # Funcion para procesar M_arch y devolver matriz con mas info sobre los genes en M_arch

    # Inicializo la matriz que se devuelve
    M_out = []; 
    # Recorro M_arch
    for i in range(len(M_arch)):
        curr_row = M_arch[i]; 
        # Defino gene_id y regulation
        gene_id = curr_row[0]; 
        regulation = curr_row[1]; 
        # Uso genoma para buscar gene_id y conseguir mas info
        gene_element = genoma.gene_by_id(gene_id); 
        #Gene(gene_id='ENSMUSG00000042401', gene_name='Crtac1', biotype='protein_coding', contig='19', start=42357527, end=42506273, strand='-', genome='NCBIM37')
        # Inicializo la lista que se carga a M_out
        L_out = []; 
        # Cargo la informacion que sea importante
        L_out.append(str(gene_id)); 
        L_out.append(str(gene_element.gene_name)); 
        L_out.append(str(regulation)); 
        # Agrego L_out a M_out
        M_out.append(L_out[:]); 
    return M_out


def traducir_entrez_to_ensembl(entrez_id, genoma, dict_entrez_ensembl):
    # Traduce una serie de IDs y trata de llegar a uno o varios Ensembl IDs

    # Lista de ids encontrado que se devuelve
    L_out = []; 
    # Reviso si entrez_id se encuentra entre las keys de dict_entrez_ensembl
    if str(entrez_id) in dict_entrez_ensembl.keys():
        L_ensembl_ids = dict_entrez_ensembl[str(entrez_id)]; 
        # Recorro los ensembl_ids en L_ensembl_id
        for i in range(len(L_ensembl_ids)):
            curr_ensembl_id = L_ensembl_ids[i]; 
            # Pruebo si el ensembl_id encontrado esta en el genoma
            try:
                # Uso gene_by_id() para ver si el ensembl_id es valido
                genoma.gene_by_id(curr_ensembl_id); 
                # Si no tira error, lo cargo a L_out
                L_out.append(str(curr_ensembl_id)); 
            except:
                pass
    return L_out


def _main():
    # Funcion para probar funciones en ejecucion del archivo

    '''
    # Recorro las distintas listas de valores para correr pipeline_peaks_genes_confirmados()
    for i in range(len(L_dist)):
        dist_max_usada = L_dist[i]; 
        dist_nombre = L_nom_dist[i]; 
        for j in range(len(L_L_updown)):
            L_updown_usado = L_L_updown[j]; 
            # Defino el nombre en el archivo (down, up, updown)
            updown_nombre = ''; 
            for k in L_updown_usado:
                updown_nombre = updown_nombre+k; 
            # Corro pipeline_peaks_genes_confirmados() para humano y raton
            M_bed_anderson = pipeline_peaks_genes_confirmados(bed_anderson, nom_anderson, nom_dict_entrez_ensembl_anderson, 'anderson_'+dist_nombre+'_'+updown_nombre, 
                                                              col_entrez_id_anderson, col_log_fc_anderson, hg19, dist_max_usada, path_bed=path_in_great, 
                                                              path_rnaseq=path_rnaseq_main, path_dict=path_in_main, path_out=path_out_main, L_updown=L_updown_usado); 
            M_bed_dupays = pipeline_peaks_genes_confirmados(bed_dupays, nom_dupays, nom_dict_entrez_ensembl_dupays, 'dupays_'+dist_nombre+'_'+updown_nombre, 
                                                            col_entrez_id_dupays, col_log_fc_dupays, mm9, dist_max_usada, path_bed=path_in_great, 
                                                            path_rnaseq=path_rnaseq_main, path_dict=path_in_main, path_out=path_out_main, L_updown=L_updown_usado); 
    '''
    pipeline_genes_great(L_nombres_genes_great, L_genomas_genes_great, path_arch=path_in_L_genes, path_out=path_in_L_genes); 
    return ''


#################################### RESULTADOS ###################################

output_dump = []; 

if __name__=='__main__':
    output_dump.append(_main()); 

