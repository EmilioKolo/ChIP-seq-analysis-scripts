
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
path_analisischip_ib3 = 'C:\\Users\\emili\\Documents\\Archivos nube\\Git\\analisissecuenciaschip\\'; 
path_analisischip_main = path_analisischip_ib3; # Cambiar aca
sys.path.append(path_analisischip_main); 
from analisischip.seq import seq_data, seq_handler

# Genomas de referencia
hg19 = EnsemblRelease(75, species='human'); 
hg38 = EnsemblRelease(102, species='human'); 
mm9 = EnsemblRelease(54, species='mouse'); 
mm10 = EnsemblRelease(102, species='mouse'); 


'''
Funciones para clasificar los distintos peaks de resultados de ChIP-seq y generar archivos .fasta para MEME-ChIP

Funcion para abrir archivos .bed y guardarlos como matriz: abrir_bed()

Funcion para procesar rangos de peaks sacados de .bed: procesar_bed()

Funcion para guardar la matriz procesada en un .csv: guardar_matriz()

Funcion que corre todo en orden para generar archivo: pipeline_bed()
'''

#################################### VARIABLES ####################################

# Direcciones de archivos y outputs
path_bed_casa = 'D:\\Dropbox\\Doctorado\\3-Genes transactivados rio abajo\\0-Fuentes\\Papers ChIP-seq\\'; 
path_bed_ib3 = 'C:\\Users\\emili\\Dropbox\\Doctorado\\3-Genes transactivados rio abajo\\0-Fuentes\\Papers ChIP-seq\\'; 
path_fasta_casa = 'D:\\Archivos doctorado\\Genomas\\'; 
path_fasta_ib3 = 'X:\\Genomas\\'; 
path_out_casa = 'D:\\Archivos doctorado\\Output_dump\\PeaksClasificados\\'; 
path_out_ib3 = 'X:\\Output_dump\\PeaksClasificados\\'; 

# Nombres de archivos .bed y datos relacionados
bed_dupays = 'Dupays2015'; 
bed_anderson = 'Anderson2018-GSE89457consensus'; 
L_bed = [bed_dupays, bed_anderson]; 
nom_genoma_dupays = 'mm9'; 
genoma_dupays = mm9; 
nom_genoma_anderson = 'hg19'; 
genoma_anderson = hg19; 

# Datos para busquedas
L_dist = [1500, 10000, 50000, 100000, 1000000]; 
L_aagtg = ['AAGTG']; 
L_confirmados = ['GCAAGTG', 'GGAAGTG', 'GAAAGTG', 'ATAAGTG', 'GTAAGTG', 'CTAAGTG', 'TCAAGTG', 'TGAAGTG', 'TAAAGTG', 'TTAAGTG']; 

# Variables main()

path_bed_main = path_bed_casa; 
path_fasta_main = path_fasta_casa; 
path_out_main = path_out_casa; 

dist_main = L_dist[0]; 
L_sitios_main = L_confirmados; 


#################################### FUNCIONES ####################################


def abrir_bed(nom_arch, dir_arch='', ext='.bed', sep='\t', ignore_headers=False):
    # Abre un archivo .bed con peaks de ChIP-seq usando abrir_csv()

    # Inicializo la matriz que se devuelve
    M_out = abrir_csv(nom_arch, dir_arch=dir_arch, ext=ext, sep=sep, ignore_headers=ignore_headers); 
    return M_out


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


def buscar_genes_cerca_peak(chr_n, pos_ini, pos_end, contig, genome, dist_max, sep_sub=','):
    # Busca genes a dist_max de un peak

    # Busco genes alrededor del peak
    L_genes_cerca_raw = genome.genes_at_locus(contig, pos_ini-dist_max, pos_end+dist_max); 
    # Defino lista curada de genes cerca (por si elimino genes no protein_coding) y una lista de genes filtrados
    L_genes_cerca = []; 
    L_genes_filtrados = []; 
    # Recorro todos los genes encontrados
    for curr_gen in L_genes_cerca_raw:
        # Reviso si curr_gen esta efectivamente cerca de peak
        start_cerca_de_peak = (curr_gen.start > (pos_ini-dist_max)) and (curr_gen.start < (pos_end+dist_max)); 
        end_cerca_de_peak = (curr_gen.end > (pos_ini-dist_max)) and (curr_gen.end < (pos_end+dist_max)); 
        peak_dentro_de_gen = ((curr_gen.start < pos_ini) and (curr_gen.end > pos_end)) or ((curr_gen.end < pos_ini) and (curr_gen.start > pos_end)); 
        # Solo me quedo con los genes con start cerca del peak
        if start_cerca_de_peak:
            # Solo registro los genes con biotype protein_coding
            if curr_gen.biotype == 'protein_coding':
                L_genes_cerca.append(curr_gen); 
            # Si no es protein_coding, lo guardo en genes filtrados (por si lo anoto de alguna manera)
            else:
                L_genes_filtrados.append(curr_gen); 
        ### Display
        elif end_cerca_de_peak:
            pass
            #print()
            #print('Peak (con dist_max agregada): ' + str(chr_n) + ', ' + str(pos_ini-dist_max) + ', ' + str(pos_end+dist_max))
            #print('End cerca de peak: ' + str(curr_gen))
        elif peak_dentro_de_gen:
            pass
            #print()
            #print('Peak (con dist_max agregada): ' + str(chr_n) + ', ' + str(pos_ini-dist_max) + ', ' + str(pos_end+dist_max))
            #print('Peak dentro de gen, end/start lejos de peak: ' + str(curr_gen))
        else:
            print()
            print('Peak (con dist_max agregada): ' + str(chr_n) + ', ' + str(pos_ini-dist_max) + ', ' + str(pos_end+dist_max))
            print('ERROR RARO CON ESTE GEN: ' + str(curr_gen))
        ###
    # Reviso si hay genes encontrados o no
    if len(L_genes_cerca) > 0:
        str_genes_cerca = 'GenesCerca'; 
    elif len(L_genes_filtrados) > 0:
        str_genes_cerca = 'GenesFiltrados'; 
    else:
        str_genes_cerca = 'NoHayGenes'; 
    # Creo una lista de ids de genes cercanos
    str_genes_cerca_id = ''; 
    # Reviso los genes cerca seleccionados
    for gen_cerca in L_genes_cerca:
        str_genes_cerca_id = str_genes_cerca_id + str(gen_cerca.gene_id) + sep_sub; 
    for gen_cerca_filtrado in L_genes_filtrados:
        str_genes_cerca_id = str_genes_cerca_id + '(' + str(gen_cerca_filtrado.gene_id) + ')' + sep_sub; 
    str_genes_cerca_id = str_genes_cerca_id.rstrip(sep_sub); 
    return str_genes_cerca, str_genes_cerca_id, L_genes_cerca, L_genes_filtrados


def buscar_sitios_cerca_genes_peak(dist_max, str_genes_cerca, str_sitios_encontrados, L_genes_cerca, L_genes_filtrados, L_pos_sitios_union, sep_sub=','):
    # Busca sitios de union cerca de genes que esten cerca de un peak
    # Procesa el output de buscar_genes_cerca_peak() y buscar_sitios_peak()

    '''# Defino cada combinacion de estados
        # Genes cerca
            # Con sitios cerca de genes (1)
            # Sin sitios cerca de genes (2)
            # Sin sitios (3)
        # Genes filtrados
            # Con sitios cerca de genes (4)
            # Sin sitios cerca de genes (4)
            # Sin sitios (5)
        # Genes no cerca
            # Con sitios de union (6)
            # Sin sitios de union (7)'''
    str_estado_peak = ''; 
    L_su_cerca_genes = []; 
    # Hay genes cerca del peak y no fueron filtrados (son protein_coding)
    if str_genes_cerca == 'GenesCerca':
        # Si hay sitios de union en el peak
        if str_sitios_encontrados == 'SitiosEncontrados': 
            str_estado_peak = 'GenConSitio'; 
            # Recorro cada gen cercano
            for gen_cerca in L_genes_cerca:
                gen_start = gen_cerca.start; 
                gen_end = gen_cerca.end; 
                # Recorro cada sitio de union
                for SU_encontrado in L_pos_sitios_union:
                    pos_ini_su = SU_encontrado[0]; 
                    pos_end_su = SU_encontrado[1]; 
                    # Mido cercania por distancia a gen_start
                    dist = min(abs(gen_start-pos_ini_su), abs(gen_start-pos_end_su)); 
                    # Veo si dist es menor a dist_max
                    if dist < dist_max:
                        L_su_cerca_genes.append(gen_cerca.gene_id); 
            # Reviso si hay algun gen cerca de un sitio de union
            if len(L_su_cerca_genes) > 0:
                str_estado_peak = str_estado_peak + 'Cercano'; 
            else:
                str_estado_peak = str_estado_peak + 'Lejos'; 
        # Si no hay sitios de union en el peak
        else:
            str_estado_peak = 'GenSinSitio'; 
    # Hay genes cerca del peak pero fueron filtrados (no son protein_coding)
    elif str_genes_cerca == 'GenesFiltrados': 
        # Si hay sitios de union en el peak
        if str_sitios_encontrados == 'SitiosEncontrados': 
            str_estado_peak = 'FiltradoConSitio'; 
            # Recorro cada gen filtrado
            for gen_cerca_filtrado in L_genes_filtrados:
                gen_start = gen_cerca_filtrado.start; 
                gen_end = gen_cerca_filtrado.end; 
                # Recorro cada sitio de union
                for SU_encontrado in L_pos_sitios_union:
                    pos_ini_su = SU_encontrado[0]; 
                    pos_end_su = SU_encontrado[1]; 
                    # Mido cercania por distancia a gen_start
                    dist = min(abs(gen_start-pos_ini_su), abs(gen_start-pos_end_su)); 
                    # Veo si dist es menor a dist_max
                    if dist < dist_max:
                        L_su_cerca_genes.append(gen_cerca_filtrado.gene_id); 
            # Reviso si hay algun gen cerca de un sitio de union
            if len(L_su_cerca_genes) > 0:
                str_estado_peak = str_estado_peak + 'Cercano'; 
            else:
                str_estado_peak = str_estado_peak + 'Lejos'; 
        # Si no hay sitios de union en el peak
        else:
            str_estado_peak = 'FiltradoSinSitio'; 
    # Sin genes cerca
    else: 
        # Con sitios encontrados
        if str_sitios_encontrados == 'SitiosEncontrados':
            str_estado_peak = 'SitiosSinGen'; 
        # Sin sitios encontrados
        else:
            str_estado_peak = 'SinSitiosNiGen'; 
    # Defino str para guardar genes con sitios cerca
    str_genes_sitios_cerca = ''; 
    # Recorro L_su_cerca_genes
    for su_cerca_gen in L_su_cerca_genes:
        str_genes_sitios_cerca = str_genes_sitios_cerca + str(su_cerca_gen) + sep_sub; 
    str_genes_sitios_cerca = str_genes_sitios_cerca.rstrip(sep_sub); 
    return str_estado_peak, str_genes_sitios_cerca


def buscar_sitios_peak(chr_n, pos_ini, pos_end, genome_name, genome, L_sitios, path_fasta='', sep_sub=','):
    # Busca sitios de union en un peak

    # Inicializo el elemento seq_data que consigue las secuencias
    sequence_data = seq_data(genome_name, genome_element=genome, path_fasta=path_fasta); 
    # Consigo la secuencia del peak con seq_data._consulta_secuencia_fasta()
    seq_peak = sequence_data._consulta_secuencia_fasta(chr_n, pos_ini, pos_end); 
    # Inicializo lista de todos los sitios que esten presentes
    L_sitios_encontrados = []; 
    # Inicializo lista de posiciones de todos los sitios encontrados
    L_pos_sitios_union = []; 
    # Busco ocurrencia de cada uno de los sitios en L_sitios
    for curr_sitio in L_sitios:
        # Busco todos los sitios de union que hayan para curr_sitio
        L_curr_sitio_encontrados = sequence_data._buscar_SU_en_seq(curr_sitio, seq_peak); 
        # Si encuentra por lo menos un sitio, agrego curr_sitio a L_sitios encontrados
        if len(L_curr_sitio_encontrados) > 0:
            L_sitios_encontrados.append(str(curr_sitio)); 
        for curr_sitio_encontrado in L_curr_sitio_encontrados:
            # Defino pos_ini y pos_end respecto a chr_n
            curr_pos_ini = curr_sitio_encontrado[0]+pos_ini-1; 
            curr_pos_end = curr_sitio_encontrado[1]+pos_ini-1; 
            SU_encontrado = [curr_pos_ini, curr_pos_end, curr_sitio_encontrado[2]]; 
            L_pos_sitios_union.append(SU_encontrado[:]); 
    # Reviso si hay sitios encontrados o no
    if len(L_sitios_encontrados) > 0:
        str_sitios_encontrados='SitiosEncontrados'; 
        # Creo el string con la lista de secuencias encontradas
        str_seq_sitios = ''; 
        for seq_sitio in L_sitios_encontrados:
            # Registro cada secuencia
            str_seq_sitios = str_seq_sitios + str(seq_sitio) + sep_sub; 
        str_seq_sitios = str_seq_sitios.rstrip(sep_sub); 
    else:
        str_sitios_encontrados='NoHaySitios'; 
        str_seq_sitios = ''; 
    # Inicializo el string con la lista de posiciones de sitios
    str_pos_su = ''; 
    # Reviso las posiciones de los sitios encontrados
    for SU_encontrado in L_pos_sitios_union:
        # Registro cada posicion
        str_pos_su = str_pos_su + str(SU_encontrado[0]) + '_' + str(SU_encontrado[1]) + sep_sub; 
    str_pos_su = str_pos_su.rstrip(sep_sub); 
    # Agrego los elementos correspondientes a busqueda de secuencias a L_out
    return str_sitios_encontrados, str_seq_sitios, str_pos_su, L_pos_sitios_union


def guardar_matriz(M_out, nom_out, dir_out='.\\', ext='.csv', sep=';', L_head=[]):
    # Guarda los contenidos de la matriz M_out en el archivo nom_out ubicado en dir_out

    # Defino la direccion del archivo en base a dir_out y nom_out
    dirarch = os.path.join(dir_out, nom_out + ext); 
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
            str_head.rstrip(sep); 
            # Agrego end of line y guardo en F_out
            F_out.write(str_head + '\n'); 
        # Recorro M_out
        for L_out in M_out:
            curr_str = ''; 
            # Recorro L_out
            for i in L_out:
                curr_str = curr_str + str(i) + sep; 
            # Elimino la ultima ocurrencia de sep
            curr_str.rstrip(sep); 
            # Agrego end of line y guardo en F_out
            F_out.write(curr_str + '\n'); 
    return M_out


def pipeline_bed(nom_bed, nom_genoma, genoma, dist_max, L_sitios, nom_out='', L_col=[0,1,2], path_bed='', path_fasta='', path_out='.\\', L_header=[]):
    # Corre abrir_bed(), procesar_bed() y guardar_matriz() en orden

    # Cargo los peaks en M_bed
    M_bed = abrir_bed(nom_bed, dir_arch=path_bed); 
    # Proceso M_bed
    M_out = procesar_bed(M_bed, nom_genoma, genoma, dist_max, L_sitios, L_col=L_col, path_fasta=path_fasta); 
    # Solo guardo si nom_out tiene texto
    if nom_out != '':
        M_out = guardar_matriz(M_out, nom_out, dir_out=path_out, L_head=L_header); 
    return M_out


def procesar_bed(M_bed, genome_name, genome, dist_max, L_sitios, L_col=[0,1,2], path_fasta=''):
    # Funcion central que recorre M_bed, procesa los peaks y devuelve un output para guardar en .csv

    ### Display
    len_bed = len(M_bed); 
    cont_display = 0; 
    display_rate = 25; 
    ###
    # Inicializo la matriz que se devuelve
    M_out = []; 
    # Defino id_chr, id_ini, id_end con L_col
    id_chr, id_ini, id_end = L_col[0], L_col[1], L_col[2]; 
    # Recorro la matriz extraida con abrir_bed()
    for peak in M_bed:
        # Cargo chr_n, pos_ini y pos_end con los datos de L_col
        chr_n = peak[id_chr]; 
        pos_ini = min(int(peak[id_ini]), int(peak[id_end])); 
        pos_end = max(int(peak[id_ini]), int(peak[id_end])); 
        # Uso procesar_peak() para hacer todo lo necesario con la info del peak
        L_out = procesar_peak(chr_n, pos_ini, pos_end, genome_name, genome, dist_max, L_sitios, path_fasta=path_fasta); 
        # Agrego la info que no se agarro con L_col al final de L_out
        for i in range(len(peak)):
            if not (i in L_col):
                L_out.append(peak[i]); 
        # Agrego L_out a M_out
        M_out.append(L_out[:]); 
        ### Display
        cont_display += 1; 
        if cont_display%display_rate==0:
            print('Progreso: ' + str(cont_display) + '/' + str(len_bed))
        ###
    return M_out


def procesar_peak(chr_n, pos_ini, pos_end, genome_name, genome, dist_max, L_sitios, path_fasta=''):
    # Funcion para hacer todo lo necesario con la info del peak
    # Busca sitios de union (L_sitios), genes cercanos (genome, dist_max) y sitios de union cerca de los genes cercanos

    # Defino un separador de subdivision
    sep_sub = ','; 
    # Inicializo la lista que se devuelve
    L_out = []; 
    # Defino contig en base a chr_n
    if chr_n[:3] == 'chr':
        contig = chr_n[3:]; 
    else:
        contig = chr_n; 
    # Agrego los elementos basicos a L_out
    L_out.append(chr_n); 
    L_out.append(contig); 
    L_out.append(pos_ini); 
    L_out.append(pos_end); 

    ## Busqueda de sitios de union
    str_sitios_encontrados, str_seq_sitios, str_pos_su, L_pos_sitios_union = buscar_sitios_peak(chr_n, pos_ini, pos_end, genome_name, genome, L_sitios, path_fasta=path_fasta, sep_sub=sep_sub); 
    '''# Inicializo el elemento seq_data que consigue las secuencias
    sequence_data = seq_data(genome_name, genome_element=genome, path_fasta=path_fasta); 
    # Consigo la secuencia del peak con seq_data._consulta_secuencia_fasta()
    seq_peak = sequence_data._consulta_secuencia_fasta(chr_n, pos_ini, pos_end); 
    # Inicializo lista de todos los sitios que esten presentes
    L_sitios_encontrados = []; 
    # Inicializo lista de posiciones de todos los sitios encontrados
    L_pos_sitios_union = []; 
    # Busco ocurrencia de cada uno de los sitios en L_sitios
    for curr_sitio in L_sitios:
        # Busco todos los sitios de union que hayan para curr_sitio
        L_curr_sitio_encontrados = sequence_data._buscar_SU_en_seq(curr_sitio, seq_peak); 
        # Si encuentra por lo menos un sitio, agrego curr_sitio a L_sitios encontrados
        if len(L_curr_sitio_encontrados) > 0:
            L_sitios_encontrados.append(str(curr_sitio)); 
        for curr_sitio_encontrado in L_curr_sitio_encontrados:
            # Defino pos_ini y pos_end respecto a chr_n
            curr_pos_ini = curr_sitio_encontrado[0]+pos_ini-1; 
            curr_pos_end = curr_sitio_encontrado[1]+pos_ini-1; 
            SU_encontrado = [curr_pos_ini, curr_pos_end, curr_sitio_encontrado[2]]; 
            L_pos_sitios_union.append(SU_encontrado[:]); 
    # Reviso si hay sitios encontrados o no
    if len(L_sitios_encontrados) > 0:
        str_sitios_encontrados='SitiosEncontrados'; 
        # Creo el string con la lista de secuencias encontradas
        str_seq_sitios = ''; 
        for seq_sitio in L_sitios_encontrados:
            # Registro cada secuencia
            str_seq_sitios = str_seq_sitios + str(seq_sitio) + sep_sub; 
        str_seq_sitios = str_seq_sitios.rstrip(sep_sub); 
    else:
        str_sitios_encontrados='NoHaySitios'; 
        str_seq_sitios = ''; 
    # Inicializo el string con la lista de posiciones de sitios
    str_pos_su = ''; 
    # Reviso las posiciones de los sitios encontrados
    for SU_encontrado in L_pos_sitios_union:
        # Registro cada posicion
        str_pos_su = str_pos_su + str(SU_encontrado[0]) + '_' + str(SU_encontrado[1]) + sep_sub; 
    str_pos_su = str_pos_su.rstrip(sep_sub);'''
    # Agrego los elementos correspondientes a busqueda de secuencias a L_out
    L_out.append(str_sitios_encontrados); 
    L_out.append(str_seq_sitios); 
    L_out.append(str_pos_su); 

    ## Busqueda de genes cercanos
    str_genes_cerca, str_genes_cerca_id, L_genes_cerca, L_genes_filtrados = buscar_genes_cerca_peak(chr_n, pos_ini, pos_end, contig, genome, dist_max, sep_sub=sep_sub); 
    '''# Busco genes alrededor del peak
    L_genes_cerca_raw = genome.genes_at_locus(contig, pos_ini-dist_max, pos_end+dist_max); 
    # Defino lista curada de genes cerca (por si elimino genes no protein_coding) y una lista de genes filtrados
    L_genes_cerca = []; 
    L_genes_filtrados = []; 
    # Recorro todos los genes encontrados
    for curr_gen in L_genes_cerca_raw:
        # Reviso si curr_gen esta efectivamente cerca de peak
        start_cerca_de_peak = (curr_gen.start > (pos_ini-dist_max)) and (curr_gen.start < (pos_end+dist_max)); 
        end_cerca_de_peak = (curr_gen.end > (pos_ini-dist_max)) and (curr_gen.end < (pos_end+dist_max)); 
        peak_dentro_de_gen = ((curr_gen.start < pos_ini) and (curr_gen.end > pos_end)) or ((curr_gen.end < pos_ini) and (curr_gen.start > pos_end)); 
        # Solo me quedo con los genes con start cerca del peak
        if start_cerca_de_peak:
            # Solo registro los genes con biotype protein_coding
            if curr_gen.biotype == 'protein_coding':
                L_genes_cerca.append(curr_gen); 
            # Si no es protein_coding, lo guardo en genes filtrados (por si lo anoto de alguna manera)
            else:
                L_genes_filtrados.append(curr_gen); 
        ### Display
        elif end_cerca_de_peak:
            pass
            #print()
            #print('Peak (con dist_max agregada): ' + str(chr_n) + ', ' + str(pos_ini-dist_max) + ', ' + str(pos_end+dist_max))
            #print('End cerca de peak: ' + str(curr_gen))
        elif peak_dentro_de_gen:
            pass
            #print()
            #print('Peak (con dist_max agregada): ' + str(chr_n) + ', ' + str(pos_ini-dist_max) + ', ' + str(pos_end+dist_max))
            #print('Peak dentro de gen, end/start lejos de peak: ' + str(curr_gen))
        else:
            print()
            print('Peak (con dist_max agregada): ' + str(chr_n) + ', ' + str(pos_ini-dist_max) + ', ' + str(pos_end+dist_max))
            print('ERROR RARO CON ESTE GEN: ' + str(curr_gen))
        ###
    # Reviso si hay genes encontrados o no
    if len(L_genes_cerca) > 0:
        str_genes_cerca = 'GenesCerca'; 
    elif len(L_genes_filtrados) > 0:
        str_genes_cerca = 'GenesFiltrados'; 
    else:
        str_genes_cerca = 'NoHayGenes'; 
    # Creo una lista de ids de genes cercanos
    str_genes_cerca_id = ''; 
    # Reviso los genes cerca seleccionados
    for gen_cerca in L_genes_cerca:
        str_genes_cerca_id = str_genes_cerca_id + str(gen_cerca.gene_id) + sep_sub; 
    for gen_cerca_filtrado in L_genes_filtrados:
        str_genes_cerca_id = str_genes_cerca_id + '(' + str(gen_cerca_filtrado.gene_id) + ')' + sep_sub; 
    str_genes_cerca_id = str_genes_cerca_id.rstrip(sep_sub);'''
    # Agrego los elementos correspondientes a busqueda de genes a L_out
    L_out.append(str_genes_cerca); 
    L_out.append(str_genes_cerca_id); 

    ## Busqueda de sitios de union cerca de genes cercanos
    str_estado_peak, str_genes_sitios_cerca = buscar_sitios_cerca_genes_peak(dist_max, str_genes_cerca, str_sitios_encontrados, L_genes_cerca, L_genes_filtrados, L_pos_sitios_union, sep_sub=sep_sub); 
    '''str_estado_peak = ''; 
    L_su_cerca_genes = []; 
    # Hay genes cerca del peak y no fueron filtrados (son protein_coding)
    if str_genes_cerca == 'GenesCerca':
        # Si hay sitios de union en el peak
        if str_sitios_encontrados == 'SitiosEncontrados': 
            str_estado_peak = 'GenConSitio'; 
            # Recorro cada gen cercano
            for gen_cerca in L_genes_cerca:
                gen_start = gen_cerca.start; 
                gen_end = gen_cerca.end; 
                # Recorro cada sitio de union
                for SU_encontrado in L_pos_sitios_union:
                    pos_ini_su = SU_encontrado[0]; 
                    pos_end_su = SU_encontrado[1]; 
                    # Mido cercania por distancia a gen_start
                    dist = min(abs(gen_start-pos_ini_su), abs(gen_start-pos_end_su)); 
                    # Veo si dist es menor a dist_max
                    if dist < dist_max:
                        L_su_cerca_genes.append(gen_cerca.gene_id); 
            # Reviso si hay algun gen cerca de un sitio de union
            if len(L_su_cerca_genes) > 0:
                str_estado_peak = str_estado_peak + 'Cercano'; 
            else:
                str_estado_peak = str_estado_peak + 'Lejos'; 
        # Si no hay sitios de union en el peak
        else:
            str_estado_peak = 'GenSinSitio'; 
    # Hay genes cerca del peak pero fueron filtrados (no son protein_coding)
    elif str_genes_cerca == 'GenesFiltrados': 
        # Si hay sitios de union en el peak
        if str_sitios_encontrados == 'SitiosEncontrados': 
            str_estado_peak = 'FiltradoConSitio'; 
            # Recorro cada gen filtrado
            for gen_cerca_filtrado in L_genes_filtrados:
                gen_start = gen_cerca_filtrado.start; 
                gen_end = gen_cerca_filtrado.end; 
                # Recorro cada sitio de union
                for SU_encontrado in L_pos_sitios_union:
                    pos_ini_su = SU_encontrado[0]; 
                    pos_end_su = SU_encontrado[1]; 
                    # Mido cercania por distancia a gen_start
                    dist = min(abs(gen_start-pos_ini_su), abs(gen_start-pos_end_su)); 
                    # Veo si dist es menor a dist_max
                    if dist < dist_max:
                        L_su_cerca_genes.append(gen_cerca_filtrado.gene_id); 
            # Reviso si hay algun gen cerca de un sitio de union
            if len(L_su_cerca_genes) > 0:
                str_estado_peak = str_estado_peak + 'Cercano'; 
            else:
                str_estado_peak = str_estado_peak + 'Lejos'; 
        # Si no hay sitios de union en el peak
        else:
            str_estado_peak = 'FiltradoSinSitio'; 
    # Sin genes cerca
    else: 
        # Con sitios encontrados
        if str_sitios_encontrados == 'SitiosEncontrados':
            str_estado_peak = 'SitiosSinGen'; 
        # Sin sitios encontrados
        else:
            str_estado_peak = 'SinSitiosNiGen'; 
    # Defino str para guardar genes con sitios cerca
    str_genes_sitios_cerca = ''; 
    # Recorro L_su_cerca_genes
    for su_cerca_gen in L_su_cerca_genes:
        str_genes_sitios_cerca = str_genes_sitios_cerca + str(su_cerca_gen) + sep_sub; 
    str_genes_sitios_cerca = str_genes_sitios_cerca.rstrip(sep_sub); '''
    # Agrego los elementos correspondientes a busqueda de sitios de union cerca de genes a L_out
    L_out.append(str_estado_peak); 
    L_out.append(str_genes_sitios_cerca); 
    return L_out


def _main():
    # Funcion para probar funciones en ejecucion del archivo

    # Cargo un archivo .bed
    M_bed_dupays = abrir_bed(bed_dupays, dir_arch=path_bed_main); 
    # Genero una matriz bed de prueba con cantidad a definir de peaks
    M_bed_test = M_bed_dupays[:]; 
    #shuffle(M_bed_test); 
    #M_bed_test = M_bed_test[:20]; 
    # Uso procesar_bed() 
    M_out = procesar_bed(M_bed_test, nom_genoma_dupays, genoma_dupays, dist_main, L_sitios_main, L_col=[0,1,2], path_fasta=path_fasta_main); 
    # Defino header de procesar_bed
    L_header = ['chr_n', 'contig', 'pos_ini', 'pos_end', 'sitios', 'seq_sitios', 'pos_sitios', 'genes', 'id_genes', 'class', 'id_genes_su_cerca', 'peak_id', 'peak_score']; 
    # Guardo M_out en archivo .csv
    M_out = guardar_matriz(M_out, 'DupaysClassPeaks_SitiosConf_dist1500', dir_out=path_out_main, L_head=L_header); 
    return M_out


#################################### RESULTADOS ###################################

output_dump = []; 

if __name__=='__main__':
    output_dump.append(_main()); 

