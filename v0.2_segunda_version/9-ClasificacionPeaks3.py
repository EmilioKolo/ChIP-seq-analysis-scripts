
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
VERSION ESPECIFICA PARA ANDERSON CON RESULTADOS RNASEQ TRADUCIDOS

Funciones para clasificar los distintos peaks de resultados de ChIP-seq y generar archivos .fasta para MEME-ChIP
Funcion principal: pipeline_bed()
Incluye funciones para: 
    Buscar genes cercanos a picos en archivo .bed
    Buscar sitios de union por lista de sitios confirmados o por matriz PSSM
    Confirmar sitios de union por genes en RNA-seq
'''


#################################### VARIABLES ####################################

# Direcciones de archivos y outputs
path_dropbox_casa = 'D:\\Users\\Admin\\Dropbox\\'; 
path_dropbox_ib3 = 'C:\\Users\\emili\\Dropbox\\'; 
path_fasta_casa = 'D:\\Archivos doctorado\\Genomas\\'; 
path_fasta_ib3 = 'X:\\Genomas\\'; 
path_output_dump_casa = 'D:\\Archivos doctorado\\Output_dump\\'; 
path_output_dump_ib3 = 'X:\\Output_dump\\'; 

'''path_bed_casa = 'D:\\Users\\Admin\\Dropbox\\Doctorado\\3-Genes transactivados rio abajo\\0-Fuentes\\Papers ChIP-seq\\'; 
path_bed_ib3 = 'C:\\Users\\emili\\Dropbox\\Doctorado\\3-Genes transactivados rio abajo\\0-Fuentes\\Papers ChIP-seq\\'; 
path_fasta_casa = 'D:\\Archivos doctorado\\Genomas\\'; 
path_fasta_ib3 = 'X:\\Genomas\\'; 
path_out_casa = 'D:\\Archivos doctorado\\Output_dump\\PeaksClasificados\\'; 
path_out_ib3 = 'X:\\Output_dump\\PeaksClasificados\\'; 
path_rnaseq_casa = 'D:\\Users\\Admin\\Dropbox\\Doctorado\\3-Genes transactivados rio abajo\\1-Interpretacion de datos ChIP-seq\\'; 
path_rnaseq_ib3 = 'C:\\Users\\emili\\Dropbox\\Doctorado\\3-Genes transactivados rio abajo\\1-Interpretacion de datos ChIP-seq\\'; 
path_pwm_casa = 'D:\\Users\\Admin\\Dropbox\\Doctorado\\3-Genes transactivados rio abajo\\1-Interpretacion de datos ChIP-seq\\PWM\\'; 
path_pwm_ib3 = 'C:\\Users\\emili\\Dropbox\\Doctorado\\3-Genes transactivados rio abajo\\1-Interpretacion de datos ChIP-seq\\PWM\\'; 
'''

# Nombres de archivos .bed y datos relacionados
bed_dupays = 'Dupays2015'; 
bed_anderson = 'Anderson2018-GSE89457consensus'; 
bed_test_dupays = 'Dupays2015test'; 
L_bed = [bed_dupays, bed_anderson]; 
nom_genoma_dupays = 'mm9'; 
genoma_dupays = mm9; 
nom_genoma_anderson = 'hg19'; 
genoma_anderson = hg19; 
nom_anderson = 'Anderson'; 
nom_dupays = 'Dupays'; 

nombre_csv_rnaseq_humano = 'Anderson2018_RNAseq'; # .csv
nombre_csv_rnaseq_mouse = 'Dupays2015_RNAseq'; # .csv
nombre_rnaseq_mouse = '9-lista_genes_rnaseq_dupays'; # .txt
nombre_rnaseq_human = '9-lista_genes_rnaseq_anderson'; # .txt

# Datos para busquedas
#L_dist = [100*1000, 1000*1000, 50*1000, 10*1000, 1500]; 
L_dist = [10*1000, 100*1000, 1000*1000]; 
L_aagtg = ['AAGTG']; 
L_confirmados = ['GCAAGTG', 'GGAAGTG', 'GAAAGTG', 'ATAAGTG', 'GTAAGTG', 'CTAAGTG', 'TCAAGTG', 'TGAAGTG', 'TAAAGTG', 'TTAAGTG']; 
#L_arch_hocomoco_pcm_mouse = ['NKX25_MOUSE.H11MO.0.A.pcm', 'TBX20_MOUSE.H11MO.0.C.pcm', 'MEIS1_MOUSE.H11MO.1.A.pcm', 'TGIF1_MOUSE.H11MO.0.A.pcm', 
#                             'HAND1_MOUSE.H11MO.0.C.pcm', 'MAF_MOUSE.H11MO.1.A.pcm', 'GATA1_MOUSE.H11MO.1.A.pcm', 'GATA6_MOUSE.H11MO.0.A.pcm', 'GATA4_MOUSE.H11MO.0.A.pcm']; 
#L_nombres_pssm_mouse = ['Nkx25', 'Tbx20', 'Meis1', 'Tgif1', 'Hand1', 'Maf', 'Gata1', 'Gata6', 'Gata4']; 
#L_arch_hocomoco_pcm_mouse = ['NKX25_MOUSE.H11MO.0.A.pcm']; 
#L_nombres_pssm_mouse = ['Nkx25']; 
L_arch_hocomoco_pcm_mouse = []; 
L_nombres_pssm_mouse = []; 
#L_arch_hocomoco_pcm_human = ['NKX25_HUMAN.H11MO.0.B.pcm', 'TBX20_HUMAN.H11MO.0.D.pcm', 'MEIS1_HUMAN.H11MO.1.B.pcm', 'TGIF1_HUMAN.H11MO.0.A.pcm', 'HAND1_HUMAN.H11MO.0.D.pcm', 
#                             'MAF_HUMAN.H11MO.1.B.pcm', 'GATA1_HUMAN.H11MO.1.A.pcm', 'GATA6_HUMAN.H11MO.0.A.pcm', 'GATA4_HUMAN.H11MO.0.A.pcm']; 
#L_nombres_pssm_human = ['NKX25', 'TBX20', 'MEIS1', 'TGIF1', 'HAND1', 'MAF', 'GATA1', 'GATA6', 'GATA4']; 
#L_arch_hocomoco_pcm_human = ['NKX25_HUMAN.H11MO.0.B.pcm']; 
#L_nombres_pssm_human = ['NKX25']; 
L_arch_hocomoco_pcm_human = []; 
L_nombres_pssm_human = []; 
#L_arch_hocomoco_pcm = []; 
#L_nombres_pssm_main = []; 
'''Equivalencias archivos pcm
NKX25_MOUSE.H11MO.0.A.pcm	NKX25_HUMAN.H11MO.0.B
TBX20_MOUSE.H11MO.0.C.pcm	TBX20_HUMAN.H11MO.0.D
MEIS1_MOUSE.H11MO.1.A.pcm	MEIS1_HUMAN.H11MO.1.B
TGIF1_MOUSE.H11MO.0.A.pcm	TGIF1_HUMAN.H11MO.0.A
HAND1_MOUSE.H11MO.0.C.pcm	HAND1_HUMAN.H11MO.0.D
 MAF_MOUSE.H11MO.1.A.pcm	MAF_HUMAN.H11MO.1.B
GATA1_MOUSE.H11MO.1.A.pcm	GATA1_HUMAN.H11MO.1.A
GATA6_MOUSE.H11MO.0.A.pcm	GATA6_HUMAN.H11MO.0.A
GATA4_MOUSE.H11MO.0.A.pcm	GATA4_HUMAN.H11MO.0.A
'''

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
path_pwm_main = path_dropbox_main + 'Doctorado\\3-Genes transactivados rio abajo\\1-Interpretacion de datos ChIP-seq\\PWM_human\\'; 
path_out_main = path_output_dump_main + '\\PeaksClasificados\\'; 

# Defino si trabajo con anderson o dupays
nom_genoma_usado = 'human'; 
test_state = False; 
modificador_nomout = '_rnaseq_notrad_pssm'; 

if nom_genoma_usado.lower() in ['human', 'humano']:
    nom_main = nom_anderson; 
    bed_main = bed_anderson;  
    nom_genoma_main = nom_genoma_anderson; 
    genoma_main = genoma_anderson; 
    L_nombres_pssm_main = L_nombres_pssm_human; 
    L_arch_hocomoco_pcm = L_arch_hocomoco_pcm_human; 
    path_pwm_main = path_dropbox_main + 'Doctorado\\3-Genes transactivados rio abajo\\1-Interpretacion de datos ChIP-seq\\PWM_human\\'; 
    nombre_rnaseq_main = nombre_rnaseq_human; 
elif nom_genoma_usado.lower() in ['mouse', 'raton']:
    nom_main = nom_dupays; 
    bed_main = bed_dupays;  
    nom_genoma_main = nom_genoma_dupays; 
    genoma_main = genoma_dupays; 
    L_nombres_pssm_main = L_nombres_pssm_mouse; 
    L_arch_hocomoco_pcm = L_arch_hocomoco_pcm_mouse; 
    path_pwm_main = path_dropbox_main + 'Doctorado\\3-Genes transactivados rio abajo\\1-Interpretacion de datos ChIP-seq\\PWM_mouse\\'; 
    nombre_rnaseq_main = nombre_rnaseq_mouse; 

dist_main = L_dist[0]; 
L_dist_main = L_dist; 
#L_dist_main = [1000000]; 
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


def abrir_L_pssm(L_nom_arch, path_arch='', pseudocounts=0.5):
    # Extrae una lista de matrices de counts y devuelve sus PSSM correspondientes

    # Inicializo la lista que se devuelve
    L_pssm = []; 
    # Recorro L_nom_arch
    for nom_arch in L_nom_arch:
        # Uso abrir_pssm() con solo_pssm=True para extraer los PSSM
        L_pssm.append(abrir_pssm(nom_arch, path_arch=path_arch, solo_pssm=True, pseudocounts=pseudocounts)); 
    return L_pssm


def abrir_pssm(nom_arch, path_arch='', solo_pssm=True, pseudocounts=0.5):
    # Extrae la matriz de counts de nom_arch y la pasa a PSSM

    # Extraigo el motif con extraer_motif_pcm()
    m = extraer_motif_pcm(nom_arch, path_arch=path_arch); 
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


def abrir_txt(nom_arch, path_arch='', ext='.txt'):
    # Abre un archivo txt y devuelve cada linea como un string en una lista

    # Inicializo la lista que se devuelve
    L_out = []; 
    # Defino la direccion del archivo con nom_arch y dir_arch
    if path_arch=='':
        filepath = nom_arch + ext; 
    else:
        filepath = os.path.join(path_arch, nom_arch + ext); 
    # Abro el archivo en modo reading
    with open(filepath, 'r') as F_txt:
        print('Archivo ' + nom_arch + ' abierto.')
        # Recorro cada fila
        for curr_line in F_txt.readlines():
            # Cargo curr_line sin end of line en L_out
            L_out.append(curr_line.rstrip()); 
    return L_out


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
            #print()
            #print('Peak (con dist_max agregada): ' + str(chr_n) + ', ' + str(pos_ini-dist_max) + ', ' + str(pos_end+dist_max))
            #print('End cerca de peak: ' + str(curr_gen))
            pass
        elif peak_dentro_de_gen:
            #print()
            #print('Peak (con dist_max agregada): ' + str(chr_n) + ', ' + str(pos_ini-dist_max) + ', ' + str(pos_end+dist_max))
            #print('Peak dentro de gen, end/start lejos de peak: ' + str(curr_gen))
            pass
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


def buscar_genes_rnaseq(L_genes_cerca, L_genes_rnaseq, str_genes_sitios_cerca, sep_sub=','):
    # Inicializo una lista de genes rnaseq en el peak y genes rnaseq con sitios de union cerca
    L_genes_rnaseq_cerca_su = []; 
    L_genes_rnaseq_cerca = []; 
    # Rearmo la lista de IDs de union y de genes cerca
    L_genes_su_cerca = str_genes_sitios_cerca.split(sep_sub); 
    # Recorro L_genes_su_cerca
    for gen_su_cerca in L_genes_su_cerca:
        # Si alguno de los genes con sitio de union cerca esta en L_genes_rnaseq
        if gen_su_cerca in L_genes_rnaseq:
            # Reviso no haberlo agregado repetido ni que sea gen vacio
            if (not (gen_su_cerca in L_genes_rnaseq_cerca_su)) and (gen_su_cerca!=''):
                # Lo agrego a L_genes_rnaseq_cerca_su
                L_genes_rnaseq_cerca_su.append(str(gen_su_cerca)); 
    # Recorro L_genes_cerca
    for gen_cerca in L_genes_cerca:
        curr_gen = gen_cerca.gene_id; 
        # Si alguno de los genes cerca de peaks esta en L_genes_rnaseq
        if curr_gen in L_genes_rnaseq:
            # Reviso no haberlo agregado repetido ni que sea gen vacio
            if (not (curr_gen in L_genes_rnaseq_cerca)) and (curr_gen!=''):
                # Lo agrego a L_genes_rnaseq_cerca
                L_genes_rnaseq_cerca.append(str(curr_gen)); 
    # Defino string con genes rnaseq con sitios de union cerca y con genes rnaseq en peak
    str_genes_rnaseq_cerca_su = ''; 
    str_genes_rnaseq_cerca = ''; 
    # Agrego cada uno de los genes en rnaseq con sitio de union cerca
    for gen_rnaseq_cerca_su in L_genes_rnaseq_cerca_su:
        if gen_rnaseq_cerca_su != '':
            str_genes_rnaseq_cerca_su = str_genes_rnaseq_cerca_su + str(gen_rnaseq_cerca_su) + sep_sub; 
        else:
            print('gen_rnaseq_cerca vacio en ' + str(L_genes_rnaseq_cerca_su))
    while len(str_genes_rnaseq_cerca_su)>0 and str_genes_rnaseq_cerca_su[-1]==sep_sub:
        str_genes_rnaseq_cerca_su = str_genes_rnaseq_cerca_su.rstrip(sep_sub); 
    # Agrego cada uno de los genes en rnaseq cerca de peak
    for gen_rnaseq_cerca in L_genes_rnaseq_cerca:
        if gen_rnaseq_cerca != '':
            str_genes_rnaseq_cerca = str_genes_rnaseq_cerca + str(gen_rnaseq_cerca) + sep_sub; 
        else:
            print('gen_rnaseq_cerca vacio en ' + str(L_genes_rnaseq_cerca))
    while len(str_genes_rnaseq_cerca)>0 and str_genes_rnaseq_cerca[-1]==sep_sub:
        str_genes_rnaseq_cerca = str_genes_rnaseq_cerca.rstrip(sep_sub); 
    # Devuelvo str_genes_rnaseq_cerca y str_genes_rnaseq_cerca_su para agregar a L_out
    return str_genes_rnaseq_cerca, str_genes_rnaseq_cerca_su


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


def buscar_sitios_cerca_pssm(dist_max, L_genes_cerca, L_genes_filtrados, L_su_pssm, sep_sub=','):
    # Busca sitios de union (encontrados con pssm) cerca de genes que esten cerca de un peak

    # Inicializo las variables que se devuelven
    str_genes_su_cerca = ''; 
    str_su_genes_cerca = ''; 
    str_scores_su_cerca = ''; 
    # Inicializo variables de chequeo
    score_cerca_max = 0; 
    L_genes_su_cerca = []; 
    # Recorro L_su_pssm (si no tiene ningun elemento la funcion no hace nada)
    for su_pssm in L_su_pssm: 
        # Defino posiciones iniciales/finales y score del sitio
        pos_ini_su = su_pssm[0]; 
        pos_end_su = su_pssm[1]; 
        score_pssm = su_pssm[2]; 
        # Recorro los genes cercanos no filtrados
        for gen_cerca in L_genes_cerca:
            # Defino posiciones iniciales/finales del gen
            gen_start = gen_cerca.start; 
            gen_end = gen_cerca.end; 
            # Mido cercania por distancia a gen_start
            dist = min(abs(gen_start-pos_ini_su), abs(gen_start-pos_end_su)); 
            # Veo si dist es menor a dist_max
            if dist < dist_max:
                # Reviso que el gene_id no se haya agregado antes
                if not (str(gen_cerca.gene_id) in L_genes_su_cerca):
                    # Agrego el gen a str_genes_su_cerca
                    str_genes_su_cerca = str_genes_su_cerca + str(gen_cerca.gene_id) + sep_sub; 
                    # Agrego el gen a L_genes_su_cerca asi no se agrega de nuevo
                    L_genes_su_cerca.append(str(gen_cerca.gene_id)); 
                # Agrego el sitio a str_su_genes_cerca y el score a str_scores_su_cerca
                str_su_genes_cerca = str_su_genes_cerca + str(pos_ini_su) + '_' + str(pos_end_su) + sep_sub; 
                str_scores_su_cerca = str_scores_su_cerca + str(score_pssm) + sep_sub; 
                # Reviso para ver el score maximo
                if score_pssm > score_cerca_max:
                    score_cerca_max = score_pssm+0; 
    # Agrego score_cerca_max al principio de str_scores_su_cerca
    if score_cerca_max > 0:
        str_scores_su_cerca = str(score_cerca_max) + sep_sub + str_scores_su_cerca; 
    # Elimino el ultimo sep_sub de los strings
    str_genes_su_cerca = str_genes_su_cerca.rstrip(sep_sub); 
    str_su_genes_cerca = str_su_genes_cerca.rstrip(sep_sub); 
    str_scores_su_cerca = str_scores_su_cerca.rstrip(sep_sub); 
    return str_genes_su_cerca, str_su_genes_cerca, str_scores_su_cerca


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


def buscar_sitios_peak_pssm(chr_n, pos_ini, pos_end, genome_name, genome, pssm, path_fasta='', sep_sub=',', score_cutoff=5.0):
    # Busca sitios de union en un peak usando PSSM

    # Inicializo el elemento seq_data que consigue las secuencias y hace el analisis de pssm
    sequence_data = seq_data(genome_name, genome_element=genome, path_fasta=path_fasta); 
    # Consigo la secuencia del peak con seq_data._consulta_secuencia_fasta()
    seq_peak = sequence_data._consulta_secuencia_fasta(chr_n, pos_ini, pos_end); 
    # Busco todos los sitios de union que hayan para pssm
    L_sitios_pssm = sequence_data._buscar_PSSM_en_seq(pssm, seq_peak, score_cutoff=score_cutoff, pos_ini_ref=pos_ini-1); # pos_ini-1 por como suma _buscar_PSSM_en_seq()
    # Inicializo las variables que se devuelven
    str_sitios_pssm = ''; 
    str_scores_pssm = ''; 
    L_su_pssm = []; 
    score_max = 0; 
    # Defino el largo de pssm
    len_pssm = len(pssm.consensus); 
    # Recorro L_sitios_pssm
    for sitio_pssm in L_sitios_pssm:
        # Defino posicion inicial/final y score de sitio_pssm
        pos_ini_sitio = sitio_pssm[0]+1; 
        pos_end_sitio = sitio_pssm[0]+len_pssm; 
        score_pssm = sitio_pssm[1]; 
        # Reviso para ver el score maximo
        if score_pssm > score_max:
            score_max = score_pssm+0; 
        # Agrego los datos a str_sitios_pssm y str_scores_pssm
        str_sitios_pssm = str_sitios_pssm + str(pos_ini_sitio) + '_' + str(pos_end_sitio) + sep_sub; 
        str_scores_pssm = str_scores_pssm + str(score_pssm) + sep_sub; 
        # Agrego datos del sitio a L_su_pssm
        L_su_pssm.append([pos_ini_sitio, pos_end_sitio, score_pssm]); 
    # Agrego score_max al principio de str_scores_pssm
    if score_max > 0:
        str_scores_pssm = str(score_max) + sep_sub + str_scores_pssm; 
        # Elimino el ultimo sep_sub de los strings
        str_sitios_pssm = str_sitios_pssm.rstrip(sep_sub); 
        str_scores_pssm = str_scores_pssm.rstrip(sep_sub); 
    return str_sitios_pssm, str_scores_pssm, L_su_pssm


def cargar_dict(nom_dict, path_dict='.\\', ext='.csv', sep=';'):
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


def extraer_motif_pcm(nom_arch, path_arch=''):
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
    #attributes_general = ['mgi_symbol', 'ensembl_gene_id', 'entrezgene_accession', 'entrezgene_id']; 
    attributes_general = ['ensembl_gene_id', 'entrezgene_accession', 'entrezgene_id']; 
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
        #gene_symbol = line[0]; 
        ensembl_gene_id = line[0]; 
        entrez_accession = line[1]; 
        entrez_gene_id = line[2]; 
        # Reviso que entrez_accession o entrez_gene_id no sean valores vacios y los cargo al dict
        if len(entrez_gene_id):
            entrez_to_ensembl[entrez_gene_id] = ensembl_gene_id; 
        if len(entrez_accession):
            entrez_to_ensembl[entrez_accession] = ensembl_gene_id; 
    return entrez_to_ensembl


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
            str_head=str_head.rstrip(sep); 
            # Agrego end of line y guardo en F_out
            F_out.write(str_head + '\n'); 
        # Recorro M_out
        for L_out in M_out:
            curr_str = ''; 
            # Recorro L_out
            for i in L_out:
                curr_str = curr_str + str(i) + sep; 
            # Elimino la ultima ocurrencia de sep
            curr_str=curr_str.rstrip(sep); 
            # Agrego end of line y guardo en F_out
            F_out.write(curr_str + '\n'); 
    return M_out


def guardar_dict(dict_out, nom_out, path_out='.\\', ext='.csv', sep=';'):
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


def pipeline_bed_anderson(nom_bed, nom_genoma, genoma, dist_max, L_sitios=[], L_pssm=[], L_nombres_pssm=[], nom_out='', nom_rnaseq='', 
                 L_col=[0,1,2], path_bed='', path_fasta='', path_out='.\\', path_rnaseq='', extra_header=[]):
    # Corre abrir_bed(), procesar_bed() y guardar_matriz() en orden

    # Cargo los peaks en M_bed
    M_bed = abrir_bed(nom_bed, dir_arch=path_bed); 
    # Defino los genes upregulados en RNA-seq si hay nom_rnaseq
    if nom_rnaseq != '':
        L_genes_rnaseq = abrir_txt(nom_rnaseq, path_arch=path_rnaseq); 
    else:
        L_genes_rnaseq = []; 
    # Proceso M_bed
    M_out, L_headers_out = procesar_bed(M_bed, nom_genoma, genoma, dist_max, L_sitios=L_sitios, L_pssm=L_pssm, L_nombres_pssm=L_nombres_pssm, 
                                        L_col=L_col, L_genes_rnaseq=L_genes_rnaseq, path_fasta=path_fasta); 
    # Solo guardo si nom_out tiene texto
    if nom_out != '':
        M_out = guardar_matriz(M_out, nom_out, dir_out=path_out, L_head=L_headers_out + extra_header); 
    return M_out


def procesar_bed(M_bed, genome_name, genome, dist_max, L_sitios=[], L_pssm=[], L_nombres_pssm=[], L_col=[0,1,2], L_genes_rnaseq=[], path_fasta='', display_rate=20):
    # Funcion central que recorre M_bed, procesa los peaks y devuelve un output para guardar en .csv

    ### Display
    len_bed = len(M_bed); 
    cont_display = 0; 
    ###
    # Inicializo la matriz que se devuelve
    M_out = []; 
    # Inicializo la lista de headers
    L_headers_out = []; 
    # Defino id_chr, id_ini, id_end con L_col
    id_chr, id_ini, id_end = L_col[0], L_col[1], L_col[2]; 
    # Recorro la matriz extraida con abrir_bed()
    for peak in M_bed:
        # Cargo chr_n, pos_ini y pos_end con los datos de L_col
        chr_n = peak[id_chr]; 
        pos_ini = min(int(peak[id_ini]), int(peak[id_end])); 
        pos_end = max(int(peak[id_ini]), int(peak[id_end])); 
        # Uso procesar_peak() para hacer todo lo necesario con la info del peak
        L_out, curr_headers = procesar_peak(chr_n, pos_ini, pos_end, genome_name, genome, dist_max, L_sitios=L_sitios, L_pssm=L_pssm, 
                                            L_nombres_pssm=L_nombres_pssm, L_genes_rnaseq=L_genes_rnaseq, path_fasta=path_fasta); 
        # Reviso que L_headers_out no este vacia
        if L_headers_out==[]:
            # Inicializo L_headers out si esta vacia
            L_headers_out = curr_headers[:]; 
        elif L_headers_out!=curr_headers:
            print('WARNING: Error con headers de peak ' + str(peak) + ': ' + str(curr_headers) + ' diferentes de ' + str(L_headers_out))
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
    return M_out, L_headers_out


def procesar_peak(chr_n, pos_ini, pos_end, genome_name, genome, dist_max, L_sitios=[], L_pssm=[], L_nombres_pssm=[], L_genes_rnaseq=[], path_fasta=''):
    # Funcion para hacer todo lo necesario con la info del peak
    # Busca sitios de union (L_sitios), genes cercanos (genome, dist_max) y sitios de union cerca de los genes cercanos

    # Defino un separador de subdivision
    sep_sub = ','; 
    # Inicializo la lista que se devuelve
    L_out = []; 
    # Inicializo una lista de headers
    L_headers_out = []; 
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
    L_headers_out = L_headers_out + ['chr_n', 'contig', 'pos_ini', 'pos_end']; 
    ## Busqueda de genes cercanos
    str_genes_cerca, str_genes_cerca_id, L_genes_cerca, L_genes_filtrados = buscar_genes_cerca_peak(chr_n, pos_ini, pos_end, contig, genome, dist_max, sep_sub=sep_sub); 
    # Agrego los elementos correspondientes a busqueda de genes a L_out
    L_out.append(str_genes_cerca); 
    L_out.append(str_genes_cerca_id); 
    L_headers_out = L_headers_out + ['genes', 'id_genes']; 
    ## Busqueda de sitios de union
    # L_sitios
    if len(L_sitios) > 0:
        str_su_enc_lista, str_seq_su_lista, str_pos_su_lista, L_pos_su_lista = buscar_sitios_peak(chr_n, pos_ini, pos_end, genome_name, genome, L_sitios, 
                                                                                                  path_fasta=path_fasta, sep_sub=sep_sub); 
        # Agrego los elementos correspondientes a busqueda de secuencias a L_out
        L_out.append(str_su_enc_lista); 
        L_out.append(str_seq_su_lista); 
        L_out.append(str_pos_su_lista); 
        L_headers_out = L_headers_out + ['sitios_lista', 'seq_sitios', 'pos_sitios']; 

        ## Busqueda de sitios de union (lista) cerca de genes cercanos
        str_estado_peak, str_genes_sitios_cerca = buscar_sitios_cerca_genes_peak(dist_max, str_genes_cerca, str_su_enc_lista, L_genes_cerca, 
                                                                                 L_genes_filtrados, L_pos_su_lista, sep_sub=sep_sub); 
        # Agrego los elementos correspondientes a busqueda de sitios de union cerca de genes a L_out
        L_out.append(str_estado_peak); 
        L_out.append(str_genes_sitios_cerca); 
        L_headers_out = L_headers_out + ['OLD_class', 'id_genes_su_lista_cerca']; 
    # PSSM
    if len(L_pssm) > 0:
        # Chequeo para ver que se encontro NKX2-5
        nkx25_no_encontrado = True; 
        # Recorro L_pssm
        for p in range(len(L_pssm)):
            pssm = L_pssm[p]; 
            nombre_pssm = L_nombres_pssm[p]; 
            str_sitios_pssm, str_scores_pssm, L_su_pssm = buscar_sitios_peak_pssm(chr_n, pos_ini, pos_end, genome_name, genome, pssm, 
                                                                                  path_fasta=path_fasta, sep_sub=sep_sub); 
            # Agrego los elementos correspondientes a L_out
            L_out.append(str_sitios_pssm); 
            L_out.append(str_scores_pssm); 
            L_headers_out = L_headers_out + ['sitios_pssm_' + nombre_pssm, 'scores_pssm_' + nombre_pssm]; 
            # Reviso si nombre_pssm es NKX25
            if nombre_pssm.upper() in ['NKX25', 'NKX2-5', 'NKX2.5']:
                ## Busqueda de sitios de union de NKX2-5 (pssm) cerca de genes cercanos
                str_genes_su_cerca, str_su_genes_cerca, str_scores_su_cerca = buscar_sitios_cerca_pssm(dist_max, L_genes_cerca, L_genes_filtrados, L_su_pssm, sep_sub=sep_sub); 
                # Agrego los elementos correspondientes a L_out
                L_out.append(str_genes_su_cerca); 
                L_out.append(str_su_genes_cerca); 
                L_out.append(str_scores_su_cerca); 
                L_headers_out = L_headers_out + ['id_genes_su_pssm_nkx25_cerca', 'sitios_pssm_gen_cerca', 'scores_genes_su_pssm_nkx25_cerca']; 
                # Actualizo nkx25_no_encontrado una vez que se corrio buscar_sitios_cerca_pssm()
                nkx25_no_encontrado = False; 
        if nkx25_no_encontrado:
            print('ERROR: no se encontro NKX2-5 entre la lista de PSSM.')
    ## Busqueda de genes en lista RNA-seq (si es que hay lista)
    if len(L_genes_rnaseq) > 0:
        str_genes_rnaseq_cerca, str_genes_rnaseq_cerca_su = buscar_genes_rnaseq(L_genes_cerca, L_genes_rnaseq, str_genes_sitios_cerca, sep_sub=sep_sub); 
        # Agrego str_genes_rnaseq_cerca y str_genes_rnaseq_cerca_su a L_out
        L_out.append(str_genes_rnaseq_cerca); 
        L_out.append(str_genes_rnaseq_cerca_su); 
        L_headers_out = L_headers_out + ['id_genes_rnaseq', 'id_genes_rnaseq_cerca_su']; 
    return L_out, L_headers_out


def seleccionar_genes_rnaseq(nom_rnaseq, genoma, dict_entrez_ensembl, path_rnaseq=''):
    # Selecciona la lista de genes que aparezcan en un archivo de resultados de RNA-seq
    # El archivo de resultados tiene que tener headers y ser un .csv

    # Abro el archivo RNA-seq
    M_rnaseq, L_headers_rnaseq = abrir_csv_headers(nom_rnaseq, dir_arch=path_rnaseq); 
    # Agarro los ids de las columnas que me interesan
    if 'anderson' in nom_rnaseq.lower():
        entrez_id_col = L_headers_rnaseq.index('ENTREZ ID'); 
        gene_symbol_col = L_headers_rnaseq.index('Symbol'); 
        ensembl_id_col = ''; 
        gene_id_col = L_headers_rnaseq.index('ENTREZ ID'); 
        aliases_col = ''; 
        log_fc_col = L_headers_rnaseq.index('log Fold Change'); 
    elif 'dupays' in nom_rnaseq.lower():
        entrez_id_col = L_headers_rnaseq.index('Entrez ID'); 
        gene_symbol_col = L_headers_rnaseq.index('Gene Symbol'); 
        ensembl_id_col = L_headers_rnaseq.index('Ensembl ID'); 
        gene_id_col = L_headers_rnaseq.index('Gene ID'); 
        aliases_col = L_headers_rnaseq.index('Aliases'); 
        log_fc_col = L_headers_rnaseq.index('Log FC ([hypo] vs [co])'); 
    else:
        print('WARNING: ' + nom_rnaseq.lower() + ' no contiene dupays ni anderson. Se asume formato dupays')
        entrez_id_col = L_headers_rnaseq.index('Entrez ID'); 
        gene_symbol_col = L_headers_rnaseq.index('Gene Symbol'); 
        ensembl_id_col = L_headers_rnaseq.index('Ensembl ID'); 
        gene_id_col = L_headers_rnaseq.index('Gene ID'); 
        aliases_col = L_headers_rnaseq.index('Aliases'); 
        log_fc_col = L_headers_rnaseq.index('Log FC ([hypo] vs [co])'); 
    # Inicializo la matriz de genes que se devuelve
    M_genes = []; 
    ### Display
    cont_fail = 0; 
    ###
    # Recorro M_rnaseq
    for i in range(len(M_rnaseq)):
        curr_L = M_rnaseq[i]; 
        # Contadores para estatus de busqueda de id
        curr_buscando = True; 
        curr_status = 0; 
        # Extraigo la info de las columnas que me interesan
        curr_gene_id = curr_L[gene_id_col]; 
        curr_gene_symbol = curr_L[gene_symbol_col]; 
        if aliases_col != '':
            curr_aliases = curr_L[aliases_col].split('|'); 
        else:
            curr_aliases = []; 
        curr_entrez_id = curr_L[entrez_id_col]; 
        if ensembl_id_col != '':
            curr_ensembl_id = curr_L[ensembl_id_col].split('|'); 
        else:
            curr_ensembl_id = []; 
        curr_log_fc = curr_L[log_fc_col]; 
        # Busco curr_entrez_id en dict_entrez_ensembl
        if curr_entrez_id in dict_entrez_ensembl.keys():
            id_ensembl = dict_entrez_ensembl[curr_entrez_id]; 
            # Pruebo si el ensembl_id encontrado esta en el genoma
            try:
                genoma.gene_by_id(id_ensembl); 
                M_genes.append([str(curr_entrez_id), str(id_ensembl), curr_log_fc*1]); 
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
                # Recorro L_gene_symbol
                for g in L_gene_symbol:
                    M_genes.append([str(curr_entrez_id), str(g), curr_log_fc*1]); 
                curr_buscando = False; 
        # Se puede hacer una ultima revision 
        if curr_buscando:
            cont_fail = cont_fail + 1; 
            print()
            print('>>>>> No se encontro Ensembl ID para ' + str(curr_L))
            print()
    print('Cantidad de genes buscados: ' + str(len(M_rnaseq)))
    print('Cantidad de genes encontrados: ' + str(len(M_rnaseq)-cont_fail))
    print('Cantidad de genes NO encontrados: ' + str(cont_fail))
    return M_genes


def _main():
    # Funcion para probar funciones en ejecucion del archivo

    # 
    M_out = []; 

    ## Scripts parseo RNAseq
    L_head_rnaseq_out=['EntrezID', 'EnsemblID', 'LogFC']; 
    '''# Script para generar la lista de genes confirmados en RNA-seq de raton
    dict_entrez_ensembl_raton = generar_dict_entrez_ensembl('mm9', biomart_server='www.ensembl.org/biomart'); 
    M_genes_rnaseq_raton = seleccionar_genes_rnaseq(nombre_csv_rnaseq_mouse, mm9, dict_entrez_ensembl_raton, path_rnaseq=path_rnaseq_main); 
    M_genes_rnaseq_raton = guardar_matriz(M_genes_rnaseq_raton, '17-lista_genes_rnaseq_raton', dir_out=path_rnaseq_main, L_head=L_head_rnaseq_out); '''
    '''# Script para generar la lista de genes confirmados en RNA-seq de humano
    dict_entrez_ensembl_humano = generar_dict_entrez_ensembl('hg19', biomart_server='www.ensembl.org/biomart'); 
    M_genes_rnaseq_humano = seleccionar_genes_rnaseq(nombre_csv_rnaseq_humano, hg19, dict_entrez_ensembl_humano, path_rnaseq=path_rnaseq_main); 
    M_genes_rnaseq_humano = guardar_matriz(M_genes_rnaseq_humano, '17-lista_genes_rnaseq_humano', dir_out=path_rnaseq_main, L_head=L_head_rnaseq_out);''' 

    '''## Scripts pipeline_bed
    # Defino L_pssm
    L_pssm = abrir_L_pssm(L_arch_hocomoco_pcm, path_arch=path_pwm_main, pseudocounts=0.5); 
    # Defino variables en comun en todas las corridas
    L_col = [0,1,2]; 
    extra_header = ['peak_id']; 
    # Corro pipeline_bed para varias distancias
    for curr_dist in L_dist_main:
        # Defino variables necesarias para pipeline_bed
        nom_out = nom_main + 'ClassPeaks_SitiosConf_dist' + str(curr_dist) + modificador_nomout + '_prueba'*test_state; 
        # Corro pipeline_bed
        M_out = pipeline_bed_anderson(bed_main, nom_genoma_main, genoma_main, curr_dist, L_sitios=L_sitios_main, L_pssm=L_pssm, nom_out=nom_out, 
                                      nom_rnaseq=nombre_rnaseq_main, L_col=L_col, path_bed=path_bed_main, path_fasta=path_fasta_main, path_out=path_out_main, 
                                      path_rnaseq=path_rnaseq_main, extra_header=extra_header, L_nombres_pssm=L_nombres_pssm_main); 
    '''
    return M_out


def _test_retry_generar_dict(max_retries=100, path_usado='.\\'):
    # Variable success
    falta_dict = True; 
    retries = 0; 
    # Pruebo generar el diccionario
    while falta_dict and retries<max_retries:
        try:
            # Si funciona, salgo del ciclo while
            dict_entrez_ensembl = generar_dict_entrez_ensembl(nom_genoma_main); 
            falta_dict = False; 
        except:
            retries += 1; 
            print('Carga dict_entrez_ensembl fallo. Intentos: ' + str(retries))
    if not falta_dict:
        try:
            dict_entrez_ensembl = guardar_dict(dict_entrez_ensembl, nom_out='dict_entrez_ensembl', path_out=path_usado); 
            dict_cargado = cargar_dict('dict_entrez_ensembl', path_dict=path_usado); 
            print('dict_entrez_ensembl==dict_cargado:' + str(dict_entrez_ensembl==dict_cargado))
            print('dict_entrez_ensembl.keys()==dict_cargado.keys():'+str(dict_entrez_ensembl.keys()==dict_cargado.keys()))
        except:
            print('ERROR: Fallo guardado de dict_entrez_ensembl')
    return dict_entrez_ensembl


def _test_dict_entrez(nom_rnaseq, genoma, nombre_genoma, path_rnaseq=''):
    # Funcion para probar el diccionario de entrez a ensembl

    # Armo dict entrez_ensembl
    dict_entrez_ensembl = generar_dict_entrez_ensembl(nombre_genoma); 
    # Abro el archivo RNA-seq
    M_rnaseq, L_headers_rnaseq = abrir_csv_headers(nom_rnaseq, dir_arch=path_rnaseq); 
    # Agarro los ids de las columnas que me interesan
    gene_id_col = L_headers_rnaseq.index('Gene ID'); 
    gene_symbol_col = L_headers_rnaseq.index('Gene Symbol'); 
    aliases_col = L_headers_rnaseq.index('Aliases'); 
    entrez_id_col = L_headers_rnaseq.index('Entrez ID'); 
    ensembl_id_col = L_headers_rnaseq.index('Ensembl ID'); 
    # Inicializo el archivo con la lista de genes que se devuelven
    L_genes = []; 
    # Inicio contadores
    cont_ok = 0; 
    cont_no_encontrado = 0; 
    cont_no_enc_arreglado = 0; 
    cont_fail = 0; 
    cont_fail_arreglado = 0; 
    cont_error = 0; 
    cont_separados = 0; 
    L_entrez_id_no_arreglados = []; 
    L_entrez_id_separados = ['241197', '100040462', '620695', '629957', '100503611', '207618', '100322896', '100039133', '78748', '100040259', '73410', 
                             '382492', '213980', '209550', '621837', '100042986', '236451', '223650', '271375', '278795', '627788', '100041206', '668382']; 
    # Recorro M_rnaseq
    for i in range(len(M_rnaseq)):
        curr_buscando = True; 
        curr_status = 0; 
        curr_L = M_rnaseq[i]; 
        # Extraigo la info de las columnas que me interesan
        curr_gene_id = curr_L[gene_id_col]; 
        curr_gene_symbol = curr_L[gene_symbol_col]; 
        curr_aliases = curr_L[aliases_col].split('|'); 
        curr_entrez_id = curr_L[entrez_id_col]; 
        curr_ensembl_id = curr_L[ensembl_id_col].split('|'); 
        # Reviso si curr_entrez_id esta en L_entrez_id_separados
        if curr_entrez_id in L_entrez_id_separados:
            curr_buscando = False; 
            r = '### Gen no procesado: '; 
            r = r + 'entrez_id=' + str(curr_entrez_id) + '; '; 
            r = r + 'gene_symbol=' + str(curr_gene_symbol) + '; '; 
            r = r + 'aliases=' + str(curr_aliases); 
            print(r)
            try:
                print(genoma.gene_ids_of_gene_name(curr_gene_symbol))
            except:
                print('No se encontro ' + curr_gene_symbol)
            for j in curr_aliases:
                if j != '':
                    try:
                        print(genoma.gene_ids_of_gene_name(j))
                    except:
                        print('No se encontro ' + j)
            curr_buscando = False; 
            cont_separados += 1; 
        # Busco curr_entrez_id en dict_entrez_ensembl
        elif curr_entrez_id in dict_entrez_ensembl.keys():
            id_ensembl = dict_entrez_ensembl[curr_entrez_id]; 
            # Pruebo si el ensembl_id encontrado esta en el genoma
            try:
                genoma.gene_by_id(id_ensembl); 
                L_genes.append(str(id_ensembl)); 
                cont_ok += 1; 
                curr_buscando = False; 
            except:
                cont_no_encontrado += 1; 
                curr_status = 1; 
        else:
            cont_fail += 1; 
            curr_status = 2; 
        # Paso por buscar_gene_id() si no se encontro en dict_entrez_ensembl
        if curr_buscando: 
            L_gene_symbol = buscar_gene_id(curr_gene_symbol, curr_aliases, curr_ensembl_id, genoma); 
            if len(L_gene_symbol) > 0:
                # Si hay mas de un id, tiro warning
                if len(L_gene_symbol) > 1:
                    print()
                    print('WARNING: mas de un gene symbol asociado a ' + str(curr_gene_symbol))
                    print()
                L_genes = L_genes + L_gene_symbol[:]; 
                if curr_status == 1:
                    cont_no_enc_arreglado += 1; 
                elif curr_status == 2:
                    cont_fail_arreglado += 1; 
                else:
                    print()
                    print('ERROR cont_buscando=True // curr_status=' + str(curr_status) + ' CON gene_symbol')
                    print()
                    cont_error += 1; 
            elif curr_status == 1:
                L_entrez_id_no_arreglados.append(str(curr_entrez_id))
                print('### ID ' + str(id_ensembl) + ' no encontrado en genoma. No arreglado. Entrez id: ' + str(curr_entrez_id))
            elif curr_status == 2:
                L_entrez_id_no_arreglados.append(str(curr_entrez_id))
                print('### ' + str(curr_entrez_id) + ' no encontrado en diccionario. No arreglado.')
            else:
                print()
                print('ERROR cont_buscando=True // curr_status=' + str(curr_status) + ' SIN gene_symbol')
                print()
    print('cont_ok: ' + str(cont_ok))
    print('cont_no_encontrado: ' + str(cont_no_encontrado))
    print('cont_no_encontrado arreglado: ' + str(cont_no_enc_arreglado))
    print('cont_fail: ' + str(cont_fail))
    print('cont_fail arreglado: ' + str(cont_fail_arreglado))
    print('cont_error: ' + str(cont_error))
    print('cont_separados: ' + str(cont_separados))
    return L_genes


def _test_busq_entrez_id_faltantes():
    L_entrez_id = [241197, 100040462, 620695, 629957, 100503611, 207618, 100322896, 100039133, 78748, 100040259, 73410, 
                   382492, 213980, 209550, 621837, 100042986, 236451, 223650, 271375, 278795, 627788, 100041206, 668382]; 
    dict_entrez_faltantes = {'241197':'ENSMUSG00000092572', '100040462':'ENSMUSG00000090272', '620695':'ENSMUSG00000087006', '207618':'ENSMUSG00000092094', 
                             '100322896':'ENSMUSG00000090326', '78748':'ENSMUSG00000098132', '73410':'ENSMUSG00000096655', '213980':'ENSMUSG00000090173', 
                             '209550':'ENSMUSG00000086022', '236451':'ENSMUSG00000091649', '223650':'ENSMUSG00000118671', '271375':'ENSMUSG00000090176', 
                             '278795':'ENSMUSG00000090291', '668382':'ENSMUSG00000084842', '100041206':'PREDICTED', '627788':'PREDICTED', '621837':'PREDICTED', 
                             '100042986':'PREDICTED', '382492':'PREDICTED', '100040259':'PREDICTED', '100039133':'PREDICTED', '629957':'PREDICTED', 
                             '100503611':'PREDICTED'}; 
    for entrez_id in dict_entrez_faltantes.keys():
        ensembl_id = dict_entrez_faltantes[entrez_id]; 
        if ensembl_id != 'PREDICTED':
            try:
                mm9.gene_by_id(ensembl_id)
                print('OK: ' + ensembl_id + ' encontrado en mm9.')
            except:
                try: 
                    mm10.gene_by_id(ensembl_id); 
                    print('UPS: Se encontro ' + ensembl_id + ' solo en mm10.')
                except:
                    print('ERROR: No se encontro ' + ensembl_id + ' en mm9 ni mm10.')



#################################### RESULTADOS ###################################

output_dump = []; 

if __name__=='__main__':
    output_dump.append(_main()); 

