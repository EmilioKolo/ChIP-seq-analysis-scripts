
ib3 = False; 
http_proxy_ib3 = 'http://proxy.fcen.uba.ar:8080'; 
https_proxy_ib3 = 'https://proxy.fcen.uba.ar:8080'; 
ftp_proxy_ib3 = 'ftp://proxy.fcen.uba.ar:8080'; 

# Generales
import os 
import time 
import math 
import sys 
import pandas as pd 
import numpy as np 

# Para debugging
from random import shuffle 

# Especificos bioinfo
from pyensembl import EnsemblRelease 
from Bio import Entrez, SeqIO, motifs, Seq 
Entrez.email = 'ekolomenski@gmail.com'; 
Entrez.api_key = '9bccbd4949f226a5930868df36c211e8b408'; 
import biomart 

# Importo analisischip
path_analisischip_casa = 'C:\\Users\\Admin\\Documents\\Scripts\\Git\\AnalisisSecuenciasChIP\\'; 
path_analisischip_ib3 = 'C:\\Users\\emili\\Documents\\Archivos nube\\Git\\analisissecuenciaschip\\'; 
if ib3:
    path_analisischip_main = path_analisischip_ib3; 
else:
    path_analisischip_main = path_analisischip_casa; 
sys.path.insert(1, path_analisischip_main); 
#from analisischip.seq import seq_data, seq_handler 

# Graficos
import matplotlib.pyplot as plt 
import seaborn as sns 

# Genomas de referencia
hg19 = EnsemblRelease(75, species='human'); 
hg38 = EnsemblRelease(102, species='human'); 
mm9 = EnsemblRelease(54, species='mouse'); 
mm10 = EnsemblRelease(102, species='mouse'); 


'''
### FALTA
# Abrir archivos .bed
# Pasar a formato usado en .csv
    # chr_n, contig, pos_ini, pos_end
    # Sitios de union encontrados (por secuencia y por matriz)
        # Genes cercanos
        # Sitios de union otros TF cercanos
'''


#################################### VARIABLES ####################################

# Direcciones de archivos y outputs para ib3 y casa
path_git_casa = 'C:\\Users\\Admin\\Documents\\Scripts\\Git\\ChIP-seq-analysis-scripts\\v1.0_scripts_simplificados\\'; 
path_git_ib3 = 'C:\\Users\\emili\\Documents\\Archivos nube\\Git\\ChIP-seq-analysis-scripts\\v1.0_scripts_simplificados\\'; 
path_fasta_casa = 'D:\\ArchivosDoctorado\\Genomas\\'; 
path_fasta_ib3 = 'X:\\Genomas\\'; 
path_output_dump_casa = 'D:\\ArchivosDoctorado\\Output_dump\\'; 
path_output_dump_ib3 = 'X:\\Output_dump\\'; 

### Variables main()

# Path main depende de si estoy en ib3 o casa
if ib3:
    path_git_main = path_git_ib3; 
    path_fasta_main = path_fasta_ib3; 
    path_output_dump_main = path_output_dump_ib3; 
else:
    path_git_main = path_git_casa; 
    path_fasta_main = path_fasta_casa; 
    path_output_dump_main = path_output_dump_casa; 
path_in_main = path_git_main; ### Cambiar si es necesario
path_out_main = path_output_dump_main + 'output_git\\'; 
# Path que dependen de la especie
path_pwm_human = path_git_main + 'PWM_human\\'; 
path_pwm_mouse = path_git_main + 'PWM_mouse\\'; 

# Variables para pipeline
dist_max_main = 1000000; 
L_confirmados = ['GCAAGTG', 'GGAAGTG', 'GAAAGTG', 'ATAAGTG', 'GTAAGTG', 'CTAAGTG', 'TCAAGTG', 'TGAAGTG', 'TAAAGTG', 'TTAAGTG']; 



#################################### FUNCIONES ####################################


def _main():
    '''Funcion para probar scripts en ejecucion del archivo'''

    pssm_usado = ''; 
    score_max_pssm = 0; 
    genoma_usado = hg19; 
    nom_genoma_usado = 'hg19'; 
    M_peaks, M_su, M_genes = pipeline_generador('Anderson2018', 'anderson_test', genoma_usado, nom_genoma_usado, pssm_usado, score_max_pssm, 
                                                path_bed=path_git_main, path_out=path_out_main, path_fasta=path_fasta_main, dist_max_gen=dist_max_main, 
                                                L_su=L_confirmados, test_mode=100); 
    ### FALTA
    # Armar pruebas de pipeline_generador
    ###
    return ''


def pipeline_generador(nom_bed, nom_out_base, genoma_ensembl, nombre_genoma, pssm, score_max, path_bed='', path_out='', path_fasta='', dist_max_gen=1000000, L_su=[], test_mode=0):
    '''Pipeline desde abrir archivos .bed y resultados RNA-seq hasta generar tablas de sitios y genes que se usan por otros scripts.
    genoma_ensembl es el elemento genoma de ensembl, usado para buscar genes cerca de peaks en el archivo .bed.
    dist_max_gen es la distancia maxima a la que se buscan genes cercanos.
    L_su es una lista de secuencias consideradas sitios de union.
    pssm es una matriz de puntaje para secuencias (FORMATO UNIFICADO POR UNA FUNCION)
    score_max es el score maximo aceptado para una secuencia correcta de pssm (uso 90% del puntaje maximo posible).
    test_mode define si se hace una prueba con una seleccion de peaks (se hace si es un numero distinto de 0)'''

    # Abro archivo .bed
    M_peaks = abrir_bed(nom_bed, path_arch=path_bed); 
    # Obtengo la secuencia de todos los archivos fasta del genoma con elementos SeqIO
    dict_chr_n = seqio_chr_n(M_peaks, nombre_genoma, path_fasta); 
    ### Reviso si hago test
    if test_mode>0:
        shuffle(M_peaks); 
        M_peaks = M_peaks[:test_mode]; 
    ###
    # Inicializo la matriz de sitios de union
    M_su = []; 
    # Inicializo la matriz de genes
    M_genes = []; 
    # Recorro M_peaks para trabajar sobre cada peak
    for i in range(len(M_peaks)):
        curr_peak = M_peaks[i]; 
        # Reviso que chr_n este en keys de dict_chr_n
        if curr_peak[0] in dict_chr_n.keys():
            # Uso genes_cercanos_peak() para buscar genes cercanos
            M_genes = M_genes + genes_cercanos_peak(curr_peak[0], curr_peak[1], curr_peak[2], dist_max_gen, genoma_ensembl); 
            # Uso secuencia_peak() para obtener la secuencia
            seq_peak = secuencia_peak(dict_chr_n[curr_peak[0]], curr_peak[1], curr_peak[2]); 
            # Uso sitios_union_lista() para buscar sitios de union por lista
            M_su = M_su + sitios_union_lista(seq_peak, curr_peak[0], L_su, int(curr_peak[1])-1); 
            # Uso sitios_union_pssm() para buscar sitios de union por pssm
            M_su = M_su + sitios_union_pssm(seq_peak, curr_peak[0], pssm, score_max, int(curr_peak[1])-1); 
            ### FALTA
            # Abrir un archivo e ir guardando seq_peak (para no tener que volver a buscarlas)
            ###
    ### FALTA
    # Recorrer M_su para chequeos (capaz es innecesario)
    # Guardar L_peaks?? (M_peaks con secuencias?)
    ###
    M_su = guardar_matriz(nom_out_base+'_sitios_union', M_su, path_out=path_out); 
    M_genes = guardar_matriz(nom_out_base+'_genes', M_genes, path_out=path_out); 
    return M_peaks, M_su, M_genes


def genes_cercanos_peak(chr_n, pos_ini, pos_end, dist_max, genoma_ensembl):
    '''Funcion para buscar los genes cercanos al rango definido por chr_n, pos_ini y pos_end en genoma_ensembl.
    Devuelve listas de genes, donde cada gen es una lista con ID, chr_n, pos0, tipo (prot_coding/otro) (ver si agrego algo mas).'''

    # Inicializo M_out
    M_out = []; 
    # Defino contig
    contig = chrn_to_contig(chr_n); 
    # Busco genes alrededor del peak con genoma_ensembl.genes_at_locus()
    L_genes_cerca_raw = genoma_ensembl.genes_at_locus(contig, pos_ini-dist_max, pos_end+dist_max); 
    # Recorro L_genes_cerca_raw
    for i in range(len(L_genes_cerca_raw)):
        curr_gen = L_genes_cerca_raw[i]; 
        # Inicializo L_gen, que se agrega a M_out
        L_gen = []; 
        # Defino pos0 con definir_pos0()
        pos0 = definir_pos0(curr_gen.start, curr_gen.end, curr_gen.strand); 
        # Cargo datos relevantes a L_gen
        L_gen.append(curr_gen.gene_id); # ID del gen para volver a buscarlo
        L_gen.append(chr_n); # chr_n para tenerlo en formato igual que el resto de scripts
        L_gen.append(pos0); # pos0 determina el +1 del gen
        L_gen.append(curr_gen.biotype); # Me importan mas los protein_coding que el resto
        # Cargo L_gen a M_out
        M_out.append(L_gen[:]); 
    return M_out


def secuencia_peak(record_seq, pos_ini, pos_end):
    '''Funcion para obtener la secuencia del rango definido por chr_n, pos_ini y pos_end.
    nombre_genoma esta pensado para funcionar con hg19, hg38, mm9 y mm10. 
    Devuelve un string con una secuencia de ADN.'''

    # Extraigo la secuencia de record_seq
    seq_out = record_seq[pos_ini-1:pos_end]; 
    return seq_out


def sitios_union_lista(seq_peak, chr_n, L_sitios, pos_ini_ref):
    '''Funcion para encontrar todas las ocurrencias de cada una de las secuencias en L_sitios dentro de seq_peak.
    Devuelve listas de sitios de union, donde cada sitio incluye chr_n, pos_ini, pos_end, seq y fuente="lista" (ver si agrego algo mas).'''

    # Inicializo la lista de sitios de union que se devuelve
    L_su = []; 
    # Recorro L_sitios
    for i in range(len(L_sitios)):
        curr_su = L_sitios[i]; 
        # Agrego la lista de sitios de union de sitio_union_str() para cada sitio de union en L_sitios
        L_su = L_su + sitio_union_str(seq_peak, curr_su, chr_n, pos_ini_ref=pos_ini_ref); 
    return L_su


def sitios_union_pssm(seq_peak, chr_n, pssm, score_cutoff, pos_ini_ref):
    '''Funcion para encontrar todas las secuencias con puntaje mayor a score_max dentro de seq_peak.
    Devuelve listas de sitios de union, donde cada sitio incluye chr_n, pos_ini, pos_end, seq, fuente="pssm" y score (ver si agrego algo mas).'''

    # Inicializo la lista de sitios de union que se devuelve
    L_su = []; 
    ### FALTA
    # Ver como se hace en seq.py (funcion _buscar_PSSM_en_seq())
    ###
    # Defino el largo de la secuencia buscada
    len_pssm = len(pssm.consensus); 
    # Veo que el largo de seq_referencia sea mayor o igual al largo de pssm.consensus
    if len(seq_peak) >= len_pssm:
        # Uso try por si se agarra otro error
        try:
            # Uso pssm.search() para obtener una lista de posiciones con scores mayores a score_cutoff
            for position, score in pssm.search(seq_peak, threshold=score_cutoff):
                # Reinicio curr_su
                curr_su = []; 
                # Defino la secuencia encontrada
                seq_encontrada = seq_peak[position:position+len_pssm]; 
                # Defino pos_out en base a position y pos_ini_ref
                if position < 0:
                    pos_ini = pos_ini_ref + len(seq_peak) + position; 
                    forward = False; 
                else:
                    pos_ini = position+pos_ini_ref; 
                    forward = True; 
                # Agrego chr_n, pos_ini y pos_end a curr_su
                curr_su = [chr_n, pos_ini, pos_ini+len_pssm-1, str(seq_encontrada), 'pssm', int(score)]; 
                # Agrego curr_su a L_su
                L_su.append(curr_su[:]); 
        except:
            print('ERROR buscando PSSM en seq ' + str(seq_peak))
    return L_su


def sitio_union_str(seq_referencia, seq_union, chr_n, pos_ini_ref=0):
    '''Funcion para buscar una secuencia de adn (seq_union) en otra secuencia mas larga (seq_referencia). 
    Devuelve una lista con cada una de las posiciones donde pega seq_union o su reverso complementario.
    Cada sitio de union se registra con chr_n, pos_ini, pos_end, seq y fuente="lista".'''

    # Inicializo la lista de sitios de union que se devuelve
    L_su = []; 
    # Inicializo el reverso complementario de seq_union
    reverso_union = complemento_secuencia(seq_union).upper(); 
    # Recorro seq_referencia
    for i in range(len(seq_referencia)-len(seq_union)+1):
        # Inicializo curr_su
        curr_su = []; 
        # Defino la secuencia que se revisa
        curr_seq = str(seq_referencia[i:i+len(seq_union)]).upper(); 
        # Si seq_union es igual a curr_seq, registro pos_ini, pos_end, seq_union
        if str(seq_union).upper() == str(curr_seq).upper():
            curr_su = [str(chr_n), pos_ini_ref+i+1, pos_ini_ref+i+len(seq_union), str(seq_union), "lista"]; 
        # Si reverso_union es igual a curr_seq, registro pos_ini, pos_end, reverso_union
        elif reverso_union == str(curr_seq).upper():
            curr_su = [pos_ini_ref+i+1, pos_ini_ref+i+len(seq_union), str(reverso_union)]; 
        # Reviso si curr_su esta registrado
        if len(curr_su) > 0:
            # Si se encontro un sitio de union, curr_su tiene largo mayor a 0 y lo registro en L_su
            L_su.append(curr_su[:]); 
    return L_su


def guardar_matriz(nom_out, M_out, path_out='', ext='.csv', sep=';'):
    '''Funcion para guardar una matriz en un archivo .csv.'''

    # Defino la direccion del archivo con nom_out y path_out
    if path_out=='':
        filepath = nom_out + ext; 
    else:
        filepath = os.path.join(path_out, nom_out + ext); 
    # Creo el archivo filepath
    with open(filepath, 'w') as f_out:
        print('Archivo ' + nom_out + ext + ' creado.')
    # Vuelvo a abrirlo en modo append para escribirlo
    with open(filepath, 'a') as f_out:
        # Recorro M_out
        for i in range(len(M_out)):
            curr_row = M_out[i]; 
            # Inicializo el string que se escribe a f_out
            row_str = ''; 
            # Recorro curr_row
            for j in range(len(curr_row)):
                # Agrego str(curr_row[j]) a row_str
                row_str = str(curr_row[j]) + sep; 
            # Una vez terminado, elimino la ultima ocurrencia de sep de row_str
            row_str = row_str.rstrip(sep); 
            # Agrego end of line y guardo en f_out
            f_out.write(row_str + '\n'); 
    return M_out


def complemento_secuencia(seq):
    '''Funcion que recibe una secuencia y devuelve el reverso complementario.'''

    # Defino el diccionario de traduccion de ADN
    dict_adn = {'T':'A','U':'A','A':'T','C':'G','G':'C','Y':'R','R':'Y','K':'M','M':'K',
                'W':'W','S':'S','D':'H','H':'D','V':'B','B':'V','N':'N'}; 
    # Inicializo la variable que se devuelve
    ret_seq = ''; 
    # Recorro seq de atras para adelante
    for i in range(len(seq)):
        # Defino el nucleotido n, siendo el de la posicion -i-1 (python cuenta para atras con indice en 1)
        n = seq[-i-1]; 
        # Reviso que el nucleotido a traducir este en las keys de los diccionarios (ambos tienen mismas keys)
        if not (n in dict_adn.keys()):
            # Si N no esta entre las keys, devuelvo nucleotido 'N'
            print('Nucleotido "' + str(n) + '" no interpretado. Se devuelve N.'); 
            ret = 'N'; 
        else:
            # Si esta entre las keys, uso dict_adn para complementar n
            ret = dict_adn[n]; 
        # Agrego ret a ret_seq
        ret_seq = ret_seq + ret; 
    return ret_seq


def chrn_to_contig(chr_n):
    '''Funcion para pasar de chr_n a contig.
    Recibe un string de formato "chr"+N, donde N puede ser un numero o letras (X, Y, M).
    Devuelve un string. con solo N'''
    return chr_n[3:]


def seqio_chr_n(peaks, nombre_genoma, path_fasta=''):
    '''Funcion para extraer todos los archivos .fasta de contigs que aparezcan en la matriz de peaks.'''

    # Inicializo el diccionario que se devuelve
    dict_out = {}; 
    # Inicializo una lista de chr_n unicos en matriz peaks
    L_chr_n = []; 
    # Recorro peaks
    for i in range(len(peaks)):
        curr_peak = peaks[i]; 
        # Defino chr_n como curr_peak[0]
        curr_chr_n = curr_peak[0]; 
        # Si curr_chr_n no esta en L_chr_n, lo agrego
        if not (str(curr_chr_n) in L_chr_n):
            L_chr_n.append(str(curr_chr_n)); 
    # Recorro L_chr_n
    for i in L_chr_n:
        chr_n = L_chr_n[i]; 
        # Defino el nombre y direccion del archivo a buscar
        nom_fasta = nombre_genoma + '_' + chr_n + '.fasta'; 
        # Defino fulldir_fasta con nom_fasta y path_fasta
        if path_fasta=='':
            fulldir_fasta = nom_fasta; 
        else:
            fulldir_fasta = os.path.join(path_fasta, nom_fasta); 
        # Pruebo abrir el archivo
        try:
            # Abro el fasta con SeqIO
            record = SeqIO.read(fulldir_fasta, 'fasta'); 
            # Agrego record.seq a dict_out
            dict_out[chr_n] = record.seq; # Asumo que los chr_n en L_chr_n no estan repetidos
        except:
            print('ERROR: chr_n ' + str(chr_n) + ' no se pudo abrir. nom_fasta: ' + str(nom_fasta) + ' ; path_fasta: ' + str(path_fasta))
    return dict_out


def definir_pos0(start, end, strand):
    '''Funcion que recibe el start, end y strand de un elemento gen de EnsemblRelease y devuelve el +1.'''

    # Veo si strand es + o -
    if strand == '+':
        ret = start; 
    elif strand == '-':
        ret = end; 
    # Si no es ninguno de los dos, devuelvo gene_start por defecto
    else:
        print('WARNING: Strand ' + str(strand) + ' no reconocido. Se devuelve start como +1.')
        ret = start; 
    return ret


def abrir_bed(nom_arch, path_arch='', sep='\t', ext='.bed', col_num=3):
    '''Funcion que abre archivos .bed y devuelve una matriz con una lista de peaks.
    col_num determina el numero de columnas que se devuelven. Por defecto, solo recopilo chr_n, pos_ini y pos_end.'''
    
    # Inicializo la matriz que se devuelve
    M_out = []; 
    # Defino la direccion del archivo con nom_arch y path_arch
    if path_arch=='':
        filepath = nom_arch + ext; 
    else:
        filepath = os.path.join(path_arch, nom_arch + ext); 
    # Abro el archivo filepath
    with open(filepath, 'r') as f_out:
        print('Archivo ' + nom_arch + ' abierto.')
        # Recorro cada fila
        for curr_line in f_out.readlines():
            # Transformo la linea en lista
            L_line = curr_line.rstrip().split(sep=sep); 
            # Cargo L_line en M_out
            M_out.append(L_line[:col_num]); 
    return M_out


#################################### RESULTADOS ###################################

output_dump = []; 

if __name__=='__main__':
    output_dump.append(_main()); 

