
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
# X Agregar info de genes cercanos a peaks y sitios de union
# X Abrir resultados RNAseq
# X Seleccionar genes RNAseq
    # * Ver en 9-ClasificacionPeaks3.py, 8-TraduccionRatonHumano.py, 8-ConteoGenesDupaysRNAseq.py y 7-ClasificacionPeaks2.py
    # X Extraer genes RNAseq
    # X Pasarlos a ID de Ensembl
    # X Seleccionar genes en M_genes por lista RNAseq y agregar columna
    # X Agregar info de updown a M_genes
# Generar secuencias para MEME-ChIP
    # Ver 4-GeneracionFastaMEME.py en Scripts_0_16
# Agregar cosas relacionadas a estudiar varios factores de transcripcion
###
'''


#################################### VARIABLES ####################################

# Direcciones de archivos y outputs para ib3 y casa
path_git_casa = 'C:\\Users\\Admin\\Documents\\Scripts\\Git\\ChIP-seq-analysis-scripts\\v1.0_scripts_simplificados\\'; 
path_git_ib3 = 'C:\\Users\\emili\\Documents\\Archivos nube\\Git\\ChIPseqanalysisscripts\\v1.0_scripts_simplificados\\'; 
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
nom_pssm_nkx25_human = 'NKX25_HUMAN.H11MO.0.B.pcm'; 
nom_pssm_nkx25_mouse = 'NKX25_MOUSE.H11MO.0.A.pcm'; 
score_mult = 0.9; # Multiplicador del score maximo de pssm que se usa como cutoff
organism = 'human'; # human o mouse
test_used = 0; # 0 para correr todo completo, numeros mayores a 0 para correr subsets de largo test_used



#################################### FUNCIONES ####################################


def _main():
    '''Funcion para probar scripts en ejecucion del archivo'''

    if organism.lower()=='human' or organism.lower()=='humano':
        pssm_usado = abrir_pssm(nom_pssm_nkx25_human, path_arch=path_pwm_human); 
        genoma_usado = hg19; 
        nom_genoma_usado = 'hg19'; 
        nom_bed_usado = 'Anderson2018'; 
        nom_rnaseq_usado = 'Anderson2018_RNAseq_source'; 
        id_col_usado = 1; 
        updown_col_usado = 3; # Columna log(fold change)
        l_translate = [0, hg19]; 
        nom_output = 'anderson_full'; 
    elif organism.lower()=='mouse' or organism.lower()=='raton':
        pssm_usado = abrir_pssm(nom_pssm_nkx25_mouse, path_arch=path_pwm_mouse); 
        genoma_usado = mm9; 
        nom_genoma_usado = 'mm9'; 
        nom_bed_usado = 'Dupays2015'; 
        nom_rnaseq_usado = 'Dupays2015_RNAseq_source'; 
        id_col_usado = 12; 
        updown_col_usado = 2; # Columna log(fold change)
        l_translate = []; 
        nom_output = 'dupays_full'; 
    else:
        print('Organismo ' + str(organism) + ' no reconocido.')
        return ''
    score_cutoff_pssm = pssm_usado.max*score_mult; 
    
    #M_peaks, M_su, M_genes = pipeline_generador(nom_bed_usado, nom_rnaseq_usado, nom_output, genoma_usado, nom_genoma_usado, pssm_usado, score_cutoff_pssm, 
    #                                            path_bed=path_git_main, path_out=path_out_main, path_fasta=path_fasta_main, dist_max_gen=dist_max_main, 
    #                                            L_su=L_confirmados, test_mode=test_used, id_col_rnaseq=id_col_usado, updown_col_rnaseq=updown_col_usado, l_translate=l_translate); 

    ### FALTA
    # pipeline_meme_chip()
    # Scripts para otros TF
    ###

    return ''


def pipeline_generador(nom_bed, nom_rnaseq, nom_out_base, genoma_ensembl, nombre_genoma, pssm, score_cutoff, path_bed='', path_out='', path_fasta='', dist_max_gen=1000000, L_su=[], 
                       id_col_rnaseq=0, updown_col_rnaseq=1, l_translate=[], test_mode=0):
    '''Pipeline desde abrir archivos .bed y resultados RNA-seq hasta generar tablas de sitios y genes que se usan por otros scripts.
    genoma_ensembl es el elemento genoma de ensembl, usado para buscar genes cerca de peaks en el archivo .bed.
    dist_max_gen es la distancia maxima a la que se buscan genes cercanos.
    L_su es una lista de secuencias consideradas sitios de union.
    pssm es una matriz de puntaje para secuencias (FORMATO UNIFICADO POR UNA FUNCION)
    score_cutoff es el score minimo aceptado para una secuencia correcta de pssm (uso 90% del puntaje maximo posible).
    l_translate se usa para traducir nombres de genes extraidos de archivo nom_rnaseq. Si esta vacio, no se traduce. Sino, tiene que tener columna de nombre de gen y genoma usado para confirmar Ensembl ID.
    test_mode define si se hace una prueba con una seleccion de peaks (se hace si es un numero distinto de 0).'''

    # Abro archivo .bed
    M_peaks = abrir_bed(nom_bed, path_arch=path_bed); 
    # Abro los resultados de RNA-seq
    dict_rnaseq = abrir_rnaseq(nom_rnaseq, path_arch=path_bed, sep=';', ext='.csv', id_col=id_col_rnaseq, updown_col=updown_col_rnaseq, traducir=l_translate); 
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
    # Defino largo de M_peaks
    len_peaks=len(M_peaks); 
    # Recorro M_peaks para trabajar sobre cada peak
    for i in range(len_peaks):
        curr_peak = M_peaks[i]; 
        # Reviso que chr_n este en keys de dict_chr_n
        if curr_peak[0] in dict_chr_n.keys():
            # Uso genes_cercanos_peak() para buscar genes cercanos
            M_genes_cerca = genes_cercanos_peak(curr_peak[0], int(curr_peak[1]), int(curr_peak[2]), dist_max_gen, genoma_ensembl); 
            # Agrego info de genes up/down regulados a M_genes_cerca y extraigo info para M_peaks y M_su
            M_genes_cerca, n_genes_total, n_genes_up, n_genes_down = agregar_info_rnaseq(M_genes_cerca, dict_rnaseq); 
            # Uso secuencia_peak() para obtener la secuencia
            seq_peak = secuencia_peak(dict_chr_n[curr_peak[0]], int(curr_peak[1]), int(curr_peak[2])); 
            # Uso sitios_union_lista() para buscar sitios de union por lista
            M_su_lista = sitios_union_lista(seq_peak, curr_peak[0], L_su, int(curr_peak[1])-1); 
            # Uso sitios_union_pssm() para buscar sitios de union por pssm
            M_su_pssm = sitios_union_pssm(seq_peak, curr_peak[0], pssm, score_cutoff, int(curr_peak[1])-1); 
            # Recorro M_su_lista y M_su_pssm para agregar info de genes up/down regulados
            for j in range(len(M_su_lista)):
                M_su_lista[j].append(''); # Agrego una columna vacia porque pssm incluye score
                M_su_lista[j].append(len(M_genes_cerca)); 
                M_su_lista[j].append(int(n_genes_total)); 
                M_su_lista[j].append(int(n_genes_up)); 
                M_su_lista[j].append(int(n_genes_down)); 
            for j in range(len(M_su_pssm)):
                M_su_pssm[j].append(len(M_genes_cerca)); 
                M_su_pssm[j].append(int(n_genes_total)); 
                M_su_pssm[j].append(int(n_genes_up)); 
                M_su_pssm[j].append(int(n_genes_down)); 
            # Agrego info de genes cercanos y sitios de union a M_peaks
            M_peaks[i].append(str(len(M_genes_cerca))); 
            M_peaks[i].append(str(len(M_su_lista)+len(M_su_pssm))); 
            M_peaks[i].append(str(len(M_su_lista))); 
            M_peaks[i].append(str(len(M_su_pssm))); 
            M_peaks[i].append(len(M_genes_cerca)); 
            # Agrego info de genes up/down regulados a M_peaks
            M_peaks[i].append(int(n_genes_total)); 
            M_peaks[i].append(int(n_genes_up)); 
            M_peaks[i].append(int(n_genes_down)); 
            # Agrego M_genes_cerca a M_genes
            M_genes = M_genes + M_genes_cerca; 
            # Agrego M_su_lista a M_su
            M_su = M_su + M_su_lista; 
            # Agrego M_su_pssm a M_su
            M_su = M_su + M_su_pssm
        ### Display
        if ((i+1)%100==0) or i==0:
            print('Avance: ' + str(i+1) + ' de ' + str(len_peaks))
        ###
    # Elimino genes duplicados
    M_genes = eliminar_duplicados(M_genes); 
    # Defino l_head para las distintas matrices
    l_head_su = ['chr_n', 'pos_ini', 'pos_end', 'seq', 'source', 'score_pssm', 'n_genes_cerca', 'n_updown_total', 'n_upreg', 'n_downreg']; 
    l_head_genes = ['gene_id', 'chr_n', 'pos0', 'biotype', 'fold_change']; 
    l_head_peaks = ['chr_n', 'pos_ini', 'pos_end', 'n_genes', 'n_su', 'n_su_lista', 'n_su_pssm', 'n_genes_cerca', 'n_updown_total', 'n_upreg', 'n_downreg']; 
    # Guardo las matrices con guardar_matriz()
    M_su = guardar_matriz(nom_out_base+'_sitios_union', M_su, path_out=path_out, l_head=l_head_su); 
    M_genes = guardar_matriz(nom_out_base+'_genes', M_genes, path_out=path_out, l_head=l_head_genes); 
    M_peaks = guardar_matriz(nom_out_base+'_peaks', M_peaks, path_out=path_out, l_head=l_head_peaks); 
    return M_peaks, M_su, M_genes


def pipeline_meme_chip():
    '''Genera un archivo .fasta con las secuencias de picos de ChIP-seq para ser mandados a MEME-ChIP.
    Se seleccionan picos de acuerdo a distintos criterios.'''
    ### FALTA
    # Abrir archivos de output_git (peaks o sitios de union?)
    # Seleccionar sitios/peaks
    # Unificar largos (ver en scripts)
    # Obtener secuencias (ver en scripts)
    # Guardar en archivo .fasta (ver en scripts)
    ###
    pass


### Funciones principales del pipeline


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
        # Solo agrego los protein_coding (Puesto aca para que sea facil modificarlo)
        if curr_gen.biotype=='protein_coding':
            # Cargo L_gen a M_out
            M_out.append(L_gen[:]); 
    return M_out


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
                curr_su = [chr_n, pos_ini, pos_ini+len_pssm-1, str(seq_encontrada), 'pssm', float(score)]; 
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
            curr_su = [str(chr_n), pos_ini_ref+i+1, pos_ini_ref+i+len(seq_union), str(reverso_union), "lista"]; 
        # Reviso si curr_su esta registrado
        if len(curr_su) > 0:
            # Si se encontro un sitio de union, curr_su tiene largo mayor a 0 y lo registro en L_su
            L_su.append(curr_su[:]); 
    return L_su


### Funciones simples de procesamiento


def agregar_info_rnaseq(M_genes, dict_rnaseq):
    '''Revisar que genes de M_genes estan en dict_rnaseq.
    dict_rnaseq es un diccionario de genes con ensembl_id como keys y float(fold_change) como valores
    M_genes es una lista de genes en formato [ensembl_id, chr_n, pos0, biotype]
    Se agrega, al final de cada gen, una columna con fold_change (0 si no esta en dict_rnaseq)
    Tambien devuelve n_genes_total (up- y downregulados), n_genes_up, n_genes_down'''

    # Inicializo valores que se devuelven
    n_genes_total = 0; 
    n_genes_up = 0; 
    n_genes_down = 0; 
    # Inicializo matriz que se devuelve
    M_out = []; 
    # Inicializo una lista de genes RNA-seq en base
    # Recorro M_genes
    for i in range(len(M_genes)):
        curr_gen = M_genes[i]; 
        # Defino curr_gen_name
        curr_gen_name = curr_gen[0]; 
        # Inicializo la lista que se agrega a M_out
        l_out = curr_gen[:]; 
        # Veo si curr_gen_name esta en dict_rnaseq.keys()
        if curr_gen_name in dict_rnaseq.keys():
            updown = dict_rnaseq[curr_gen_name]; 
            # Si esta, agrego info de fold_change a l_out
            l_out.append(updown*1); 
            # Agrego 1 a n_genes_total
            n_genes_total+=1; 
            # Veo si hay que agregar +1 a up o down
            if updown>0:
                n_genes_up+=1; 
            elif updown<0:
                n_genes_down+=1; 
            else:
                print('WARNING: Gen ' + str(curr_gen_name) + ' con fold change=0')
        else:
            # Si el gen no esta en dict_rnaseq, agrego string vacio a l_out
            l_out.append(''); 
        M_out.append(l_out[:]); 
    return M_out, n_genes_total, n_genes_up, n_genes_down


def chrn_to_contig(chr_n):
    '''Funcion para pasar de chr_n a contig.
    Recibe un string de formato "chr"+N, donde N puede ser un numero o letras (X, Y, M).
    Devuelve un string. con solo N'''
    return chr_n[3:]


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


def eliminar_duplicados(l_in, col_index=0):
    '''Funcion que elimina todos los elementos duplicados de una lista.
    Busca por col_index y asume que l_in es una lista de listas.'''

    # Inicializo la lista que se devuelve
    l_out = []; 
    l_index = []; 
    len_l_in = len(l_in); 
    # Recorro l_in
    for i in range(len_l_in):
        curr_elem = l_in[i]; 
        # Veo si curr_elem esta en l_out
        if not (curr_elem[col_index] in l_index):
            # Agrego curr_elem a l_out y curr_elem[col_index] a l_index
            l_out.append(l_in[i]); 
            l_index.append(curr_elem[col_index]); 
        ### Display
        if (i+1)%1000==0 or i==0:
            print('Progreso: ' + str(i+1) + ' de ' + str(len_l_in))
        ###
    return l_out


def secuencia_peak(record_seq, pos_ini, pos_end):
    '''Funcion para obtener la secuencia del rango definido por chr_n, pos_ini y pos_end.
    Devuelve un string con una secuencia de ADN.'''

    # Extraigo la secuencia de record_seq
    seq_out = record_seq[pos_ini-1:pos_end]; 
    return seq_out


### Funciones para abrir y guardar archivos


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


def abrir_pssm(nom_arch, path_arch='', solo_pssm=True, pseudocounts=0.5):
    '''Funcion para abrir una matriz de pesos para un sitio de union.'''

    # Defino dir_arch
    if path_arch == '':
        dir_arch = str(nom_arch); 
    else:
        dir_arch = os.path.join(path_arch, nom_arch); 
    # Abro pcm con motifs.parse
    handle = open(dir_arch); 
    record = motifs.parse(handle, 'pfm-four-columns'); 
    handle.close(); 
    # Veo que record tenga un solo motif
    if len(record) == 1:
        m = record[0]; 
    # Si tiene mas de un motif, tiro error, hago print a todos los motif y agarro el primero
    elif len(record) > 1:
        print('WARNING: ' + nom_arch + ' tiene ' + str(len(record)) + ' motifs. Se devuelve solo el primero.')
        m = record[0]; 
        print('Motifs en record:')
        for i in record:
            print(i.name)
            print(i.counts)
    # Si no hay ningun motif, tiro error y dejo m como string vacio
    else:
        print('ERROR: No se encontro ningun motif en ' + nom_arch + '. Se devuelve string vacio.')
        m = ''; 
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


def abrir_rnaseq(nom_arch, path_arch='', sep='\t', ext='.csv', id_col=0, updown_col=1, traducir=[], header=True):
    '''Funcion que abre archivos con output de estudios de RNA-seq y devuelve listas de genes con numeros de up- o downregulacion.
    id_col determina la columna del identificador del gen.
    updown_col determina la columna del numero correspondiente a up- o downregulacion.
        Se asume que contiene log(fold change) y que se obtiene el fold change con signo(updown)*2^(abs(updown))
    traducir se usa para traducir nombres de genes extraidos de archivo nom_rnaseq. Si es una lista o string vacio, no se traduce. 
    Si traducir tiene elementos, tiene que tener el numero de la columna con el nombre del gen y genoma usado para confirmar Ensembl ID (col_gene_name, genoma).'''

    # Inicializo la matriz que se devuelve
    dict_out = {}; 
    # Inicializo booleano skip_header
    skip_header = header; 
    # Defino la direccion del archivo con nom_arch y path_arch
    if path_arch=='':
        filepath = nom_arch + ext; 
    else:
        filepath = os.path.join(path_arch, nom_arch + ext); 
    ### Display
    # Inicializo contador de genes encontrados y genes totales
    cont_encontrados = 0; 
    cont_totales = 0; 
    ###
    # Abro el archivo filepath
    with open(filepath, 'r') as f_out:
        print('Archivo ' + nom_arch + ' abierto.')
        # Recorro cada fila
        for curr_line in f_out.readlines():
            # Veo si salteo header
            if skip_header:
                skip_header=False; 
            else:
                # Transformo la linea en lista
                l_line = curr_line.rstrip().split(sep=sep); 
                # Defino id y updown
                main_id_raw = str(l_line[id_col]); 
                try:
                    fc_abs = fc_from_log_fc(float(l_line[updown_col])); 
                except:
                    print('ERROR con fc_abs')
                    print(l_line)
                # Separo main_id por '|'
                for main_id in main_id_raw.split('|'):
                    # Veo si traducir contiene elementos, si no contiene, guardo ids normalmente
                    if len(traducir)==0:
                        # Veo que main_id tenga algo
                        if main_id != '':
                            # Veo si main_id esta en dict_out.keys()
                            if not (main_id in dict_out.keys()):
                                # Agrego fold change transformado con fc_from_log_fc()
                                dict_out[main_id] = fc_abs*1; 
                            elif not (dict_out[main_id] == fc_abs*1):
                                print('ERROR: id repetido y valores de fc diferentes.')
                                print('fc_abs = ' + str(fc_abs))
                                print('dict_out[' + str(main_id) + '] = ' + str(dict_out[main_id]))
                                if len(main_id_raw.split('|'))==1:
                                    dict_out[main_id] = fc_abs*1; 
                            else:
                                print('WARNING: id ' + str(main_id) + ' repetido en dict_out.keys().')
                            # Sumo 1 a gen encontrado si se encontro un ID
                            cont_encontrados+=1; 
                        # Sumo 1 a genes totales aunque no se haya encontrado un ID
                        cont_totales+=1; 
                    # Si contiene elementos, trato de buscar los ensembl IDs correspondientes a los nombres de genes
                    else:
                        # Inicializo un booleano de genes encontrados
                        gen_encontrado = False; 
                        # Defino el ID de entrez y el nombre del gen
                        entrez_id = main_id; # Por ahora no se usa
                        gene_name = l_line[traducir[0]]; 
                        genoma = traducir[1]; 
                        # Uso .gene_ids_of_gene_name() para ver si el gen aparece en el genoma
                        try:
                            l_genes = genoma.gene_ids_of_gene_name(gene_name); 
                            # Veo si hay genes en l_genes
                            for curr_gen in l_genes:
                                gen_encontrado = True; 
                                # Veo si curr_gen esta en dict_out.keys()
                                if not (curr_gen in dict_out.keys()):
                                    # Agrego fold change transformado con fc_from_log_fc()
                                    dict_out[curr_gen] = fc_abs*1; 
                                elif not (dict_out[curr_gen] == fc_abs*1):
                                    print('ERROR: id repetido y valores de fc diferentes.')
                                    print('fc_abs = ' + str(fc_abs))
                                    print('dict_out[' + str(curr_gen) + '] = ' + str(dict_out[curr_gen]))
                                else:
                                    print('WARNING: id ' + str(curr_gen) + ' repetido en dict_out.keys().')
                        except:
                            # Tiro error si no se encontro nada
                            print('No se encontro ningun gen con gene_name ' + str(gene_name))
                        # Sumo 1 a gen encontrado si se encontro un ID
                        if gen_encontrado:
                            cont_encontrados+=1; 
                        # Sumo 1 a genes totales aunque no se haya encontrado un ID
                        cont_totales+=1; 
    ### Display
    print('Genes totales: ' + str(cont_totales))
    print('Genes traducidos: ' + str(cont_encontrados))
    print('Genes encontrados: ' + str(len(dict_out.keys())))
    ###
    return dict_out


def fc_from_log_fc(fc_float, base=2):
    '''Funcion para devolver el fold change a partir de log(fold change).
    Se asume que se obtiene el fold change con signo(fc_float)*base^(abs(fc_float))'''

    # Defino el exponente como abs(fc_float)
    exponent_fc = abs(fc_float); 
    # Defino fc_abs exponenciando la base y multiplicando por el signo de fc_float
    fc_abs = sign(fc_float)*(base**exponent_fc); 
    return fc_abs


def guardar_matriz(nom_out, M_out, path_out='', ext='.csv', sep=';', l_head=[]):
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
        # Si l_head tiene elementos, los uso de titulos para la tabla
        if len(l_head)>0:
            # Inicializo el string que se escribe a f_out
            row_str = ''; 
            # Recorro l_head
            for k in range(len(l_head)):
                # Agrego str(l_head[k]) a row_str
                row_str = row_str + str(l_head[k]) + sep; 
            # Una vez terminado, elimino la ultima ocurrencia de sep de row_str
            row_str = row_str.rstrip(sep); 
            # Agrego end of line y guardo en f_out
            f_out.write(row_str + '\n'); 
        # Recorro M_out
        for i in range(len(M_out)):
            curr_row = M_out[i]; 
            # Inicializo el string que se escribe a f_out
            row_str = ''; 
            # Recorro curr_row
            for j in range(len(curr_row)):
                # Agrego str(curr_row[j]) a row_str
                row_str = row_str + str(curr_row[j]) + sep; 
            # Una vez terminado, elimino la ultima ocurrencia de sep de row_str
            row_str = row_str.rstrip(sep); 
            # Agrego end of line y guardo en f_out
            f_out.write(row_str + '\n'); 
    return M_out


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
    for i in range(len(L_chr_n)):
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


def sign(n, zero=0, nan=0):
    '''Funcion signo, hecha a mano porque python no la tiene.
    zero es el valor que se devuelve si n==0.
    nan es el valor que se devuelve si n no es un numero.'''
    # Inicializo el valor que se devuelve
    ret = ''; 
    # Uso try para que se rompa si n no es un numero
    try:
        # Veo si n es igual a 0
        if n==0:
            ret = zero; 
        else:
            # Si n no es igual a cero, uso esta funcion para evitar if/else
            ret = (int(n>0)*2)-1; 
    except:
        # Si n no es un numero, devuelvo nan
        ret = nan; 
    return ret


#################################### RESULTADOS ###################################

output_dump = []; 

if __name__=='__main__':
    output_dump.append(_main()); 

