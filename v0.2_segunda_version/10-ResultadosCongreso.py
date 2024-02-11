
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
Funciones para abrir tablas en PeaksClasificados y procesarlas para producir resultados

* Tabla con picos ChIP-seq
	* Genes cercanos (principalmente confirmados por RNA-seq)
	* Clasificados en regulados directamente/indirectamente (segun si tienen sitio de union NKX2-5)
	* No regulados no le damos mucha bola
	* Sitios de union usar 90% del score maximo
'''

#################################### VARIABLES ####################################

# Variables de testeo
test_state = False; 
verbose_main = True; 
display_interval_main = 1000; 

# Valores de corte
score_fraction_main = 0.9; 
pseudocounts_main = 0.5; 

# Direcciones de archivos y outputs
path_dropbox_casa = 'D:\\Users\\Admin\\Dropbox\\'; 
path_dropbox_ib3 = 'C:\\Users\\emili\\Dropbox\\'; 
path_fasta_casa = 'D:\\Archivos doctorado\\Genomas\\'; 
path_fasta_ib3 = 'X:\\Genomas\\'; 
path_output_dump_casa = 'D:\\Archivos doctorado\\Output_dump\\'; 
path_output_dump_ib3 = 'X:\\Output_dump\\'; 

# Nombres de archivos usados
#L_pc_anderson = ['AndersonClassPeaks_SitiosConf_dist1000000_rnaseq_pssm', 'AndersonClassPeaks_SitiosConf_dist100000_rnaseq_pssm', 
#                 'AndersonClassPeaks_SitiosConf_dist50000_rnaseq_pssm']; 
#L_pc_dupays = ['DupaysClassPeaks_SitiosConf_dist1000000_rnaseq_pssm', 'DupaysClassPeaks_SitiosConf_dist100000_rnaseq_pssm', 
#               'DupaysClassPeaks_SitiosConf_dist50000_rnaseq_pssm', 'DupaysClassPeaks_SitiosConf_dist1500_rnaseq_pssm']; 
### NOMBRES PARA CORRER SOLO 10k
L_pc_anderson = ['AndersonClassPeaks_SitiosConf_dist10000_rnaseq_notrad_pssm']; 
L_pc_dupays = ['DupaysClassPeaks_SitiosConf_dist10000_rnaseq_notrad_pssm']; 

# Nombres de genomas
nom_genoma_dupays = 'mm9'; 
genoma_dupays = mm9; 
nom_genoma_anderson = 'hg19'; 
genoma_anderson = hg19; 

# Nombres de archivos de matrices de peso y factores de transcripcion asociados
L_arch_pcm_raton = ['NKX25_MOUSE.H11MO.0.A.pcm', 'TBX20_MOUSE.H11MO.0.C.pcm', 'MEIS1_MOUSE.H11MO.1.A.pcm', 'TGIF1_MOUSE.H11MO.0.A.pcm', 
                    'HAND1_MOUSE.H11MO.0.C.pcm', 'MAF_MOUSE.H11MO.1.A.pcm', 'GATA1_MOUSE.H11MO.1.A.pcm', 'GATA6_MOUSE.H11MO.0.A.pcm', 'GATA4_MOUSE.H11MO.0.A.pcm']; 
L_nombres_pssm_raton = ['Nkx25', 'Tbx20', 'Meis1', 'Tgif1', 'Hand1', 'Maf', 'Gata1', 'Gata6', 'Gata4']; 
L_arch_pcm_humano = ['NKX25_HUMAN.H11MO.0.B.pcm', 'TBX20_HUMAN.H11MO.0.D.pcm', 'MEIS1_HUMAN.H11MO.1.B.pcm', 'TGIF1_HUMAN.H11MO.0.A.pcm', 'HAND1_HUMAN.H11MO.0.D.pcm', 
                     'MAF_HUMAN.H11MO.1.B.pcm', 'GATA1_HUMAN.H11MO.1.A.pcm', 'GATA6_HUMAN.H11MO.0.A.pcm', 'GATA4_HUMAN.H11MO.0.A.pcm']; 
L_nombres_pssm_humano = ['NKX25', 'TBX20', 'MEIS1', 'TGIF1', 'HAND1', 'MAF', 'GATA1', 'GATA6', 'GATA4']; 

# Nombres archivos out
#L_out_anderson = ['Anderson1M', 'Anderson100k', 'Anderson50k']; 
#L_out_dupays = ['Dupays1M', 'Dupays100k', 'Dupays50k', 'Dupays1500']; 
### NOMBRES PARA CORRER SOLO 10k
L_out_anderson = ['Anderson10k']; 
L_out_dupays = ['Dupays10k']; 

### Variables main()

# Direcciones de archivos y outputs
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
path_pwm_human = path_dropbox_main + 'Doctorado\\3-Genes transactivados rio abajo\\1-Interpretacion de datos ChIP-seq\\PWM_human\\'; 
path_pwm_mouse = path_dropbox_main + 'Doctorado\\3-Genes transactivados rio abajo\\1-Interpretacion de datos ChIP-seq\\PWM_mouse\\'; 
path_pc_main = path_output_dump_main + 'PeaksClasificados\\'; 
path_out_main = path_output_dump_main + 'ArchivosCongreso\\'; 

# Nombres de archivos usados y paths definidos por genoma
nom_pc_main = L_pc_dupays[0]; 
path_pwm_main = path_pwm_mouse; 

# Nombres de genomas
nom_genoma_main = nom_genoma_dupays; 
genoma_main = genoma_dupays; 

# Nombres de archivos de matrices de peso y factores de transcripcion asociados
L_nombres_pssm_main = L_nombres_pssm_raton; 
L_arch_pcm_main = L_arch_pcm_raton; 

# Nombre del output
nom_out_main = L_pc_dupays[0]; 

# Nombres de archivos de testeo
#pc_test = 'testClassPeaks'; 
pc_test = L_pc_dupays[0]; 
nom_genoma_test = nom_genoma_dupays; 
genoma_test = genoma_dupays; 
nom_out_test = 'testDupays1M'; 
L_nombres_pssm_test = L_nombres_pssm_raton; 
L_arch_pcm_test = L_arch_pcm_raton; 
path_pwm_test = path_pwm_mouse; 
# Variables de testeo si es test_state
if test_state: 
    nom_pc_main = pc_test; 
    nom_genoma_main = nom_genoma_test; 
    genoma_main = genoma_test; 
    nom_out_main = nom_out_test; 
    L_nombres_pssm_main = L_nombres_pssm_test; 
    L_arch_pcm_main = L_arch_pcm_test; 
    path_pwm_main = path_pwm_test; 


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


def buscar_score_cutoff_pssm(L_nom_arch, path_arch='', pseudocounts=0.5, score_fraction=0.9):
    # Funcion para buscar la lista de los cutoff para los scores de una lista de pssm

    # Defino L_pssm usando abrir_L_pssm()
    L_pssm = abrir_L_pssm(L_nom_arch, path_arch=path_arch, pseudocounts=pseudocounts); 
    # Inicializo la lista que se devuelve
    L_out = []; 
    # Recorro L_pssm
    for p in range(len(L_pssm)):
        curr_pssm = L_pssm[p]; 
        # Agrego curr_pssm.max*score_fraction como threshold, 0.9 es valor por defecto
        pssm_threshold = curr_pssm.max * score_fraction; 
        L_out.append(pssm_threshold*1.0); 
    return L_out


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


def generar_tabla_final(nom_pc, L_scores_pssm, L_nombres_tf, path_pc='', sub_sep=',', verbose=False, display_interval=100):
    # Funcion que abre el archivo con peaks clasificados y lo procesa a una matriz
    # Simplifica la informacion de genes cercanos (menos la de RNA-seq)
    # Selecciona picos con genes regulados (RNA-seq) y analiza sitios de union de NKX2-5
    # Para otros TF se queda con los sitios con score mayor al 90% del score maximo de la matriz

    # Abro nom_pc en path_pc y saco headers (por las dudas) y matriz
    M_peaks_pc = abrir_csv(nom_pc, dir_arch=path_pc); 
    headers_pc = abrir_csv_headers(nom_pc, dir_arch=path_pc); 
    ### Display
    len_M_peaks = len(M_peaks_pc); 
    ###
    # Inicializo la matriz y la lista de headers que se devuelven
    M_out = []; 
    L_header = []; 
    ## Defino informacion de las columnas que uso
    # Identificadores de peak
    chr_n_col = headers_pc.index('chr_n'); 
    contig_col = headers_pc.index('contig'); 
    peak_ini_col = headers_pc.index('pos_ini'); 
    peak_end_col = headers_pc.index('pos_end'); 
    peak_id_col = headers_pc.index('peak_id'); 
    # Genes cerca
    genes_col = headers_pc.index('genes'); 
    id_genes_col = headers_pc.index('id_genes'); 
    genes_su_cerca_col = headers_pc.index('id_genes_su_lista_cerca'); 
    id_genes_rnaseq_col = headers_pc.index('id_genes_rnaseq'); 
    rnaseq_su_cerca_col = headers_pc.index('id_genes_rnaseq_cerca_su'); 
    # Sitios de union NKX2-5
    sitios_col = headers_pc.index('sitios_lista'); 
    seq_sitios_col = headers_pc.index('seq_sitios'); 
    pos_sitios_col = headers_pc.index('pos_sitios'); 
    id_genes_pssm_nkx25_cerca_col = headers_pc.index('id_genes_su_pssm_nkx25_cerca'); 
    sitios_pssm_nkx25_cerca_col = headers_pc.index('sitios_pssm_gen_cerca'); 
    scores_pssm_nkx25_cerca_col = headers_pc.index('scores_genes_su_pssm_nkx25_cerca'); 
    # Recorro M_peaks_pc
    for i in range(len(M_peaks_pc)):
        curr_peak_pc = M_peaks_pc[i]; 
        ### Display
        if verbose and (i+1)%display_interval==0:
            print('### Progreso: ' + str(i+1) + '/' + str(len_M_peaks))
        # Inicializo la lista que se agrega a M_out
        L_out = []; 
        # Agrego los valores posicionales
        L_out.append(str(curr_peak_pc[chr_n_col])); 
        L_out.append(str(curr_peak_pc[contig_col])); 
        L_out.append(str(curr_peak_pc[peak_ini_col])); 
        L_out.append(str(curr_peak_pc[peak_end_col])); 
        # Solo hago esto la primera vuelta
        if i == 0: 
            L_header = L_header + ['chr_n', 'contig', 'peak_ini', 'peak_end']; 
        # Extraigo la informacion relevante a genes cercanos con seleccionar_genes_pc()
        L_info_seleccionar_genes = [curr_peak_pc[genes_col], curr_peak_pc[id_genes_col], curr_peak_pc[genes_su_cerca_col], 
                                    curr_peak_pc[id_genes_rnaseq_col], curr_peak_pc[rnaseq_su_cerca_col]]; 
        str_gen_cerca, str_id_genes_rnaseq, str_gen_rnaseq_cerca, str_genes_rnaseq_su_cerca = seleccionar_genes_pc(L_info_seleccionar_genes, sub_sep=sub_sep); 
        # Agrego info de genes a L_out
        L_out.append(str(str_gen_cerca)); 
        L_out.append(str(str_id_genes_rnaseq)); 
        L_out.append(str(str_gen_rnaseq_cerca)); 
        # Solo hago esto la primera vuelta
        if i == 0: 
            L_header = L_header + ['genes_cerca', 'genes_rnaseq', 'genes_rnaseq_cerca']; 
        # Booleano para ver si se encontro nkx25
        nkx25_no_encontrado = True;  
        # Recorro cada uno de los pssm
        for j in range(len(L_scores_pssm)):
            # Defino nombre del TF y limite de score aceptado
            curr_tf = L_nombres_tf[j]; 
            curr_score_limit = L_scores_pssm[j]; 
            # Extraigo la informacion relevante al factor de transcripcion actual
            curr_pos_pssm_tf_col = headers_pc.index('sitios_pssm_' + curr_tf); 
            curr_score_pssm_tf_col = headers_pc.index('scores_pssm_' + curr_tf); 
            curr_pos_pssm_tf = curr_peak_pc[curr_pos_pssm_tf_col]; 
            curr_score_pssm_tf = curr_peak_pc[curr_score_pssm_tf_col]; 
            # Extraigo la informacion relevante a sitios de union con seleccionar_su_pssm()
            str_su_pssm_tf, str_pos_pssm_tf, str_score_pssm_tf = seleccionar_su_pssm_pc(curr_pos_pssm_tf, curr_score_pssm_tf, curr_score_limit, sub_sep=sub_sep); 
            # Reviso si curr_tf es NKX2-5
            if curr_tf.upper() in ['NKX25', 'NKX2-5', 'NKX2.5']:
                # Si lo encuentro, defino las variables relevantes y actualizo nkx25_no_encontrado
                nkx25_no_encontrado = False; 
                sitios_lista = curr_peak_pc[sitios_col]; # 'sitios_lista'
                seq_sitios_lista = curr_peak_pc[seq_sitios_col]; # 'seq_sitios'
                pos_sitios_lista = curr_peak_pc[pos_sitios_col]; # 'pos_sitios'
                genes_su_pssm_nkx25_cerca = curr_peak_pc[id_genes_pssm_nkx25_cerca_col]; # 'id_genes_su_pssm_nkx25_cerca'
                pos_pssm_nkx25_cerca = curr_peak_pc[sitios_pssm_nkx25_cerca_col]; # 'sitios_pssm_gen_cerca'
                score_pssm_nkx25_cerca = curr_peak_pc[scores_pssm_nkx25_cerca_col]; # 'scores_genes_su_pssm_nkx25_cerca'
                # Confirmo que el pico tenga sitio de union lista 
                su_lista = sitios_lista=='SitiosEncontrados'; 
                # Confirmo que el pico tenga sitio de union pssm para NKX2-5
                su_pssm = str_su_pssm_tf=='X'; 
                # Defino picos con sitio confirmado de NKX2-5
                if su_lista or su_pssm:
                    str_sitio_confirmado_nkx25 = 'X'; 
                else:
                    str_sitio_confirmado_nkx25 = ''; 
                # Defino picos con sitio lista de NKX2-5
                if su_lista:
                    str_su_lista = 'X'; 
                else:
                    str_su_lista = ''; 
                # Agrego str_sitio_confirmado_nkx25, str_su_lista, seq_sitios_lista y pos_sitios_lista a L_out
                L_out.append(str(str_sitio_confirmado_nkx25)); 
                L_out.append(str(str_su_lista)); 
                L_out.append(str(seq_sitios_lista)); 
                L_out.append(str(pos_sitios_lista)); 
                # Solo hago esto la primera vuelta
                if i == 0: 
                    L_header = L_header + [curr_tf+'_su_confirmado', 'su_'+curr_tf+'_lista', 'seq_su_lista', 'pos_su_lista']; 
                ### FALTA (para optimizar mas adelante)
                # Confirmar que alguno de los sitios este efectivamente cerca del gen
                ###
                # Defino la clasificacion del pico como "No confirmado", "Regulado directamente" o "Regulado indirectamente"
                if str_gen_rnaseq_cerca=='X':
                    # Genes regulados
                    if su_lista or su_pssm:
                        str_class = '1-Regulado directamente'; 
                    else:
                        str_class = '2-Regulado indirectamente'; 
                else:
                    str_class = '3-No confirmado'; 
                # Agrego str_class en posicion 4 de L_out y L_header
                L_out = L_out[:4] + [str(str_class)] + L_out[4:]; 
                # Solo hago esto la primera vuelta
                if i == 0:  
                    L_header = L_header[:4] + ['peak_class'] + L_header[4:]; 
            # Agrego str_su_pssm_tf despues de la posicion donde va str_class
            if nkx25_no_encontrado:
                num_usado = 4; 
            else:
                num_usado = 5; 
            # Agrego str_su_pssm_tf en la posicion num_usado+j
            L_out = L_out[:num_usado+j] + [str(str_su_pssm_tf)] + L_out[num_usado+j:]; 
            # Solo hago esto la primera vuelta
            if i == 0:  
                L_header = L_header[:num_usado+j] + [curr_tf + '_pssm'] + L_header[num_usado+j:]; 
            # Agrego el resto de la info de sitios de union pssm al final de L_out
            L_out.append(str(str_pos_pssm_tf)); 
            L_out.append(str(str_score_pssm_tf)); 
            # Solo hago esto la primera vuelta
            if i == 0: 
                L_header = L_header + ['pos_pssm_' + curr_tf, 'score_pssm_' + curr_tf]; 
        if nkx25_no_encontrado:
            print('WARNING: No se encontro NKX2-5 en la lista L_nombres_tf.')
            print(L_nombres_tf)
        # Agrego peak id al final
        L_out.append(str(curr_peak_pc[peak_id_col])); 
        # Solo hago esto la primera vuelta
        if i == 0: 
            L_header = L_header + ['peak_id']; 
        # Agrego L_out a M_out
        M_out.append(L_out[:]); 
    return M_out, L_header


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


def pipeline_tabla_final(nom_pc, nom_out, L_scores_pssm, L_nombres_tf, path_pc='', path_out='.\\', verbose=False, display_interval=100):
    # Funcion que corre generar_tabla_final() y guarda el output

    # Uso generar_tabla_final() para crear M_out
    M_out, L_header = generar_tabla_final(nom_pc, L_scores_pssm, L_nombres_tf, path_pc=path_pc, verbose=verbose, display_interval=display_interval); 
    # Guardo M_out en archivo nom_out dentro de path_out
    M_out = guardar_tabla_final(nom_out, M_out, path_out=path_out, L_head=L_header); 
    return M_out


def pipeline_varias_tablas(L_nom_pc, L_nom_out, L_arch_pcm, L_nom_pssm, path_pc='', path_out='.\\', path_pcm='', pseudocounts=0.5, score_fraction=0.9, 
                           verbose=False, display_interval=500):
    # Pipeline para correr varias veces pipeline_tabla_final() 
    # Funciona con fuentes derivadas del mismo .bed

    # Inicializo la matriz que se devuelve
    M_out = []; 
    # Extraigo la lista de scores de cutoff para pssm (compartida entre las corridas)
    L_scores_pssm = buscar_score_cutoff_pssm(L_arch_pcm, path_arch=path_pcm, pseudocounts=pseudocounts, score_fraction=score_fraction); 
    # Recorro cada uno de los elementos en L_nom_pc y L_nom_out
    for i in range(len(L_nom_pc)):
        curr_nom_pc = L_nom_pc[i]; 
        curr_nom_out = L_nom_out[i]; 
        # Corro pipeline_tabla_final()
        M_pipeline = pipeline_tabla_final(curr_nom_pc, curr_nom_out, L_scores_pssm, L_nom_pssm, path_pc=path_pc, path_out=path_out, 
                                          verbose=verbose, display_interval=display_interval); 
        # Agrego M_pipeline a M_out
        M_out.append(M_pipeline[:]); 
    return M_out


def seleccionar_genes_pc(L_info, sub_sep=','):
    # Funcion para seleccionar la informacion relevante a genes cerca de picos de ChIP-seq clasificados
    # Recibe clasificacion por genes cercanos, genes cercanos, genes con sitios (lista) cerca, genes confirmados RNA-seq, genes confirmados RNA-seq con sitios (lista) cerca
    # Devuelve un string ('X'/'F'/'') si hay gen cercano (filtrado o no) o no hay ningun gen cerca
    # Devuelve lista de genes confirmados por RNA-seq
    # Devuelve booleano ('X'/'') si hay gen confirmado por RNA-seq con sitio cercano

    # Extraigo los elementos de L_info
    class_genes = L_info[0]; 
    str_id_genes = L_info[1]; 
    str_id_genes_su_cerca = L_info[2]; 
    str_id_genes_rnaseq = L_info[3]; 
    str_id_genes_rnaseq_su_cerca = L_info[4]; 
    # Reviso class_genes para ver si hay genes cercanos
    if class_genes=='GenesCerca':
        str_gen_cerca = 'X'; 
    elif class_genes=='GenesFiltrados':
        str_gen_cerca = 'F'; 
    elif class_genes=='NoHayGenes':
        str_gen_cerca = ''; 
    else:
        print('ERROR: class_genes tiene un valor inesperado: "' + str(class_genes) + '", se devuelve str_gen_cerca vacio.')
        str_gen_cerca = ''; 
    # Hago una lista con genes de RNA-seq (PARA USAR SI MODIFICO LA FUNCION)
    #L_genes_rnaseq = str_id_genes_rnaseq.rstrip().split(sep=sub_sep); 
    # Reviso str_id_genes_rnaseq para ver si hay genes confirmados por RNA-seq cerca del pico
    if len(str_id_genes_rnaseq):
        str_gen_rnaseq_cerca = 'X'; 
    else:
        str_gen_rnaseq_cerca = ''; 
    # Reviso str_id_genes_rnaseq_su_cerca para ver si hay genes confirmados por RNA-seq cerca de sitios de union confirmados
    if len(str_id_genes_rnaseq_su_cerca):
        str_genes_rnaseq_su_cerca = 'X'; 
    else:
        str_genes_rnaseq_su_cerca = ''; 
    return str_gen_cerca, str_id_genes_rnaseq, str_gen_rnaseq_cerca, str_genes_rnaseq_su_cerca


def seleccionar_su_pssm_pc(str_pos_pssm, str_score_pssm, score_limit, sub_sep=',', max_first=True): 
    # Funcion para seleccionar la informacion relevante a sitios de union por pssm en picos de ChIP-seq clasificados
    # Recibe las posiciones encontradas por pssm, los scores de las posiciones y el limite de score aceptable
    # sub_sep define separador interno de scores y posiciones
    # max_first define si el primer valor es el score maximo
    # Devuelve un booleano ('X'/'') si hay sitio de union que pase score_limit o ninguno lo pasa
    # Devuelve dos strings con los sitios de union que pasan score_limit (sus scores y sus posiciones)

    # Inicializo los strings que se devuelven
    str_su_pssm_out = 'NA'; 
    str_pos_pssm_out = ''; 
    str_score_pssm_out = ''; 
    # Divido str_pos_pssm y str_score_pssm usando sub_sep
    L_pos_pssm_temp = str_pos_pssm.split(sep=sub_sep); 
    L_score_pssm_temp = str_score_pssm.split(sep=sub_sep); 
    # Inicializo las listas para eliminar posiciones repetidas
    L_pos_pssm = []; 
    L_score_pssm = []; 
    # Elimino cualquier aparicion de '' (string vacio) como elemento de L_pos_pssm_temp y L_score_pssm_temp
    if '' in L_pos_pssm_temp:
        ### Display
        if L_pos_pssm_temp != [''] or L_score_pssm_temp != ['']:
            print('WARNING: string vacio aparece en L_pos_pssm_temp, caso no testeado:')
            print(L_pos_pssm_temp)
            print(L_score_pssm_temp)
        ###
        L_pos_pssm_temp.remove(''); 
        L_score_pssm_temp.remove(''); 
        ### Display
        if L_pos_pssm_temp != [] or L_score_pssm_temp != []:
            print('Listas pasadas por remove()')
            print(L_pos_pssm_temp)
            print(L_score_pssm_temp)
        ###
    else:
        L_score_pssm_temp = L_score_pssm_temp[int(max_first):]; 
    # Transformo los valores en float() y elimino repetidos
    for k in range(len(L_pos_pssm_temp)):
        curr_pos_pssm = L_pos_pssm_temp[k]; 
        curr_score_pssm = float(L_score_pssm_temp[k]); 
        # Veo que la posicion no se encuentre en L_pos_pssm
        if not (str(curr_pos_pssm) in L_pos_pssm):
            # Veo que curr_score_pssm sea mayor que score_limit
            if curr_score_pssm > score_limit:
                # Si curr_score_pssm es mayor a score_limit, agrego posicion y score a listas
                L_pos_pssm.append(str(curr_pos_pssm)); 
                L_score_pssm.append(str(curr_score_pssm)); 
                # Tambien las agrego a strings de output y defino str_su_pssm_out como 'X'
                str_pos_pssm_out = str_pos_pssm_out + str(curr_pos_pssm) + sub_sep; 
                str_score_pssm_out = str_score_pssm_out + str(curr_score_pssm) + sub_sep; 
                str_su_pssm_out = 'X'; 
        # Si la posicion se encuentra en L_pos_pssm, veo si es misma combinacion pos/score
        else:
            # Inicializo un booleano para ver si alguna de las combinaciones pos/score se repite
            repetido = False; 
            # Busco la lista de posiciones en L_pos_pssm donde se encuentra str(curr_pos_pssm)
            indices_curr_pos = [n for n, x in enumerate(L_pos_pssm) if x == str(curr_pos_pssm)]; 
            # Veo si curr_score_pssm para alguno de los n en indices_curr_pos da mismo score en L_score_pssm
            for index_curr_pos in indices_curr_pos:
                # Defino el score correspondiente a curr_pos_pssm para este index
                corresponding_score = L_score_pssm[index_curr_pos]; 
                # Reviso si corresponding_score es igual a curr_score_pssm
                if corresponding_score == curr_score_pssm: 
                    # Si es asi, defino repetido como True
                    repetido = True; 
            # Si despues de recorrer todos los indices_curr_pos repetido sigue siendo False, reviso si el score es mayor a score_limit
            if not repetido:
                # Veo que curr_score_pssm sea mayor que score_limit
                if curr_score_pssm > score_limit:
                    # Agrego los valores a listas y strings
                    L_pos_pssm.append(str(curr_pos_pssm)); 
                    L_score_pssm.append(str(curr_score_pssm)); 
                    str_pos_pssm_out = str_pos_pssm_out + str(curr_pos_pssm) + sub_sep; 
                    str_score_pssm_out = str_score_pssm_out + str(curr_score_pssm) + sub_sep; 
                    str_su_pssm_out = 'X'; 
    # Si despues de recorrer la lista de sitios, no se encuentra ninguno con score mayor a score_limit, cambio str_su_pssm_out a ''
    if str_su_pssm_out == 'NA':
        str_su_pssm_out = ''; 
    # Tambien elimino las comas al final de los otros dos strings
    str_pos_pssm_out = str_pos_pssm_out.rstrip(sub_sep); 
    str_score_pssm_out = str_score_pssm_out.rstrip(sub_sep); 
    return str_su_pssm_out, str_pos_pssm_out, str_score_pssm_out


def _main():
    # Funcion para probar funciones en ejecucion del archivo

    # 
    M_out = []; 
    # Corro pipeline_varias_tablas() para raton y humano
    M_raton = pipeline_varias_tablas(L_pc_dupays, L_out_dupays, L_arch_pcm_raton, L_nombres_pssm_raton, path_pc=path_pc_main, path_out=path_out_main, 
                                     path_pcm=path_pwm_mouse, pseudocounts=pseudocounts_main, score_fraction=score_fraction_main, verbose=verbose_main, 
                                     display_interval=display_interval_main); 
    M_humano = pipeline_varias_tablas(L_pc_anderson, L_out_anderson, L_arch_pcm_humano, L_nombres_pssm_humano, path_pc=path_pc_main, path_out=path_out_main, 
                                      path_pcm=path_pwm_human, pseudocounts=pseudocounts_main, score_fraction=score_fraction_main, verbose=verbose_main, 
                                      display_interval=display_interval_main); 
    M_out.append(M_raton); 
    M_out.append(M_humano); 
    return M_out


#################################### RESULTADOS ###################################

output_dump = []; 

if __name__=='__main__':
    output_dump.append(_main()); 

