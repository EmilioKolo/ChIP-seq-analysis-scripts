
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
Funcion para seleccionar picos de ChIP-seq y hacer archivos .fasta para mandar a MEME-ChIP

### FALTA:
# Agarrar clasificacion de picos (ver si hay tabla en archivos del congreso)
# Agarrar secuencias con seq_data
# Generar archivos .fasta correspondientes
    # Picos sin sitios de union NKX2-5
    # Picos con sitios de union NKX2-5
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
nom_in_dupays = '13-Dupays100k'; 
nom_in_anderson = '13-Anderson100k'; 


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
path_in_main = path_dropbox_main + 'Doctorado\\3-Genes transactivados rio abajo\\1-Interpretacion de datos ChIP-seq\\'; 
path_out_main = path_output_dump_main + 'FastaMEME\\PostCorrientes\\'; 
# Path que dependen de la especie
path_pwm_human = path_dropbox_main + 'Doctorado\\3-Genes transactivados rio abajo\\1-Interpretacion de datos ChIP-seq\\PWM_human\\'; 
path_pwm_mouse = path_dropbox_main + 'Doctorado\\3-Genes transactivados rio abajo\\1-Interpretacion de datos ChIP-seq\\PWM_mouse\\'; 

# Variables para generar .fasta
largo_main = 1500; 


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


def guardar_fasta(M_fasta, nom_out, path_out=''):
    # Guarda los elementos de M_fasta en archivo nom_out .fasta ubicado en path_out
    # Pensado para funcionar con M_fasta sacado de secuencias_L_rangos

    # Defino dir_arch
    if path_out == '':
        dir_arch = nom_out + '.fasta'; 
    else:
        dir_arch = os.path.join(path_out, nom_out + '.fasta'); 
    # Abro el archivo para crearlo si no existe o sobreescribirlo si ya existe
    with open(dir_arch, 'w') as F_out:
        print('Archivo ' + nom_out + '.fasta creado.')
    # Vuelvo a abrir el archivo en modo append
    with open(dir_arch, 'a') as F_out:
        # Recorro M_fasta
        for i in M_fasta:
            # Guardo cada fila en F_out
            F_out.write(str(i)); 
    return M_fasta


def secuencias_L_rangos(M_peaks, genome_name, L_col_rangos=[0,1,2], path_fasta='', largo_max=0, largo_min=10, recortar_rango=False, largo_target=500, 
                        verbose=False, rango_print=25):
    # Obtiene las secuencias de una lista de rangos
    # Recibe una matriz con la info filtrada por seleccionar_rango() y agarra las posiciones L_col_rangos para armar el rango

    ### Display
    cont_print = 0; 
    len_print = len(M_peaks); 
    if verbose:
        print('Inicializando recopilacion de ' + str(len_print) + ' secuencias.')
    ###
    # Inicializo la matriz que se devuelve
    M_fasta = []; 
    # Inicializo un contador de id para el principio del fasta
    cont_fasta = 0; 
    # Recorro M_peaks
    for i in range(len(M_peaks)):
        L_peak = M_peaks[i]; 
        # Defino el rango en base a L_col_rangos
        curr_rango = []; 
        for j in L_col_rangos:
            curr_rango.append(L_peak[j]); 
        # Si L_col es lista vacia, uso [0,1,2]
        if len(L_col_rangos) == 0: 
            curr_rango = [L_peak[0], int(L_peak[1]), int(L_peak[2])]; 
        # Si L_col_rangos es de largo menor a 3, tiro warning y actuo como si la lista estuviera vacia
        elif len(L_col_rangos) < 3: 
            print('ERROR: L_col_rangos tienen menos de 3 elementos, se usa [0,1,2].')
            curr_rango = [L_peak[0], int(L_peak[1]), int(L_peak[2])]; 
        # Si L_col_rangos fue usado correctamente, transformo posiciones 1 y 2 en int() para pos_ini y pos_end
        else:
            curr_rango[1] = int(curr_rango[1]); 
            curr_rango[2] = int(curr_rango[2]); 
        # Recorto el rango si recortar_rango es True
        if recortar_rango:
            new_ini, new_end = _rango_medio(curr_rango[1], curr_rango[2], largo_target); 
            curr_rango[1] = new_ini; 
            curr_rango[2] = new_end; 
        # Solo reviso largo_max si es un valor mayor a 0
        elif largo_max > 0:
            # Recorto los rangos mayores a largo_max
            if (max(curr_rango[1],curr_rango[2])-min(curr_rango[1],curr_rango[2])) > largo_max:
                ### Display
                if verbose:
                    print('Rango ' + str(curr_rango) + ' de largo mayor a ' + str(largo_max))
                ###
                new_ini, new_end = _rango_medio(curr_rango[1], curr_rango[2], largo_max); 
                curr_rango[1] = new_ini; 
                curr_rango[2] = new_end; 
        # Alargo los rangos menores a largo_min
        if (max(curr_rango[1],curr_rango[2])-min(curr_rango[1],curr_rango[2])) < largo_min:
            ### Display
            if verbose: 
                print('Rango ' + str(curr_rango) + ' de largo menor a ' + str(largo_min))
            ###
            new_ini, new_end = _rango_medio(curr_rango[1], curr_rango[2], largo_min); 
            curr_rango[1] = new_ini; 
            curr_rango[2] = new_end; 
        # Creo elemento seq_data para buscar secuencias
        seq_element = seq_data(genome_name, path_fasta=path_fasta); 
        # Consigo la secuencia del rango
        curr_seq = seq_element._consulta_secuencia_fasta(curr_rango[0], curr_rango[1], curr_rango[2]); 
        # Agrego 1 al contador fasta
        cont_fasta += 1; 
        # Defino un string formato .fasta
        str_fasta = '>ID' + str(cont_fasta) + ', ' + str(curr_rango[0]) + ', ' + str(curr_rango[1]) + ', ' + str(curr_rango[2]) + '\n' + str(curr_seq) + '\n'; 
        # Agrego str_fasta a M_fasta
        M_fasta.append(str(str_fasta)); 
        ### Display
        cont_print += 1; 
        if verbose and (((cont_print+1)%rango_print)==0):
            print('Progreso: ' + str(cont_print+1) + ' / ' + str(len_print))
        ###
    return M_fasta


def seleccionar_peaks(M_in):
    # Funcion para seleccionar los grupos de peaks para hacer los archivos fasta

    # Inicializo las matrices que se devuelven
    M_out_con_su, M_out_sin_su = [], []; 

    # Recorro la matriz de input
    for i in range(len(M_in)):
        curr_peak = M_in[i]; 
        # Inicializo la lista de datos relevantes
        L_out = []; 
        # Agrego los datos relevantes
        L_out.append(curr_peak[0]); # chr_n
        #L_out.append(curr_peak[1]); # n (sin 'chr')
        L_out.append(curr_peak[2]); # pos_ini
        L_out.append(curr_peak[3]); # pos_end
        # Veo si el peak tiene sitio de union de NKX2-5
        peak_con_su = curr_peak[5]=='X'; 
        if peak_con_su:
            # Si tiene sitio de union, agrego L_out a M_out_con_su
            M_out_con_su.append(L_out[:]); 
        else:
            # Si no tiene sitio de union, agrego L_out a M_out_sin_su
            M_out_sin_su.append(L_out[:]); 
    return M_out_con_su, M_out_sin_su


def seleccionar_peaks_rnaseq(M_in):
    # Funcion para seleccionar los grupos de peaks para hacer los archivos fasta

    # Inicializo las matrices que se devuelven
    M_out_con_su, M_out_sin_su = [], []; 

    # Recorro la matriz de input
    for i in range(len(M_in)):
        curr_peak = M_in[i]; 
        # Inicializo la lista de datos relevantes
        L_out = []; 
        # Agrego los datos relevantes
        L_out.append(curr_peak[0]); # chr_n
        #L_out.append(curr_peak[1]); # n (sin 'chr')
        L_out.append(curr_peak[2]); # pos_ini
        L_out.append(curr_peak[3]); # pos_end
        # Veo si el peak esta confirmado por RNAseq
        peak_confirmado_rnaseq = curr_peak[18]=='X'; 
        if peak_confirmado_rnaseq:
            # Veo si el peak tiene sitio de union de NKX2-5
            peak_con_su = curr_peak[5]=='X'; 
            if peak_con_su:
                # Si tiene sitio de union, agrego L_out a M_out_con_su
                M_out_con_su.append(L_out[:]); 
            else:
                # Si no tiene sitio de union, agrego L_out a M_out_sin_su
                M_out_sin_su.append(L_out[:]); 
    return M_out_con_su, M_out_sin_su


def _main():
    # Funcion para probar funciones en ejecucion del archivo

    # Abro los archivos csv con las clasificaciones de los peaks
    M_in_mouse = abrir_csv(nom_in_dupays, dir_arch=path_in_main, ignore_headers=True); 
    M_in_human = abrir_csv(nom_in_anderson, dir_arch=path_in_main, ignore_headers=True); 

    # Selecciono los grupos de peaks con y sin sitios de union de NKX2-5
    #M_con_su_mouse, M_sin_su_mouse = seleccionar_peaks(M_in_mouse); 
    #M_con_su_human, M_sin_su_human = seleccionar_peaks(M_in_human); 
    M_con_su_mouse, M_sin_su_mouse = seleccionar_peaks_rnaseq(M_in_mouse); 
    M_con_su_human, M_sin_su_human = seleccionar_peaks_rnaseq(M_in_human); 

    # Paso las 4 matrices por secuencias_L_rangos() y guardar_fasta()
    fasta_con_su_mouse = secuencias_L_rangos(M_con_su_mouse, 'mm9', path_fasta=path_fasta_main, recortar_rango=True, largo_target=largo_main, verbose=True); 
    fasta_con_su_mouse = guardar_fasta(fasta_con_su_mouse, 'ConSUDupaysLargo1500_rnaseq', path_out=path_out_main); 
    fasta_sin_su_mouse = secuencias_L_rangos(M_sin_su_mouse, 'mm9', path_fasta=path_fasta_main, recortar_rango=True, largo_target=largo_main, verbose=True); 
    fasta_sin_su_mouse = guardar_fasta(fasta_sin_su_mouse, 'SinSUDupaysLargo1500_rnaseq', path_out=path_out_main); 
    fasta_con_su_human = secuencias_L_rangos(M_con_su_human, 'hg19', path_fasta=path_fasta_main, recortar_rango=True, largo_target=largo_main, verbose=True); 
    fasta_con_su_human = guardar_fasta(fasta_con_su_human, 'ConSUAndersonLargo1500_rnaseq', path_out=path_out_main); 
    fasta_sin_su_human = secuencias_L_rangos(M_sin_su_human, 'hg19', path_fasta=path_fasta_main, recortar_rango=True, largo_target=largo_main, verbose=True); 
    fasta_sin_su_human = guardar_fasta(fasta_sin_su_human, 'SinSUAndersonLargo1500_rnaseq', path_out=path_out_main); 

    return ''


def _rango_medio(pos_ini, pos_end, largo_target, extend_end=True, verbose=False):
    # Devuelve un rango de largo largo_target alrededor del centro de pos_ini y pos_end
    
    # Defino paridad del rango dado y de largo_target (1 impar, 0 par)
    len_range = (max(pos_end,pos_ini) - min(pos_end,pos_ini))+1; 
    parity_range = len_range%2; 
    parity_target = largo_target%2; 
    ### Display
    if verbose:
        print('parity_range: ' + str(parity_range))
        print('parity_target: ' + str(parity_target))
    ###
    # Si el rango es par
    if not parity_range:
        # Defino ambos lados del punto medio
        center_ini = (pos_ini+pos_end-1)/2; 
        center_end = (pos_ini+pos_end+1)/2; 
        # Defino la mitad de largo_target redondeado
        half_target_round = ((largo_target-parity_target)/2)-1; 
        # Sumo/Resto la mitad de largo_target redondeado a ambos lados del punto medio
        if extend_end:
            pos_ini_out = center_ini - half_target_round; 
            pos_end_out = center_end + half_target_round + parity_target; 
        else:
            pos_ini_out = center_ini - half_target_round - parity_target; 
            pos_end_out = center_end + half_target_round; 
        ### Display
        if verbose:
            print('center_ini: ' + str(center_ini))
            print('center_end: ' + str(center_end))
            print('half_target_round: ' + str(half_target_round))
            print('pos_ini_out: ' + str(pos_ini_out))
            print('pos_end_out: ' + str(pos_end_out))
        ###
    # Si el rango es impar
    else:
        # Defino el punto medio
        center = (pos_ini+pos_end)/2; 
        # Defino cuanto crece el rango a cada lado del centro
        half_target_center = ((largo_target+parity_target)/2)-1; 
        # Sumo/Resto el valor definido a cada lado del centro
        if extend_end:
            pos_ini_out = center - half_target_center; 
            pos_end_out = center + half_target_center + (parity_target+1)%2; 
        else:
            pos_ini_out = center - half_target_center - (parity_target+1)%2; 
            pos_end_out = center + half_target_center; 
        ### Display
        if verbose:
            print('center: ' + str(center))
            print('half_target_center: ' + str(half_target_center))
            print('pos_ini_out: ' + str(pos_ini_out))
            print('pos_end_out: ' + str(pos_end_out))
        ###
    return int(pos_ini_out), int(pos_end_out)


#################################### RESULTADOS ###################################

output_dump = []; 

if __name__=='__main__':
    output_dump.append(_main()); 

