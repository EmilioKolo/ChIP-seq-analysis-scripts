
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
path_analisischip_ib3 = ''; 
sys.path.append(path_analisischip_casa); 
from analisischip.seq import seq_data

# Genomas de referencia
hg19 = EnsemblRelease(75, species='human'); 
hg38 = EnsemblRelease(102, species='human'); 
mm9 = EnsemblRelease(54, species='mouse'); 
mm10 = EnsemblRelease(102, species='mouse'); 


'''
Funciones para generar archivos .fasta que puedan ser usados para MEME y MEME-ChIP

Pipeline empieza con archivos generados por 3-ClasificacionPeaks.py

Funcion pipeline que corre todo en orden: pipeline_fasta()

Funcion para filtrar peaks segun si tienen sitio de union y/o genes cerca: filtrar_L_rangos()

Funcion para transformar rangos en secuencias: secuencias_L_rangos()
    Incluye opcion de hacer rangos del mismo largo: rango_medio()

Funcion para guardar salida de secuencias_L_rangos() como archivo .fasta: guardar_fasta()

### FALTA:
~ Probar funcion secuencias_L_rangos()
# Hacer funcion de testing para medir rangos con media de largo, varianza, histograma, etc.
X Probar pipeline_fasta()
# Conseguir secuencias correspondientes a rangos
    # Secuencias completas de peaks para MEME
        # Con genes cerca con sitio de union NKX2-5
        # Con genes cerca sin sitio de union NKX2-5
        # Sin genes cerca con sitio de union NKX2-5
        # Sin genes cerca sin sitio de union NKX2-5
        # Con genes cerca
        # Sin genes cerca
        # Sin sitio de union NKX2-5
    # Secuencias de mismo largo (idealmente 500, maximo 1000~1500)
        # Alrededor del medio para sin SU
        # Usar posiciones de SU para con SU
X Definir listas de peaks para correr en MEME-ChIP
X Usar abrir_csv() para agarrar archivos con peaks clasificados
X Seleccionar peaks (chr_n/contig, pos_ini, pos_end)
    X pos_ini/pos_end pueden ser modificados al agarrar
X Crear funcion pipeline_fasta() que agrupe abrir_csv(), filtrar_L_rangos(), secuencias_L_rangos() y guardar_fasta()
###
'''

#################################### VARIABLES ####################################

# Direcciones de archivos y outputs
path_peaks_casa = 'D:\\Archivos doctorado\\Output_dump\\PeaksClasificados\\'; 
path_peaks_ib3 = 'X:\\Output_dump\\PeaksClasificados\\'; 
path_out_casa = 'D:\\Archivos doctorado\\Output_dump\\FastaMEME\\'; 
path_out_ib3 = 'X:\\Output_dump\\FastaMEME\\'; 
path_fasta_casa = 'D:\\Archivos doctorado\\Genomas\\'; 
path_fasta_ib3 = 'X:\\Genomas\\'; 

# Nombres de archivos
L_arch_anderson = ['AndersonClassPeaks_SitiosConf_dist1500', 'AndersonClassPeaks_SitiosConf_dist50k', 'AndersonClassPeaks_SitiosConf_dist1M']; 
L_arch_dupays = ['DupaysClassPeaks_SitiosConf_dist1500', 'DupaysClassPeaks_SitiosConf_dist50k', 'DupaysClassPeaks_SitiosConf_dist1M']; 

## Datos de input
# Con genes cerca con sitio de union NKX2-5
dict_filtros_con_genes_con_su = {'filtrar_su':True, 'con_su':True, 'filtrar_genes':True, 'con_genes':True, 'genes_filtrados':False, 'sitios_cerca_genes':True}; # 'ConGenesConSU'
# Con genes cerca sin sitio de union NKX2-5
dict_filtros_con_genes_sin_su = {'filtrar_su':True, 'con_su':False, 'filtrar_genes':True, 'con_genes':True, 'genes_filtrados':False, 'sitios_cerca_genes':True}; # 'ConGenesSinSU'
# Sin sitio de union NKX2-5
dict_filtros_solo_sin_su = {'filtrar_su':True, 'con_su':False, 'filtrar_genes':False, 'con_genes':True, 'genes_filtrados':False, 'sitios_cerca_genes':True}; # 'SoloSinSU'
# Sin genes cerca con sitio de union NKX2-5
dict_filtros_sin_genes_con_su = {'filtrar_su':True, 'con_su':True, 'filtrar_genes':True, 'con_genes':False, 'genes_filtrados':False, 'sitios_cerca_genes':True}; # 'SinGenesConSU'
# Sin genes cerca sin sitio de union NKX2-5
dict_filtros_sin_genes_sin_su = {'filtrar_su':True, 'con_su':False, 'filtrar_genes':True, 'con_genes':False, 'genes_filtrados':False, 'sitios_cerca_genes':True}; # 'SinGenesSinSU'
# Con genes cerca
dict_filtros_solo_con_genes = {'filtrar_su':False, 'con_su':True, 'filtrar_genes':True, 'con_genes':True, 'genes_filtrados':False, 'sitios_cerca_genes':True}; # 'SoloConGenes'
# Sin genes cerca
dict_filtros_solo_sin_genes = {'filtrar_su':False, 'con_su':True, 'filtrar_genes':True, 'con_genes':False, 'genes_filtrados':False, 'sitios_cerca_genes':True}; # 'SoloSinGenes'
    
L_dict_filtros_total = [dict_filtros_con_genes_con_su, dict_filtros_con_genes_sin_su, dict_filtros_solo_sin_su, dict_filtros_sin_genes_con_su, dict_filtros_sin_genes_sin_su, 
                        dict_filtros_solo_con_genes, dict_filtros_solo_sin_genes]; 
L_dict_nom_total = ['ConGenesConSU', 'ConGenesSinSU', 'SoloSinSU', 'SinGenesConSU', 'SinGenesSinSU', 'SoloConGenes', 'SoloSinGenes']; 

L_dict_filtros = [dict_filtros_con_genes_sin_su, dict_filtros_solo_sin_su, dict_filtros_sin_genes_con_su, dict_filtros_sin_genes_sin_su, 
                  dict_filtros_solo_con_genes, dict_filtros_solo_sin_genes]; 
L_dict_nom = ['ConGenesSinSU', 'SoloSinSU', 'SinGenesConSU', 'SinGenesSinSU', 'SoloConGenes', 'SoloSinGenes']; 

# Variables main()

path_fasta_main = path_fasta_casa; 
path_peaks_main = path_peaks_casa; 
path_out_main = path_out_casa; 

rango_print_main = 100; 
L_dict_filtros_main = L_dict_filtros_total; 
L_dict_nom_main = L_dict_nom_total; 



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


def filtrar_L_rangos(M_csv, filtrar_su=True, con_su=True, filtrar_genes=True, con_genes=True, genes_filtrados=False, sitios_cerca_genes=True, L_col=[0,2,3], class_col=9):
    # Recibe matriz de archivo generado en 3-ClasificacionPeaks y selecciona de acuerdo a si tiene genes y/o sitios de union cerca
    # Devuelve lista de rangos con parte de la informacion (definido por L_col)

    # Inicializo la matriz que se devuelve
    M_out = []; 
    # Recorro M_csv
    for curr_rango in M_csv:
        # Hago el filtro por rango en filtrar_rango()
        if filtrar_rango(curr_rango, filtrar_su, con_su, filtrar_genes, con_genes, genes_filtrados, sitios_cerca_genes, class_col):
            # Defino el formato de rango a guardar
            rango_guardado = seleccionar_rango(curr_rango, L_col=L_col); 
            # Agrego rango_guardado a M_out
            M_out.append(rango_guardado[:]); 
    return M_out


def filtrar_rango(curr_rango, filtrar_su, con_su, filtrar_genes, con_genes, genes_filtrados, sitios_cerca_genes, class_col):
    # Define si curr_rango se tiene que seleccionar o filtrar
    ''' Posibles clasificaciones
GenConSitioCercano
GenConSitioLejos
GenSinSitio
FiltradoConSitioCercano
FiltradoConSitioLejos
FiltradoSinSitio
SitiosSinGen
SinSitiosNiGen
    '''
    # Inicializo el booleano que se devuelve
    seleccionado = True; 
    # Defino el valor de class para curr_rango
    class_rango = curr_rango[class_col]; 
    # Si se filtra por sitio de union y genes
    if filtrar_su and filtrar_genes:
        # Si se seleccionan rangos con sitios de union
        if con_su:
            # Si se seleccionan rangos con genes cerca y con sitios de union
            if con_genes:
                # Si se tienen en cuenta los genes filtrados
                if genes_filtrados:
                    if sitios_cerca_genes:
                        seleccionado = class_rango in ['GenConSitioCercano', 'FiltradoConSitioCercano']; 
                    else:
                        seleccionado = class_rango in ['GenConSitioCercano', 'GenConSitioLejos', 'FiltradoConSitioCercano', 'FiltradoConSitioLejos']; 
                # Si no se tienen en cuenta los genes filtrados
                else:
                    if sitios_cerca_genes:
                        seleccionado = class_rango in ['GenConSitioCercano']; 
                    else:
                        seleccionado = class_rango in ['GenConSitioCercano', 'GenConSitioLejos']; 
            # Si se seleccionan rangos sin genes cerca y con sitios de union
            else:
                seleccionado = class_rango in ['SitiosSinGen']; 
        # Si se seleccionan rangos sin sitios de union
        else:
            # Si se seleccionan rangos con genes cerca y sin sitios de union
            if con_genes:
                # Si se tienen en cuenta los genes filtrados
                if genes_filtrados:
                    seleccionado = class_rango in ['GenSinSitio', 'FiltradoSinSitio']; 
                # Si no se tienen en cuenta los genes filtrados
                else:
                    seleccionado = class_rango in ['GenSinSitio']; 
            # Si se seleccionan rangos sin genes cerca y sin sitios de union
            else:
                seleccionado = class_rango in ['SinSitiosNiGen']; 
    # Si se filtra solo por sitio de union
    elif filtrar_su:
        # Si se seleccionan rangos con sitio de union
        if con_su:
            seleccionado = class_rango in ['GenConSitioCercano', 'GenConSitioLejos', 'FiltradoConSitioCercano', 'FiltradoConSitioLejos', 'SitiosSinGen']; 
        # Si se seleccionan rangos sin sitio de union
        else:
            seleccionado = class_rango in ['GenSinSitio', 'FiltradoSinSitio', 'SinSitiosNiGen']; 
    # Si se filtra solo por genes cercanos
    elif filtrar_genes:
        # Si se seleccionan rangos con genes cerca
        if con_genes:
            # Si se tienen en cuenta los genes filtrados
            if genes_filtrados:
                seleccionado = class_rango in ['FiltradoSinSitio', 'GenSinSitio', 'GenConSitioCercano', 'GenConSitioLejos', 'FiltradoConSitioCercano', 'FiltradoConSitioLejos']; 
            # Si no se tienen en cuenta los genes filtrados
            else:
                seleccionado = class_rango in ['GenSinSitio', 'GenConSitioCercano', 'GenConSitioLejos']; 
        # Si se seleccionan rangos sin genes cerca
        else:
            seleccionado = class_rango in ['SitiosSinGen', 'SinSitiosNiGen']; 
    else:
        print('ERROR: filtrar_genes y filtrar_su son False')
        seleccionado = True; 
    return seleccionado


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


def pipeline_fasta(arch_peaks, dict_filtros, genome_name, nom_out, path_fasta='', path_peaks='', path_out='', recortar_rango=False, 
                   largo_target=500, largo_max=0, verbose=True, rango_print=25):
    # Pipeline que corre abrir_csv(), filtrar_L_rangos(), secuencias_L_rangos() y guardar_fasta() en orden

    # Abro archivos de entrada con todos los peaks
    M_bed = abrir_csv(arch_peaks, path_peaks); 
    ### Display
    if verbose:
        print('Empezando filtrar_L_rangos()')
    ###
    # Cargo las variables de dict_filtros
    filtrar_su = dict_filtros['filtrar_su']; 
    con_su = dict_filtros['con_su']; 
    filtrar_genes = dict_filtros['filtrar_genes']; 
    con_genes = dict_filtros['con_genes']; 
    genes_filtrados = dict_filtros['genes_filtrados']; 
    sitios_cerca_genes = dict_filtros['sitios_cerca_genes']; 
    # Variables predefinidas
    L_col=[0,2,3]; 
    class_col = 9; 
    # Uso filtrar_L_rangos() para seleccionar peaks de acuerdo a dict_filtros
    M_filtrado = filtrar_L_rangos(M_bed, filtrar_su=filtrar_su, con_su=con_su, filtrar_genes=filtrar_genes, con_genes=con_genes, genes_filtrados=genes_filtrados, 
                                  sitios_cerca_genes=sitios_cerca_genes, L_col=L_col, class_col=class_col); 
    ### Display
    if verbose:
        print('Empezando secuencias_L_rangos()')
    ###
    # Uso secuencias_L_rango() para conseguir las secuencias de los rangos en formato .fasta
    M_fasta = secuencias_L_rangos(M_filtrado, genome_name=genome_name, path_fasta=path_fasta, verbose=verbose, recortar_rango=recortar_rango, 
                                  largo_target=largo_target, largo_max=largo_max, rango_print=rango_print); 
    ### Display
    if verbose:
        print('Empezando guardar_fasta()')
    # Guardo M_fasta
    M_fasta = guardar_fasta(M_fasta, nom_out=nom_out, path_out=path_out); 
    return M_fasta


def fasta_total(arch_peaks, genome_name, nom_out, L_largos_target, largo_max=0, path_fasta='', path_peaks='', path_out='', verbose=True, rango_print=25):
    # Pipeline que genera archivos .fasta para todos los peaks de un archivo .bed
    # Genera uno para los largos completos y uno mas por numero en L_largos_target con largos definidos

    # Abro el archivo de entrada con todos los peaks
    M_bed = abrir_csv(arch_peaks, path_peaks); 
    # Defino L_col_rangos
    L_col_rangos = [0,2,3]; 
    ### Display
    if verbose:
        print('Empezando secuencias_L_rangos() para largo total.')
    ###
    # Uso secuencias_L_rango() para conseguir las secuencias de los rangos en formato .fasta
    M_out = secuencias_L_rangos(M_bed, genome_name=genome_name, path_fasta=path_fasta, verbose=verbose, recortar_rango=False, 
                                largo_max=largo_max, rango_print=rango_print, L_col_rangos=L_col_rangos); 
    ### Display
    if verbose:
        print('Empezando guardar_fasta() para largo total.')
    ###
    # Guardo M_out
    M_out = guardar_fasta(M_out, nom_out=nom_out + 'Total', path_out=path_out); 
    print()
    # Repito por elemento en L_largos_target
    for largo_target in L_largos_target:
        ### Display
        if verbose:
            print('Empezando secuencias_L_rangos() para largo ' + str(largo_target))
        ###
        # Uso secuencias_L_rango() para conseguir las secuencias de los rangos en formato .fasta
        M_fasta = secuencias_L_rangos(M_bed, genome_name=genome_name, path_fasta=path_fasta, verbose=verbose, recortar_rango=True, 
                                      largo_target=largo_target, rango_print=rango_print, L_col_rangos=L_col_rangos); 
        ### Display
        if verbose:
            print('Empezando guardar_fasta() para largo ' + str(largo_target))
        ###
        # Guardo M_fasta
        M_fasta = guardar_fasta(M_fasta, nom_out=nom_out + 'Largo' + str(largo_target), path_out=path_out); 
        print()
    return M_out


def rango_medio(pos_ini, pos_end, largo_target, extend_end=True, verbose=False):
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


def secuencias_L_rangos(M_bed, genome_name, L_col_rangos=[0,1,2], path_fasta='', largo_max=0, largo_min=10, recortar_rango=False, largo_target=500, 
                        verbose=False, rango_print=25):
    # Obtiene las secuencias de una lista de rangos
    # Recibe una matriz con la info filtrada por seleccionar_rango() y agarra las posiciones L_col_rangos para armar el rango

    ### Display
    cont_print = 0; 
    len_print = len(M_bed); 
    if verbose:
        print('Inicializando recopilacion de ' + str(len_print) + ' secuencias.')
    ###
    # Inicializo la matriz que se devuelve
    M_fasta = []; 
    # Inicializo un contador de id para el principio del fasta
    cont_fasta = 0; 
    # Recorro M_bed
    for i in range(len(M_bed)):
        L_bed = M_bed[i]; 
        # Defino el rango en base a L_col_rangos
        curr_rango = []; 
        for j in L_col_rangos:
            curr_rango.append(L_bed[j]); 
        # Si L_col es lista vacia, uso [0,1,2]
        if len(L_col_rangos) == 0: 
            curr_rango = [L_bed[0], int(L_bed[1]), int(L_bed[2])]; 
        # Si L_col_rangos es de largo menor a 3, tiro warning y actuo como si la lista estuviera vacia
        elif len(L_col_rangos) < 3: 
            print('ERROR: L_col_rangos tienen menos de 3 elementos, se usa [0,1,2].')
            curr_rango = [L_bed[0], int(L_bed[1]), int(L_bed[2])]; 
        # Si L_col_rangos fue usado correctamente, transformo posiciones 1 y 2 en int() para pos_ini y pos_end
        else:
            curr_rango[1] = int(curr_rango[1]); 
            curr_rango[2] = int(curr_rango[2]); 
        # Recorto el rango si recortar_rango es True
        if recortar_rango:
            new_ini, new_end = rango_medio(curr_rango[1], curr_rango[2], largo_target); 
            curr_rango[1] = new_ini; 
            curr_rango[2] = new_end; 
        # Solo reviso largo max si es un valor mayor a 0
        elif largo_max > 0:
            # Recorto los rangos mayores a largo_max
            if (max(curr_rango[1],curr_rango[2])-min(curr_rango[1],curr_rango[2])) > largo_max:
                ### Display
                if verbose:
                    print('Rango ' + str(curr_rango) + ' de largo mayor a ' + str(largo_max))
                ###
                new_ini, new_end = rango_medio(curr_rango[1], curr_rango[2], largo_max); 
                curr_rango[1] = new_ini; 
                curr_rango[2] = new_end; 
        # Alargo los rangos menores a largo_min
        if (max(curr_rango[1],curr_rango[2])-min(curr_rango[1],curr_rango[2])) < largo_min:
            ### Display
            if verbose: 
                print('Rango ' + str(curr_rango) + ' de largo menor a ' + str(largo_min))
            ###
            new_ini, new_end = rango_medio(curr_rango[1], curr_rango[2], largo_min); 
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


def seleccionar_rango(curr_rango, L_col=[0,2,3]):
    # Selecciona columnas de curr_rango y las devuelve como una lista

    # Inicializo la lista que se devuelve
    L_out = []; 
    # Recorro L_col
    for i in L_col:
        L_out.append(curr_rango[i]); 
    # Si L_col es lista vacia, agrego todo el rango
    if L_col == []:
        L_out = curr_rango[:]; 
    return L_out


def str_de_L(i, L):
    return str(L[i])


def _main():
    # Funcion para probar funciones en ejecucion del archivo

    # Inicializo M_out
    M_out = []; 

    # Corro fasta_total() para generar los archivos de todos los peaks para Dupays y Anderson
    M_fasta = fasta_total('DupaysClassPeaks_SitiosConf_dist1M','mm9','Dupays', [500,1500], path_fasta=path_fasta_main, path_peaks=path_peaks_main, 
                          path_out=path_out_main, rango_print=rango_print_main); 
    M_fasta = fasta_total('AndersonClassPeaks_SitiosConf_dist1M','hg19','Anderson', [500,1500], path_fasta=path_fasta_main, path_peaks=path_peaks_main, 
                          path_out=path_out_main, rango_print=rango_print_main); 

    '''# Recorro los 6 posibles nombres de archivo
    for i in range(6):
        # Inicializo nom_out
        nom_out = ''; 
        # Preparo lista de nombres para Dupays y Anderson
        L = ['1500', '50k', '1M']; 
        # Si es menor a 3 es Dupays
        if i < 3:
            nom_arch = L_arch_dupays[i]; 
            genome_name = 'mm9'; 
            # Agrego Dupays a nom_out
            nom_out = nom_out + 'Dupays' + str_de_L(i, L); 
        # Si es mayor o igual a 3 es Anderson
        else:
            nom_arch = L_arch_anderson[i-3]; 
            genome_name = 'hg19'; 
            # Agrego Anderson a nom_out
            nom_out = nom_out + 'Anderson' + str_de_L(i-3, L); 
        # Recorro L_dict_filtros
        for j in range(len(L_dict_filtros_main)):
            # Re-inicializo nom_out
            full_nom_out = ''; 
            # Defino dict_filtros y dict_nom
            dict_filtros = L_dict_filtros_main[j]; 
            dict_nom = L_dict_nom_main[j]; 
            # Agrego dict_nom a nom_out
            full_nom_out = dict_nom + nom_out; 
            M_fasta = pipeline_fasta(nom_arch, dict_filtros, genome_name, full_nom_out, path_fasta=path_fasta_main, path_peaks=path_peaks_main, 
                                     path_out=path_out_main, largo_max=2000, rango_print=rango_print_main, largo_target=500, recortar_rango=True); 
            print()'''

    '''# Prueba pipeline a mano
    # Abro archivos 
    M_bed = abrir_csv(L_arch_anderson[0], path_peaks_main); 
    #print(M_bed[:5])
    print('Empezando filtrar_L_rangos()')
    # Uso filtrar_L_rangos() para seleccionar peaks con genes y sitios de union cerca (pero no necesariamente igual de cerca de los genes)
    M_filtrado = filtrar_L_rangos(M_bed, sitios_cerca_genes=False); 
    print('Empezando secuencias_L_rangos()')
    # Uso secuencias_L_rango() para conseguir las secuencias de los rangos en formato .fasta
    M_fasta = secuencias_L_rangos(M_filtrado, genome_name='hg19', path_fasta=path_fasta_main, verbose=True); 
    print('Empezando guardar_fasta()')
    # Guardo M_fasta
    M_fasta = guardar_fasta(M_fasta, nom_out='testAnderson1500_sitios_cerca_genes_false_rango_variable', path_out=path_out_main); '''

    '''# Pruebo rango_medio()
    rango1 = (5,10); # len=6, center: 7,8
    rango2 = (6,12); # len=7, center: 9
    rango3 = (15,39); # len=25, center: 27
    rango4 = (12,37); # len=26, center: 24,25
    target1 = 2; 
    target2 = 10; 
    target3 = 3; 
    target4 = 15; 
    L_rango = [rango1, rango2, rango3, rango4]; 
    L_target = [target1, target2, target3, target4]; 
    _test_L_rango(L_rango, L_target); 
    #_test_rango(rango1, target1); '''

    '''# Pruebo filtrar_L_rangos()
    #M_filtrado = filtrar_L_rangos(M_bed, filtrar_su=True, con_su=True, filtrar_genes=True, con_genes=True, genes_filtrados=False, sitios_cerca_genes=True, L_col=[0,2,3,9]); 
    #print(M_filtrado[:10])
    _test_filtros(M_bed); '''

    return M_out


def _test_filtros(M_bed):
    # Prueba la funcion filtrar_L_rangos() para los distintos booleanos

    # Contador de pruebas
    cont_test = 0; 
    # Recorro cada uno de los booleanos
    for filtrar_su, filtrar_genes in [[True, True],[True, False],[False, True]]:
        if filtrar_su and filtrar_genes:
            for con_su in [True, False]:
                for con_genes in [True, False]:
                    if con_genes:
                        for genes_filtrados in [True, False]: 
                            for sitios_cerca_genes in [True, False]:
                                cont_test += 1; 
                                print('### Prueba numero ' + str(cont_test))
                                _test_un_filtro(M_bed, filtrar_su, con_su, filtrar_genes, con_genes, genes_filtrados, sitios_cerca_genes); 
                    else:
                        cont_test += 1; 
                        print('### Prueba numero ' + str(cont_test))
                        _test_un_filtro(M_bed, filtrar_su, con_su, filtrar_genes, con_genes, False, True); 
        elif filtrar_su:
            for con_su in [True, False]:
                cont_test += 1; 
                print('### Prueba numero ' + str(cont_test))
                _test_un_filtro(M_bed, filtrar_su, con_su, filtrar_genes, True, False, True); 
        elif filtrar_genes:
            for con_genes in [True, False]:
                if con_genes:
                    for genes_filtrados in [True, False]:
                        cont_test += 1; 
                        print('### Prueba numero ' + str(cont_test))
                        _test_un_filtro(M_bed, filtrar_su, True, filtrar_genes, con_genes, genes_filtrados, True); 
                else:
                    cont_test += 1; 
                    print('### Prueba numero ' + str(cont_test))
                    _test_un_filtro(M_bed, filtrar_su, True, filtrar_genes, con_genes, False, True); 
        else:
            print('ERROR? filtrar_su y filtrar_genes son False')
    return M_bed


def _test_un_filtro(M_bed, filtrar_su, con_su, filtrar_genes, con_genes, genes_filtrados, sitios_cerca_genes):
    # Prueba un solo caso de los booleanos para la funcion filtrar_L_rangos() y printea la lista de todas las clases de las distintas matrices

    # Texto para primer print
    str_ini = '>Probando\t?\nfiltrar_su:' + str(filtrar_su) + '\n\tcon_su:' + str(con_su) + '\nfiltrar_genes:' + str(filtrar_genes) + '\n\tcon_genes:'; 
    str_ini = str_ini + str(con_genes) + '\n\t\tgenes_filtrados:' + str(genes_filtrados) + '\n\t\tsitios_cerca_genes:' + str(sitios_cerca_genes); 
    print(str_ini)
    M_filtrado = filtrar_L_rangos(M_bed, filtrar_su=filtrar_su, con_su=con_su, filtrar_genes=filtrar_genes, con_genes=con_genes, 
                                  genes_filtrados=genes_filtrados, sitios_cerca_genes=sitios_cerca_genes, L_col=[0,2,3,9]); 
    # Chequeo de valores en columna class
    print('Valores unicos en columna class de M_bed:')
    print(_unique_col(M_bed, 9))
    print('Valores unicos en columna class de M_filtrado:')
    print(_unique_col(M_filtrado, 3))
    print()
    return M_filtrado


def _test_L_rango(L_rango, L_target):
    # Prueba varios rangos con varios target

    L_out = []; 
    for rango in L_rango:
        for target in L_target:
            L_out.append(_test_rango(rango,target)); 
            print()
    return L_out


def _test_rango(rango, target):
    # Prueba rango_medio() y printea output legible

    print('Probando rango ' + str(rango) + ' para largo target: ' + str(target))
    pos_ini, pos_end = rango_medio(rango[0], rango[1], target); 
    print('Rango medio: ' + str(pos_ini) + ', ' + str(pos_end))
    ret = [pos_ini, pos_end]; 
    return ret


def _unique_col(M_csv, col_sel):
    # Devuelve todas las variables (no repetidas) de la columna col_sel en la matriz M_csv

    # Inicializo la lista que se devuelve
    L_unique = []; 
    # Recorro M_csv
    for curr_row in M_csv:
        # Solo agrego las variables que no se encuentren ya en L_unique
        if not (curr_row[col_sel] in L_unique):
            L_unique.append(curr_row[col_sel]); 
    return L_unique


#################################### RESULTADOS ###################################

output_dump = []; 

if __name__=='__main__':
    output_dump.append(_main()); 

