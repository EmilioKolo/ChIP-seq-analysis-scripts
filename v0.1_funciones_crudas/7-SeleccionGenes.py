
import os
import time
import copy
from random import shuffle
import matplotlib.pyplot as plt
plt.style.use('classic')
from matplotlib import cycler
colors = cycler('color',
                ['#EE6666', '#3388BB', '#9988DD',
                 '#EECC55', '#88BB44', '#FFBBBB'])
import numpy as np
from pyensembl import EnsemblRelease
from Bio import Entrez, SeqIO
Entrez.email = 'ekolomenski@gmail.com';
Entrez.api_key = '9bccbd4949f226a5930868df36c211e8b408';
from referencias_funciones import ConsultaSecuencia, IDchr, extraer_columnas_csv, buscar_chr, buscar_en_seq_2dir_unificado, buscar_pos0, complemento_secuencia, _test_secuencia_correcta, histograma_con_lista

mm9 = EnsemblRelease(67, species='mouse');
hg19 = EnsemblRelease(102, species='human');
GRCm38 = EnsemblRelease(102, species='mouse');


'''
Contiene programas para seleccionar genes por su id de Ensembl

Funcion para buscar peaks cerca de lista de genes
    * buscar_peaks_en_lista_genes()

Funcion para buscar sitios definidos en peaks definidos
    * buscar_sitio_en_peak()

Funcion para seleccionar genes involucrados en desarrollo del corazon y hacer histogramas como en 5-HistogramasSitiosUnion.py
    * histogramas_por_sitio_union() (mismo nombre que en 5-HistogramasSitiosUnion.py)
'''

#################################### VARIABLES ####################################

nombres_archivos_chip_hg19 = ['Anderson2018-GSE89457consensus'];
nombres_archivos_chip_mm9 = ['Dupays2015', 'He2011', 'vandenBoogard2012'];
nombres_archivos_chip_mm9_extra = ['Dupays_GSM1087143_nkx_s1e14_peaks', 'Dupays_GSM1087144_nkx_s4e14_peaks'];
path_peaks_chip = '..\\0-Fuentes\\Papers ChIP-seq\\';
ext_chip = '.bed';

L_human = ['ENSG00000175206', 'ENSG00000265107', 'ENSG00000081189', 'ENSG00000171476', 
'ENSG00000141052', 'ENSG00000163485', 'ENSG00000179218', 'ENSG00000183023', 'ENSG00000117298', 
'ENSG00000164093', 'ENSG00000073146', 'ENSG00000141448'];
L_mouse = ['ENSMUSG00000041616', 'ENSMUSG00000057123', 'ENSMUSG00000005583', 'ENSMUSG00000059325', 
'ENSMUSG00000020542', 'ENSMUSG00000042429', 'ENSMUSG00000003814', 'ENSMUSG00000054640', 'ENSMUSG00000057530', 
'ENSMUSG00000028023', 'ENSMUSG00000015365', 'ENSMUSG00000005836'];

sitios_union_confirmados = ['GCAAGTG', 'GGAAGTG', 'TAAGTG', 'TCAAGTG'];
sitios_union_confirmados_v2 = ['GCAAGTG', 'GGAAGTG', 'GAAAGTG', 'ATAAGTG', 'GTAAGTG', 'CTAAGTG',
                               'TCAAGTG', 'TGAAGTG', 'TAAAGTG', 'TTAAGTG'];

L_peaks_cercanos = [['chr4', 147374182, 147374942], ['chr5', 77558912, 77560118], ['chr8', 87370330, 87371143]];
L_peaks_lejanos = [['chr17', 82136742, 82138153],['chr11', 65117066, 65117808]];
L_peaks_improbables = [['chr1', 136045697, 136046303],['chr15', 88850236, 88851312], ['chr3', 96866654, 96867281],
                       ['chr11', 65018536, 65019324], ['chr11', 65022712, 65023573], ['chr18', 11069363, 11070039]];

path_genes_cardio = '.\\7_ResultadosPaperNature\\';
L_arch_genes_cardio = ['Ensembl_ID_het_cardio_phenotype.txt', 'Ensembl_ID_hom_cardio_phenotype.txt', 'Ensembl_ID_hem_cardio_phenotype.txt'];

archivo_sitios_promotor = 'SitiosDeUnionGenoma';
archivo_sitios_lejanos = 'SitiosDeUnionLejanos';
path_sitios_genoma = '.\\sitios de union genoma\\';

L_bins_lejanos = [-11502,-11400,-11300,-11200,-11100,-11000,-10900,-10800,-10700,-10600,-10500,-10400,-10300,-10200,-10100,-10000];
L_bins_promotor = [-1502,-1400,-1300,-1200,-1100,-1000,-900,-800,-700,-600,-500,-400,-300,-200,-100,0];

curr_path = os.path.dirname(__file__);

# Variables main()
genoma_usado = mm9;
nombre_genoma_usado = 'mm9';
L_genes_main = L_mouse;
L_nombres_arch = nombres_archivos_chip_mm9_extra;
nombre_arch_usado = nombres_archivos_chip_mm9_extra[0];
sitios_union_main = sitios_union_confirmados_v2;

file_path_sitios_promotor = os.path.join(curr_path, str(path_sitios_genoma) + str(archivo_sitios_promotor));
file_path_sitios_lejanos = os.path.join(curr_path, str(path_sitios_genoma) + str(archivo_sitios_lejanos));

L_bins_main = L_bins_promotor;
nombre_arch_hist = archivo_sitios_promotor;
dir_arch_hist = file_path_sitios_promotor;
m = 0;

#################################### FUNCIONES ####################################


def abrir_resultados_bed(dir_arch, nom_arch, sep_arch='\t'):
    # Abre el archivo .bed en dir_arch y devuelve las primeras 3 columnas como matriz

    # Inicializo la matriz a devolver
    M_bed = [];
    
    # Abro el archivo
    with open(dir_arch, 'r') as F:
        print('Archivo ' + str(nom_arch) + '.bed abierto.')
        for curr_line in F:
            L_curr_line = curr_line.rstrip().split(sep_arch);
            M_bed.append(L_curr_line[:3]);
    return M_bed


def abrir_lista_ensembl_id(dir_arch, nom_arch):
    # Funcion que abre un .txt con un ensembl ID en cada fila y los devuelve como una lista de strings

    # Inicializo la lista que se devuelve
    L_out = [];

    # Abro el archivo
    with open(dir_arch, 'r') as F:
        print('Archivo ' + str(nom_arch) + '.bed abierto.')
        for curr_line in F:
            curr_id = str(curr_line.rstrip());
            L_out.append(curr_id);
    return L_out


def buscar_peaks_cerca_pos(pos_0, b_forw, chr_n, M_peaks, long_dist=5000):
    # Busca el peak mas cercano a pos_0 en la direccion rio arriba determinada por b_forw
    # M_peaks extraido de archivos .bed en formato [chr_n, pos_ini, pos_end]
    # recursive_b es para display en casos que se encuentre un peak a distancia mayor a long_dist

    # Inicializo la lista que se devuelve
    peak_out = ['chr',-1,-1];
    
    # Defino M_peaks_chr para registrar solo los peaks con el mismo chr_n
    M_peaks_chr = [];
    # Recorro M_peaks para seleccionar peaks en mismo chr_n
    for j in range(len(M_peaks)):
        curr_peak = M_peaks[j];
        if curr_peak[0] == chr_n:
            # Solo agrego los ultimos dos valores
            M_peaks_chr.append([int(curr_peak[1]),int(curr_peak[2])]);
    
    # Reviso que haya picos
    if len(M_peaks_chr)==0:
        print('No se encontraron peaks con chr_n = ' + str(chr_n) + '; se devuelve lista vacia')
        peak_out = [];

    # Ordeno M_peaks_chr
    M_peaks_chr = sorted(M_peaks_chr, key=lambda x: x[0]);
    
    # Recorro cada pico en M_peaks_chr
    k = 0;
    peak_encontrado = False;
    while (not peak_encontrado) and k < len(M_peaks_chr):
        # Si el gen es forward, se busca que el punto final del peak pase pos_0
        if b_forw:
            # Cuando el punto final de un peak pasa a pos_0, se elige uno de los dos peaks (k o k-1)
            if M_peaks_chr[k][1] > pos_0:
                peak_encontrado = True;
                # Si la posicion inicial del peak esta pasando pos0, se devuelve el peak anterior (k-1)
                if M_peaks_chr[k][0] > pos_0:
                    if k > 0:
                        peak_out = [str(chr_n)] + M_peaks_chr[k-1];
                        # Agrego la distancia de pos0 al peak
                        peak_out.append(pos_0-M_peaks_chr[k-1][1]);
                        if peak_out[-1] > long_dist:
                            print('Distancia muy grande al peak mas cercano: ' + str(peak_out[-1]))
                            print('Peak del otro lado: ' + str([str(chr_n)] + M_peaks_chr[k]))
                    # Excepto que sea el primer peak del cromosoma, en cuyo caso se devuelve M_peaks_chr[0]
                    else:
                        print('ERROR: El primer peak en el cromosoma esta pasando pos0. Se devuelve el primer peak.')
                        peak_out = [str(chr_n)] + M_peaks_chr[0];
                        # Agrego la distancia de pos0 al peak
                        peak_out.append(pos_0-M_peaks_chr[0][1]);
                # Si la posicion inicial del peak es pos0 o menor, el gen se encuentra dentro del peak
                else:
                    # Se devuelve directamente el peak k
                    peak_out = [str(chr_n)] + M_peaks_chr[k];
                    # Agrego la distancia de pos0 al peak
                    peak_out.append(0);
        # Si el gen es reverse, se busca que el punto inicial del peak este antes de pos_0 (y se recorre al reves)
        else:
            rev_k = len(M_peaks_chr)-1-k;
            # Cuando el punto inicial de un peak esta antes de pos_0, se elige uno de los dos peaks (rev_k o rev_k+1)
            if M_peaks_chr[rev_k][0] < pos_0:
                peak_encontrado = True;
                # Si la posicion final del peak esta antes de pos0, se devuelve el peak siguiente (rev_k+1)
                if M_peaks_chr[rev_k][1] < pos_0:
                    if rev_k < len(M_peaks_chr)-1:
                        peak_out = [str(chr_n)] + M_peaks_chr[rev_k+1];
                        # Agrego la distancia de pos0 al peak
                        peak_out.append(M_peaks_chr[rev_k+1][0]-pos_0);
                        # Reviso que la distancia al peak no sea mayor que long_dist
                        if peak_out[-1] > long_dist:
                            print('Distancia muy grande al peak mas cercano: ' + str(peak_out[-1]))
                            print('Peak del otro lado: ' + str([str(chr_n)] + M_peaks_chr[rev_k]))
                    # Excepto que sea el ultimo peak del cromosoma, en cuyo caso se devuelve M_peaks_chr[-1]
                    else:
                        print('ERROR: El ultimo peak en el cromosoma esta antes de pos0. Se devuelve el ultimo peak.')
                        peak_out = [str(chr_n)] + M_peaks_chr[-1];
                        # Agrego la distancia de pos0 al peak
                        peak_out.append(M_peaks_chr[-1][0]-pos_0);
                # Si la posicion final del peak es pos0 o mayor, el gen se encuentra dentro del peak
                else:
                    # Se devuelve directamente el peak rev_k
                    peak_out = [str(chr_n)] + M_peaks_chr[rev_k];
                    # Agrego la distancia de pos0 al peak
                    peak_out.append(0);
        k = k + 1;
    return peak_out


def buscar_peaks_en_lista_genes(L_gene_id, genome, genome_name, path_bed, nombre_bed, sep_bed='\t'):
    # Busca peaks de ChIP-seq en archivo nombre_chip y busca el pico mas cercano a pos0 de cada gen en L_gene_id

    # Inicializo la matriz que se devuelve
    M_out = [];
    
    # Abro el archivo nombre_bed y guardo el contenido en M_bed
    arch_path_bed = os.path.join(curr_path, path_bed + nombre_bed + '.bed');
    M_bed = abrir_resultados_bed(arch_path_bed, nombre_bed, sep_arch=sep_bed);
    print()

    # Hago el analisis para cada gen
    for i in range(len(L_gene_id)):
        print('>>>Revisando peaks cerca de gen:')
        curr_gene = genome.gene_by_id(L_gene_id[i]);
        print(curr_gene)
        print()
        # Defino pos0, forward y chr
        pos0, forward = buscar_pos0(curr_gene);
        curr_chr = buscar_chr(curr_gene);

        print('>Calculando peak mas cercano')
        # Busco los peaks mas cerca de pos0
        L_out = buscar_peaks_cerca_pos(pos0, forward, curr_chr, M_bed);
        print('Peak mas cercano: ' + str(L_out))
        print()
        print()

        # Agrego el Ensembl id del gen a L_out
        L_out.append(L_gene_id[i]);
        # Guardo L_out en M_out
        M_out.append(L_out[:]);
    return M_out


def buscar_sitio_en_peak(L_sitios, L_peak, genome_name, verbose=False, prueba_seq=False):
    # Busca una lista de sitios en una secuencia delimitada por L_peak (en formato chr_n, pos_ini, pos_end)
    # Devuelve la posicion absoluta del sitio (pos_ini, pos_end)

    # Inicializo la lista que se devuelve
    L_pos_sitios = [];

    # Determino el cromosoma a buscar en base a L_peak[0]
    chr_id = IDchr(L_peak[0][3:], genome=genome_name);
    # Determino la secuencia del peak
    seq_peak = ConsultaSecuencia(chr_id,int(L_peak[1]),int(L_peak[2]));
    # Busco cada sitio
    for curr_sitio in L_sitios:
        if verbose:
            print('Buscando sitio ' + str(curr_sitio))
        # Busco los sitios con buscar_en_seq_2dir_unificado()
        L_sitios_encontrados = buscar_en_seq_2dir_unificado(curr_sitio, seq_peak);
        # Paso los resultados a formato que prefiero
        for sitio_encontrado in L_sitios_encontrados:
            if verbose:
                print('Sitio encontrado en ' + str(sitio_encontrado))
            # Defino sitio_ini y sitio_end
            if sitio_encontrado < 0:
                sitio_ini = int(L_peak[2])+sitio_encontrado+1;
                sitio_end = sitio_ini + len(curr_sitio)-1;
            else:
                sitio_ini = int(L_peak[1])+sitio_encontrado+1;
                sitio_end = sitio_ini + len(curr_sitio)-1;
            # Creo una lista para armar el formato (seq_sitio, pos_ini, pos_end)
            L_pos_sitio_encontrado = [str(curr_sitio), sitio_ini, sitio_end];

            ##### PRUEBA:
            if prueba_seq:
                _test_secuencia_correcta(str(curr_sitio), chr_id, sitio_ini, sitio_end, verbose=verbose);
            
            # Guardo el sitio en formato correcto en L_pos_curr_sitio
            L_pos_sitios.append(L_pos_sitio_encontrado[:]);
    return L_pos_sitios


def extraer_columnas_csv(dir_archivo, nom_archivo, L_col, L_ids=[], ext='.csv', sep_arch=';'):
    # Abre un archivo .csv y devuelve las columnas en L_col como matriz
    # Si L_col es una lista vacia, se devuelven todas las columnas
    # Selecciona solo las filas con id (en columna 0) presente en L_ids

    # Creo la matriz que se devuelve
    M_out = [];
    # Abro el archivo .csv
    with open(str(dir_archivo) + str(ext), 'r') as F:
        print('Archivo ' + str(nom_archivo) + str(ext) + ' abierto.')
        #################### REVISAR ESTO ####################
        # Parte provisional del script que genera un output para revisar
        # Es posible que tenga que sacarlo si vuelvo a necesitar este script
        with open('ListaGenesCardio.csv', 'w') as F_out:
            print('!!!')
            print('VER SI HAY QUE SACAR ESTO DEL SCRIPT: Archivo ListaGenesCardio.csv creado.')
            print('!!!')
        with open('ListaGenesCardio.csv', 'w') as F_out:
            cont_out = 0;
            # Reviso cada linea de F
            for curr_line in F:
                # Paso cada linea a formato de lista
                L = curr_line.rstrip().split(sep_arch);
                # Creo una lista para registrar las columnas correspondientes
                L_out = [];
                # Si hay numeros en L_col, paso solo esas columnas
                if len(L_col) > 0:
                    for col in L_col:
                        L_out.append(L[col]);
                # Si L_col es una lista vacia, se devuelven todas las columnas
                else:
                    L_out = L[:];

                ###print('L[0]: ' + str(L[0]))
                # Reviso si L[0] se encuentra en L_ids antes de cargar
                if L_ids == []:
                    # Si no di L_ids para filtrar, guardo L_out en M_out
                    M_out.append(L_out[:]);
                elif L[0] in L_ids:
                    ###print('Presente en L_ids!')
                    # Chanchada que hice para hacer output
                    cont_out = cont_out + 1;
                    r_out = str(cont_out);
                    for h in L_out:
                        r_out = r_out + ';' + h;
                    F_out.write(r_out + '\n');
                    # Guardo L_out en M_out
                    M_out.append(L_out[:]);
    return M_out


def histogramas_por_sitio_union(dir_arch, nom_arch, L_bins_hist, L_ensembl_id, rango_hist=(-1600,100), modificador=0):
    # Crea histogramas por sitio de union con los datos extraidos de un archivo .csv
    # Funciona como la funcion de 5-HistogramasSitiosUnion.py pero selecciona IDs de ensembl que esten en L_ensembl_id

    # Extraigo la info del csv
    # La seleccion por Ensembl ID se hace en la funcion parseo_columnas_csv
    M_out = extraer_columnas_csv(dir_arch, nom_arch, [], L_ensembl_id);

    # Cargo la lista de sitios de union (L_sitios) y la matriz de posiciones para cada sitio registrado (M_hist)
    # Tambien recopilo una lista para hacer el histograma general con todos los datos (L_hist)
    L_sitios, M_hist, L_hist = parseo_columnas_csv(M_out, modificador);

    # Defino el rango a usar para los histogramas
    rango_usado = (rango_hist[0]+modificador, rango_hist[1]+modificador);
    # Hago los histogramas por sitio de union
    for i in range(len(M_hist)):
        L_pos_sitios = M_hist[i];
        print('Viendo M_hist[' + str(i) + '], correspondiente a secuencia ' + str(L_sitios[i]))
        print('Largo: ' + str(len(L_pos_sitios)))
        L_pos_sitios, dict_hist_sitios = histograma_con_lista(L_pos_sitios, str(L_sitios[i]) + '_m_' + str(modificador),
                                                              rango_usado, L_bins=L_bins_hist, use_float=False);
        print(dict_hist_sitios)
        print()

    # Hago un histograma del total de los datos
    L_hist, dict_hist_total = histograma_con_lista(L_hist, 'total_m_' + str(modificador), rango_usado, L_bins=L_bins_hist, use_float=False);
    print('Histograma total')
    print('Largo: ' + str(len(L_hist)))
    print(dict_hist_total)
    print()
    return dict_hist_total


def main():
    # Funcion para probar funciones en ejecucion del archivo

    M_out = [];
    #for archivo in L_nombres_arch:
    #    N_out = buscar_peaks_en_lista_genes(L_genes_main, genoma_usado, nombre_genoma_usado, path_peaks_chip, nombre_arch_usado);
    #    M_out.append(N_out[:]);
    
    #M_out = buscar_peaks_en_lista_genes(L_genes_main, genoma_usado, nombre_genoma_usado, path_peaks_chip, nombre_arch_usado);

    print('>>> Viendo primer tanda, peaks cercanos')
    print()
    for p in L_peaks_cercanos:
        print('> Peak: ' + str(p))
        M_p = buscar_sitio_en_peak(sitios_union_main, p, 'mm9', verbose=True);
        print()
    print()

    print('>>> Viendo segunda tanda, peaks lejanos.')
    print()
    for p in L_peaks_lejanos:
        print('> Peak: ' + str(p))
        M_p = buscar_sitio_en_peak(sitios_union_main, p, 'mm9', verbose=True);
        print()
    print()

    print('>>> Viendo segunda tanda, peaks improbables.')
    print()
    for p in L_peaks_improbables:
        print('> Peak: ' + str(p))
        M_p = buscar_sitio_en_peak(sitios_union_main, p, 'mm9', verbose=True);
        print()
    print()
    #M_out = buscar_sitio_en_peak(sitios_union_main, ['chr4', 147374182, 147374942], 'mm9', verbose=True);
    return M_out


def main_hist():
    # Funcion para usar funciones de visualizacion de resultados

    # Extraigo los genes de la lista del paper de Nature
    L_ensembl_id = [];
    for archivo in L_arch_genes_cardio:
        dir_ensembl = os.path.join(curr_path, path_genes_cardio + archivo);
        curr_L_ensembl_id = abrir_lista_ensembl_id(dir_ensembl, archivo);
        L_ensembl_id = L_ensembl_id + curr_L_ensembl_id;
    # Creo los histogramas por sitio de union
    M_out = histogramas_por_sitio_union(dir_arch_hist, nombre_arch_hist, L_bins_main, L_ensembl_id, rango_hist=(-1600,100), modificador=m);
    
    return M_out


def main_libre():
    # Para probar cosas sueltas

    M_out = [];
    # Extraigo los genes de la lista del paper de Nature
    #L_ensembl_id = [];
    #for archivo in L_arch_genes_cardio:
    #    dir_ensembl = os.path.join(curr_path, path_genes_cardio + archivo);
    #    curr_L_ensembl_id = abrir_lista_ensembl_id(dir_ensembl, archivo);
    #    L_ensembl_id = L_ensembl_id + curr_L_ensembl_id;
    # Pruebo que los genes asociados a Nkx2-5 aparezcan en la lista de genes del paper de Nature
    #M_out = _prueba_genes_nature(L_genes_main, L_ensembl_id, genoma_usado);
    return M_out


def main_test():
    # Funcion para probar funciones por partes

    M_out = [];
    
    #M_peaks_test = [['chr1', 2000, 2500],['chr1', 3000, 3500],['chr1', 1000, 1500],['chr1', 4000, 5500],
    #                ['chr2', 1200, 1700],['chr2', 6200, 6700],['chr2', 5200, 5700],['chr2', 2200, 2700]];
    #print('Iniciando prueba de buscar_peaks_cerca_pos()')
    #M_out = buscar_peaks_cerca_pos(2780, True, 'chr2', M_peaks_test, long_dist=500);
    #print('Prueba terminada, M_out:')
    #print(M_out)

    #M_out = buscar_sitio_en_peak(sitios_union_main, ['chr4', 147374182, 147374942], 'mm9', verbose=True, prueba_seq=True);

    test_arch_ensembl = L_arch_genes_cardio[0];
    dir_ensembl = os.path.join(curr_path, path_genes_cardio + test_arch_ensembl);
    L_ensembl_id = abrir_lista_ensembl_id(dir_ensembl, test_arch_ensembl);
    M_out = _prueba_extraer_columnas_csv(file_path_sitios_promotor, archivo_sitios_promotor, L_ensembl_id);
    return M_out


def parseo_columnas_csv(M_csv, m):
    # Funcion que recibe una matriz M_csv recopilando sitios de union y posiciones relativas
    # Devuelve lista de sitios de union (L_sitios), matriz de posiciones para cada sitio registrado (M_hist)
    # Tambien devuelve una lista para hacer un histograma general con todos los datos (L_hist)
    
    # Matriz con listas de posiciones relativas para cada sitio en L_sitios
    M_hist = [];
    # Lista de los distintos sitios de union en M_csv
    L_sitios = [];
    # Lista para probar un histograma general con los datos arreglados
    L_hist = [];
    # Recorro M_csv
    for j in range(len(M_csv)):
        sitio = M_csv[j];
        # Cargo el sitio de union si no esta ya en L_sitios
        if not (sitio[1] in L_sitios):
            L_sitios.append(sitio[1]);
            # Tambien inicializo una lista en M_hist
            M_hist.append([]);
        # Defino la posicion del sitio en L_sitios para hacer el append en M_hist
        pos_sitio = L_sitios.index(sitio[1]);

        # Cargo la posicion del sitio en M_hist y en L_hist, sumando el modificador
        M_hist[pos_sitio].append(int(sitio[2])+int(m));
        L_hist.append(int(sitio[2])+int(m));
    return L_sitios, M_hist, L_hist


def _prueba_extraer_columnas_csv(dir_arch, nom_arch, L_ids):
    # Prueba de la funcion extraer_columnas_csv

    # Inicializo la matriz que se devuelve
    M_out = extraer_columnas_csv(dir_arch, nom_arch, [], L_ids);
    return M_out


def _prueba_genes_nature(L_genes, L_ids, genome):
    # Prueba que los genes en L_genes aparecen en la lista extraida del paper de Nature

    # Inicializo la matriz que se devuelve
    M_out = [];
    # Recorro cada uno de los genes y hago display
    for i in range(len(L_genes)):
        curr_ensembl_id = L_genes[i];
        print('Leyendo ensembl id: ' + str(curr_ensembl_id))
        curr_gene = genome.gene_by_id(L_genes[i]);
        print('Gen encontrado: ')
        print(curr_gene)
        if curr_ensembl_id in L_ids:
            print('Gen presente en lista de IDs de Ensembl!')
        elif curr_gene.gene_id in L_ids:
            print('ADVERTENCIA: Gen presente en lista de IDs de Ensembl!')
        else:
            print('TODO MAL: Gen NO ENCONTRADO en la lista de IDs de Ensembl.')
        print()
    return M_out


#################################### RESULTADOS ###################################

output_dump = '';

if __name__=='__main__':
    #output_dump = main();
    #output_dump = main_test();
    #output_dump = main_hist();
    output_dump = main_libre();


