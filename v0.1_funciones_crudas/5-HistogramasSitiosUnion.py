

import os
import copy
import time
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
from referencias_funciones import ConsultaSecuencia, IDchr, extraer_columnas_csv

mm9 = EnsemblRelease(67, species='mouse');
hg19 = EnsemblRelease(102, species='human');
GRCm38 = EnsemblRelease(102, species='mouse');

'''
Agarro los resultados de 5-SitiosDeUnionGenoma.py y hago un histograma de distancias

Incluye funciones para hacer histogramas por sitio de union

Incluye funcion para ver probabilidad de sitio de union en secuencia de largo n
'''


#################################### VARIABLES ####################################

# Corrida original

nombre_archivo = 'SitiosDeUnionGenoma';
columna = 2;

extension = '.csv';
path_archivo = '.\\sitios de union genoma\\';
curr_path = os.path.dirname(__file__);

file_path = os.path.join(curr_path, str(path_archivo) + str(nombre_archivo));

genome = mm9;
genome_name = 'mm9';

nombre_archivo_lejano = 'SitiosDeUnionLejanos';
file_path_lejano = os.path.join(curr_path, str(path_archivo) + str(nombre_archivo_lejano));

L_bins_lejano = [-11502,-11400,-11300,-11200,-11100,-11000,-10900,-10800,-10700,-10600,-10500,-10400,-10300,-10200,-10100,-10000];
L_bins_promotor = [-1502,-1400,-1300,-1200,-1100,-1000,-900,-800,-700,-600,-500,-400,-300,-200,-100,0];

# Corrida amboslados

nombre_archivo_amboslados = 'SitiosDeUnionGenoma_AmbosLados' + '_hg19';
nombre_archivo_pos = 'SitiosDeUnionLejanos_pos' + '_hg19';
nombre_archivo_neg = 'SitiosDeUnionLejanos_neg' + '_hg19';

file_path_amboslados = os.path.join(curr_path, str(path_archivo) + str(nombre_archivo_amboslados));
file_path_pos = os.path.join(curr_path, str(path_archivo) + str(nombre_archivo_pos));
file_path_neg = os.path.join(curr_path, str(path_archivo) + str(nombre_archivo_neg));

L_bins_amboslados = [-1001, -900, -800, -700, -600, -500, -400, -300, -200, -100, 0,
                     100, 200, 300, 400, 500, 600, 700, 800, 900, 1001];
L_bins_pos = [8999, 9200, 9400, 9600, 9800, 10000, 10200, 10400, 10600, 10800, 11001];
L_bins_neg = [-11001, -10800, -10600, -10400, -10200, -10000, -9800, -9600, -9400, -9200, -8999];

#################################### FUNCIONES ####################################


def extraer_columna_csv(dir_arch, nombre_archivo, n_col, ext='.csv', sep_arch=';'):
    # Abre un archivo .csv y devuelve la columna numero n_col como lista
    # DISTINTO DE FUNCION extraer_columnas_csv en referencias_funciones.py

    # Creo la lista que se devuelve
    L_col = [];
    # Abro el archivo .csv
    with open(str(dir_arch) + str(ext), 'r') as F:
        print('Archivo ' + str(nombre_archivo) + str(ext) + ' abierto.')
        # Reviso cada linea de F
        for curr_line in F:
            # Paso cada linea a formato de lista
            L = curr_line.rstrip().split(sep_arch);
            # Guardo el valor correspondiente a la columna n_col
            L_col.append(L[n_col]);
    return L_col


def histograma_con_lista(L, nombre_histograma, rango_hist, bins=10, L_bins=[], use_float=False, guardar=True):
    # Hace un histograma con una lista de numeros (L)
    # Se asegura que los elementos sean numeros

    # Creo un diccionario para funcionar como histograma
    dict_hist = {};
    # Creo la lista para tener los numeros pasados a int o float
    L_numbers = [];
    # Convierto todos los numeros a int o float de acuerdo a use_float
    for i in range(len(L)):
        if use_float:
            L_numbers.append(float(L[i]));
        else:
            L_numbers.append(int(L[i]));
    # Registro el rango de los numeros
    min_val = min(L_numbers);
    max_val = max(L_numbers);

    # Veo si se paso una lista de bins con mas de un valor
    if len(L_bins) < 2:
        # Si no hay, defino una
        L_bins_usada = [];
        ########## ESTO TODAVIA NO ANDA
    else:
        L_bins_usada = sorted(L_bins);
    
    for i in range(len(L_bins_usada)-1):
        dict_hist[L_bins[i]] = 0;
        for j in range(len(L_numbers)):
            if L_numbers[j] > L_bins[i] and L_numbers[j] <= L_bins[i+1]:
                dict_hist[L_bins[i]] += 1;

    # Creo un histograma
    n, bins, patches = plt.hist(x=L, bins=L_bins_usada, color='#040499', rwidth=0.9, alpha=0.8);
    plt.grid(axis='y', alpha=0.6);
    plt.xlabel('Distancia al +1');
    plt.ylabel('Frecuencia');
    plt.title(nombre_histograma);
    maxfreq = n.max()
    plt.ylim(ymax=np.ceil(maxfreq / 100) * 100 if maxfreq % 100 else maxfreq + 100);
    plt.xlim(xmin=rango_hist[0], xmax=rango_hist[1]);
    if guardar:
        plt.savefig('Histograma_' + str(nombre_histograma) + '.png');
    else:
        plt.show();
    plt.close();
    return L_numbers, dict_hist


def histogramas_por_sitio_union(dir_arch, nom_arch, L_bins_hist, rango_hist=(-1600,100), modificador=0):
    # Crea histogramas por sitio de union con los datos extraidos de un archivo .csv

    # Extraigo la info del csv
    M_out = extraer_columnas_csv(dir_arch, nom_arch, []);

    # Cargo la lista de sitios de union (L_sitios) y la matriz de posiciones para cada sitio registrado (M_hist)
    # Tambien recopilo una lista para hacer el histograma general con todos los datos (L_hist)
    L_sitios, M_hist, L_hist = parseo_columnas_csv(M_out, modificador)

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
    # Para probar el histograma
    L_sitios = extraer_columna_csv(file_path, nombre_archivo, columna, ext=extension)
    L_sitios, dict_hist_sitios = histograma_con_lista(L_sitios, L_bins=[-1502,-1250,-1000,-750,-500,-250,0], use_float=False)
    print(dict_hist_sitios)
    return None


def main2():
    # Para probar probabilidad de que los sitios de union aparezcan por azar
    for i in [5,6,7]:
        print('>Probabilidades para sitio de largo ' + str(i))
        for j in [250, 500, 1500, 3000, 5000]:
            print('En secuencia de largo ' + str(j) + ': ')
            print(prob_site_in_seq(i, j))
        print()
    return None


def main3():
    # Para probar histogramas por sitio de union
    # Hay un pedazo comentado que corresponde a un chequeo si los datos estan mal cargados

    # Extraigo la info del csv
    M_out = extraer_columnas_csv(file_path, nombre_archivo, []);

    # Lista de bins para los diccionarios que definen histogramas
    L_bins_hist = [-1502,-1400,-1300,-1200,-1100,-1000,-900,-800,-700,-600,-500,-400,-300,-200,-100,0];
    #L_bins_hist = [-1502,-1250,-1000,-750,-500,-250,0];
    
    # Cargo la lista de sitios de union (L_sitios) y la matriz de posiciones para cada sitio registrado (M_hist)
    # Tambien recopilo una lista para hacer el histograma general con todos los datos (L_hist)
    L_sitios, M_hist, L_hist = parseo_columnas_csv(M_out, 0);

    # Hago los histogramas por sitio de union
    for i in range(len(M_hist)):
        L_pos_sitios = M_hist[i];
        print('Viendo M_hist[' + str(i) + '], correspondiente a secuencia ' + str(L_sitios[i]))
        #print(L_pos_sitios[:20])
        print('Largo: ' + str(len(L_pos_sitios)))
        L_pos_sitios, dict_hist_sitios = histograma_con_lista(L_pos_sitios, str(L_sitios[i]) + '_promotor', (-1600, 100), L_bins=L_bins_hist, use_float=False);
        print(dict_hist_sitios)
        print()

    # Pruebo un histograma del total de los datos
    L_hist, dict_hist_total = histograma_con_lista(L_hist, 'total_promotor', (-1600, 100), L_bins=L_bins_hist, use_float=False);
    print('Histograma total')
    print('Largo: ' + str(len(L_hist)))
    print(dict_hist_total)
    print()
    return None


def main4():
    # Para probar histogramas por sitio de union EN REGIONES LEJANAS

    # Extraigo la info del csv
    M_out = extraer_columnas_csv(file_path_lejano, nombre_archivo_lejano, []);

    # Lista de bins para los diccionarios que definen histogramas
    L_bins_hist = [-11502,-11400,-11300,-11200,-11100,-11000,-10900,-10800,
                   -10700,-10600,-10500,-10400,-10300,-10200,-10100,-10000];
    #L_bins_hist = [-11502,-11250,-11000,-10750,-10500,-10250,-10000];
    
    # Cargo la lista de sitios de union (L_sitios) y la matriz de posiciones para cada sitio registrado (M_hist)
    # Tambien recopilo una lista para hacer el histograma general con todos los datos (L_hist)
    L_sitios, M_hist, L_hist = parseo_columnas_csv(M_out, -10000);

    # Hago los histogramas por sitio de union
    for i in range(len(M_hist)):
        L_pos_sitios = M_hist[i];
        print('Viendo M_hist[' + str(i) + '], correspondiente a secuencia ' + str(L_sitios[i]))
        #print(L_pos_sitios[:20])
        print('Largo: ' + str(len(L_pos_sitios)))
        L_pos_sitios, dict_hist_sitios = histograma_con_lista(L_pos_sitios, str(L_sitios[i]) + '_lejano', (-11600, -9900), L_bins=L_bins_hist, use_float=False);
        print(dict_hist_sitios)
        print()

    # Pruebo un histograma del total de los datos
    L_hist, dict_hist_total = histograma_con_lista(L_hist, 'total_lejano', (-11600, -9900), L_bins=L_bins_hist, use_float=False);
    print('Histograma total')
    print('Largo: ' + str(len(L_hist)))
    print(dict_hist_total)
    print()
    return None


def main5():
    # Para histogramas por sitio de union en humano

    D_out_prom = histogramas_por_sitio_union(file_path + '_hg19', nombre_archivo + '_hg19', L_bins_promotor);

    D_out_lejos = histogramas_por_sitio_union(file_path_lejano + '_hg19', nombre_archivo_lejano + '_hg19', L_bins_lejano, modificador=-10000);
    
    return None


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

'''
def histogramas_por_sitio_union(dir_arch, nom_arch, L_bins_hist, rango_hist=(-1600,100), modificador=0):
    # Crea histogramas por sitio de union con los datos extraidos de un archivo .csv

    # Extraigo la info del csv
    M_out = extraer_columnas_csv(dir_arch, nom_arch, []);

    # Cargo la lista de sitios de union (L_sitios) y la matriz de posiciones para cada sitio registrado (M_hist)
    # Tambien recopilo una lista para hacer el histograma general con todos los datos (L_hist)
    L_sitios, M_hist, L_hist = parseo_columnas_csv(M_out, modificador)

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
'''

def pipeline_histogramas(dir_arch, nom_arch, L_bins, rango_hist, output_mod):
    # Crea histogramas por sitio de union y totales con los datos extraidos del archivo .csv

    # Extraigo la info del .csv
    M_csv = extraer_columnas_csv(dir_arch, nom_arch, []);

    # Cargo la lista de sitios de union (L_sitios) y la matriz de posiciones para cada sitio registrado (M_hist)
    # Tambien recopilo una lista para hacer el histograma general con todos los datos (L_hist)
    L_sitios, M_hist, L_hist = parseo_columnas_csv(M_csv, 0);

    # Defino el rango a usar para los histogramas
    rango_usado = (rango_hist[0], rango_hist[1]);

    # Hago los histogramas por sitio de union
    for i in range(len(M_hist)):
        L_pos_sitios = M_hist[i];
        print('Viendo M_hist[' + str(i) + '], correspondiente a secuencia ' + str(L_sitios[i]))
        print('Largo: ' + str(len(L_pos_sitios)))
        L_pos_sitios, dict_hist_sitios = histograma_con_lista(L_pos_sitios, str(L_sitios[i]) + '_' + str(output_mod),
                                                              rango_usado, L_bins=L_bins, use_float=False);
        print(dict_hist_sitios)
        print()

    # Hago un histograma del total de los datos
    L_hist, dict_hist_total = histograma_con_lista(L_hist, 'total_' + str(output_mod), rango_usado,
                                                   L_bins=L_bins, use_float=False);
    print('Histograma total')
    print('Largo: ' + str(len(L_hist)))
    print(dict_hist_total)
    print()
    return dict_hist_total


def prob_site_in_seq(largo_sitio, largo_seq):
    # Devuelve la probabilidad de que una secuencia de adn largo_sitio aparezca en una secuencia largo_seq
    # Recibe los largos de ambas secuencias
    # largo_sitio tiene que ser menor o igual que largo_seq
    # ES UNA APROXIMACION, la probabilidad real es un poco mas alta

    # Probabilidad de que el sitio aparezca en una secuencia de su largo
    prob_sitio = 1/(4**largo_sitio);
    # Probabilidad de que el sitio no aparezca en una secuencia de su largo
    prob_not_sitio = 1 - prob_sitio;
    # Cantidad de posiciones posibles del sitio en la secuencia largo_seq
    #### SE PUEDE DUPLICAR PARA CONSIDERAR ANTISENSE?
    posiciones_sitio = largo_seq - largo_sitio + 1
    # Probabilidad de que el sitio no aparezca en ninguna de las posiciones
    prob_not_sitio_posiciones = prob_not_sitio**posiciones_sitio
    
    return 1 - prob_not_sitio_posiciones


def _main_amboslados():
    #

    # Inicializo la variable que se devuelve
    M_out = [];

    # Corro el pipeline para resultados AmbosLados
    rango_amboslados = (min(L_bins_amboslados), max(L_bins_amboslados));
    M_out = pipeline_histogramas(file_path_amboslados, nombre_archivo_amboslados,
                                 L_bins_amboslados, rango_amboslados, 'AmbosLados');
    return M_out


def _main_pos():
    #

    # Inicializo la variable que se devuelve
    M_out = [];

    # Corro el pipeline para resultados positivos
    rango_pos = (min(L_bins_pos), max(L_bins_pos));
    M_out = pipeline_histogramas(file_path_pos, nombre_archivo_pos,
                                 L_bins_pos, rango_pos, 'pos');
    return M_out


def _main_neg():
    #

    # Inicializo la variable que se devuelve
    M_out = [];

    # Corro el pipeline para resultados negativos
    rango_neg = (min(L_bins_neg), max(L_bins_neg));
    M_out = pipeline_histogramas(file_path_neg, nombre_archivo_neg,
                                 L_bins_neg, rango_neg, 'neg');
    return M_out


#################################### RESULTADOS ###################################

output_dump = [];

if __name__=='__main__':
    #main();
    #main2();
    #main3();
    #main4();
    #main5();
    output_dump.append(_main_amboslados());
    output_dump.append(_main_pos());
    output_dump.append(_main_neg());


'''
####### RESULTADOS HISTOGRAMAS -11.5k/-10k ("lejano") y -1500/0 ("promotores")


>TOTAL

Lejano
Largo: 81986
{-11502: 5502, -11400: 5484, -11300: 5485, -11200: 5381, -11100: 5434, -11000: 5438, -10900: 5524, -10800: 5503, -10700: 5403, -10600: 5525, -10500: 5462, -10400: 5338, -10300: 5538, -10200: 5556, -10100: 5056}

Promotores
Largo: 80449
{-1502: 5657, -1400: 5485, -1300: 5457, -1200: 5400, -1100: 5429, -1000: 5468, -900: 5472, -800: 5513, -700: 5587, -600: 5298, -500: 5299, -400: 5226, -300: 5040, -200: 5089, -100: 5029}


>Secuencia GGAAGTG

Lejano
Largo: 11002
{-11502: 744, -11400: 756, -11300: 714, -11200: 680, -11100: 677, -11000: 745, -10900: 750, -10800: 731, -10700: 747, -10600: 749, -10500: 755, -10400: 766, -10300: 726, -10200: 761, -10100: 680}

Promotores
Largo: 12301
{-1502: 785, -1400: 758, -1300: 720, -1200: 771, -1100: 752, -1000: 737, -900: 750, -800: 774, -700: 735, -600: 755, -500: 746, -400: 759, -300: 795, -200: 940, -100: 1524}


>Secuencia ATAAGTG

Lejano
Largo: 7485
{-11502: 500, -11400: 483, -11300: 482, -11200: 473, -11100: 508, -11000: 462, -10900: 490, -10800: 482, -10700: 467, -10600: 512, -10500: 482, -10400: 486, -10300: 508, -10200: 524, -10100: 458}

Promotores
Largo: 6405
{-1502: 498, -1400: 407, -1300: 493, -1200: 464, -1100: 440, -1000: 466, -900: 437, -800: 505, -700: 422, -600: 415, -500: 377, -400: 408, -300: 374, -200: 421, -100: 278}


>Secuencia TTAAGTG

Lejano
Largo: 7365
{-11502: 513, -11400: 488, -11300: 494, -11200: 473, -11100: 483, -11000: 507, -10900: 474, -10800: 514, -10700: 508, -10600: 472, -10500: 476, -10400: 467, -10300: 506, -10200: 488, -10100: 460}

Promotores
Largo: 7325
{-1502: 519, -1400: 536, -1300: 464, -1200: 464, -1100: 554, -1000: 520, -900: 561, -800: 558, -700: 506, -600: 511, -500: 505, -400: 495, -300: 453, -200: 406, -100: 273}


>Secuencia TCAAGTG

Lejano
Largo: 8101
{-11502: 532, -11400: 557, -11300: 547, -11200: 520, -11100: 547, -11000: 554, -10900: 512, -10800: 548, -10700: 508, -10600: 548, -10500: 537, -10400: 541, -10300: 530, -10200: 555, -10100: 565}

Promotores
Largo: 7717
{-1502: 588, -1400: 560, -1300: 525, -1200: 502, -1100: 506, -1000: 521, -900: 528, -800: 508, -700: 573, -600: 490, -500: 552, -400: 507, -300: 476, -200: 512, -100: 369}


>Secuencia TAAAGTG

Lejano
Largo: 8274
{-11502: 524, -11400: 547, -11300: 575, -11200: 628, -11100: 563, -11000: 542, -10900: 599, -10800: 560, -10700: 544, -10600: 560, -10500: 539, -10400: 532, -10300: 566, -10200: 498, -10100: 476}

Promotores
Largo: 8160
{-1502: 602, -1400: 574, -1300: 594, -1200: 557, -1100: 544, -1000: 590, -900: 533, -800: 578, -700: 575, -600: 569, -500: 583, -400: 522, -300: 490, -200: 478, -100: 371}


>Secuencia TGAAGTG

Lejano
Largo: 8989
{-11502: 623, -11400: 599, -11300: 581, -11200: 609, -11100: 602, -11000: 614, -10900: 561, -10800: 648, -10700: 594, -10600: 632, -10500: 580, -10400: 582, -10300: 605, -10200: 562, -10100: 534}

Promotores
Largo: 8663
{-1502: 639, -1400: 611, -1300: 584, -1200: 611, -1100: 590, -1000: 532, -900: 630, -800: 590, -700: 608, -600: 538, -500: 551, -400: 584, -300: 575, -200: 512, -100: 508}


>Secuencia GTAAGTG

Lejano
Largo: 5883
{-11502: 360, -11400: 400, -11300: 454, -11200: 389, -11100: 405, -11000: 324, -10900: 460, -10800: 378, -10700: 369, -10600: 383, -10500: 401, -10400: 357, -10300: 441, -10200: 414, -10100: 348}

Promotores
Largo: 5796
{-1502: 397, -1400: 406, -1300: 393, -1200: 357, -1100: 345, -1000: 495, -900: 430, -800: 387, -700: 379, -600: 426, -500: 374, -400: 325, -300: 349, -200: 300, -100: 433}


>Secuencia GAAAGTG

Lejano
Largo: 10639
{-11502: 692, -11400: 724, -11300: 694, -11200: 649, -11100: 700, -11000: 733, -10900: 736, -10800: 715, -10700: 715, -10600: 705, -10500: 740, -10400: 691, -10300: 720, -10200: 744, -10100: 660}

Promotores
Largo: 10182
{-1502: 685, -1400: 653, -1300: 679, -1200: 711, -1100: 698, -1000: 665, -900: 670, -800: 711, -700: 769, -600: 666, -500: 665, -400: 732, -300: 656, -200: 686, -100: 536}


>Secuencia GCAAGTG

Lejano
Largo: 7427
{-11502: 523, -11400: 480, -11300: 501, -11200: 523, -11100: 504, -11000: 500, -10900: 478, -10800: 478, -10700: 478, -10600: 527, -10500: 471, -10400: 475, -10300: 496, -10200: 505, -10100: 467}

Promotores
Largo: 7552
{-1502: 517, -1400: 531, -1300: 574, -1200: 511, -1100: 511, -1000: 508, -900: 494, -800: 512, -700: 537, -600: 482, -500: 515, -400: 491, -300: 462, -200: 489, -100: 418}


>Secuencia CTAAGTG

Lejano
Largo: 6821
{-11502: 491, -11400: 450, -11300: 443, -11200: 437, -11100: 445, -11000: 457, -10900: 464, -10800: 449, -10700: 473, -10600: 437, -10500: 481, -10400: 441, -10300: 440, -10200: 505, -10100: 408}

Promotores
Largo: 6348
{-1502: 427, -1400: 449, -1300: 431, -1200: 452, -1100: 489, -1000: 434, -900: 439, -800: 390, -700: 483, -600: 446, -500: 431, -400: 403, -300: 410, -200: 345, -100: 319}
'''

#############################################################################################################################################
########################################################## RESULTADOS SIN PROCESAR ##########################################################
#############################################################################################################################################

'''
Archivo SitiosDeUnionGenoma.csv abierto.
Viendo M_hist[0], correspondiente a secuencia GGAAGTG
Largo: 12301
{-1502: 785, -1400: 758, -1300: 720, -1200: 771, -1100: 752, -1000: 737, -900: 750, -800: 774, -700: 735, -600: 755, -500: 746, -400: 759, -300: 795, -200: 940, -100: 1524}

Viendo M_hist[1], correspondiente a secuencia TGAAGTG
Largo: 8663
{-1502: 639, -1400: 611, -1300: 584, -1200: 611, -1100: 590, -1000: 532, -900: 630, -800: 590, -700: 608, -600: 538, -500: 551, -400: 584, -300: 575, -200: 512, -100: 508}

Viendo M_hist[2], correspondiente a secuencia GTAAGTG
Largo: 5796
{-1502: 397, -1400: 406, -1300: 393, -1200: 357, -1100: 345, -1000: 495, -900: 430, -800: 387, -700: 379, -600: 426, -500: 374, -400: 325, -300: 349, -200: 300, -100: 433}

Viendo M_hist[3], correspondiente a secuencia GAAAGTG
Largo: 10182
{-1502: 685, -1400: 653, -1300: 679, -1200: 711, -1100: 698, -1000: 665, -900: 670, -800: 711, -700: 769, -600: 666, -500: 665, -400: 732, -300: 656, -200: 686, -100: 536}

Viendo M_hist[4], correspondiente a secuencia ATAAGTG
Largo: 6405
{-1502: 498, -1400: 407, -1300: 493, -1200: 464, -1100: 440, -1000: 466, -900: 437, -800: 505, -700: 422, -600: 415, -500: 377, -400: 408, -300: 374, -200: 421, -100: 278}

Viendo M_hist[5], correspondiente a secuencia TCAAGTG
Largo: 7717
{-1502: 588, -1400: 560, -1300: 525, -1200: 502, -1100: 506, -1000: 521, -900: 528, -800: 508, -700: 573, -600: 490, -500: 552, -400: 507, -300: 476, -200: 512, -100: 369}

Viendo M_hist[6], correspondiente a secuencia TAAAGTG
Largo: 8160
{-1502: 602, -1400: 574, -1300: 594, -1200: 557, -1100: 544, -1000: 590, -900: 533, -800: 578, -700: 575, -600: 569, -500: 583, -400: 522, -300: 490, -200: 478, -100: 371}

Viendo M_hist[7], correspondiente a secuencia GCAAGTG
Largo: 7552
{-1502: 517, -1400: 531, -1300: 574, -1200: 511, -1100: 511, -1000: 508, -900: 494, -800: 512, -700: 537, -600: 482, -500: 515, -400: 491, -300: 462, -200: 489, -100: 418}

Viendo M_hist[8], correspondiente a secuencia TTAAGTG
Largo: 7325
{-1502: 519, -1400: 536, -1300: 464, -1200: 464, -1100: 554, -1000: 520, -900: 561, -800: 558, -700: 506, -600: 511, -500: 505, -400: 495, -300: 453, -200: 406, -100: 273}

Viendo M_hist[9], correspondiente a secuencia CTAAGTG
Largo: 6348
{-1502: 427, -1400: 449, -1300: 431, -1200: 452, -1100: 489, -1000: 434, -900: 439, -800: 390, -700: 483, -600: 446, -500: 431, -400: 403, -300: 410, -200: 345, -100: 319}

Histograma total
Largo: 80449
{-1502: 5657, -1400: 5485, -1300: 5457, -1200: 5400, -1100: 5429, -1000: 5468, -900: 5472, -800: 5513, -700: 5587, -600: 5298, -500: 5299, -400: 5226, -300: 5040, -200: 5089, -100: 5029}


Archivo SitiosDeUnionLejanos.csv abierto.
Viendo M_hist[0], correspondiente a secuencia TCAAGTG
Largo: 8101
{-11502: 532, -11400: 557, -11300: 547, -11200: 520, -11100: 547, -11000: 554, -10900: 512, -10800: 548, -10700: 508, -10600: 548, -10500: 537, -10400: 541, -10300: 530, -10200: 555, -10100: 565}

Viendo M_hist[1], correspondiente a secuencia ATAAGTG
Largo: 7485
{-11502: 500, -11400: 483, -11300: 482, -11200: 473, -11100: 508, -11000: 462, -10900: 490, -10800: 482, -10700: 467, -10600: 512, -10500: 482, -10400: 486, -10300: 508, -10200: 524, -10100: 458}

Viendo M_hist[2], correspondiente a secuencia TAAAGTG
Largo: 8274
{-11502: 524, -11400: 547, -11300: 575, -11200: 628, -11100: 563, -11000: 542, -10900: 599, -10800: 560, -10700: 544, -10600: 560, -10500: 539, -10400: 532, -10300: 566, -10200: 498, -10100: 476}

Viendo M_hist[3], correspondiente a secuencia CTAAGTG
Largo: 6821
{-11502: 491, -11400: 450, -11300: 443, -11200: 437, -11100: 445, -11000: 457, -10900: 464, -10800: 449, -10700: 473, -10600: 437, -10500: 481, -10400: 441, -10300: 440, -10200: 505, -10100: 408}

Viendo M_hist[4], correspondiente a secuencia GCAAGTG
Largo: 7427
{-11502: 523, -11400: 480, -11300: 501, -11200: 523, -11100: 504, -11000: 500, -10900: 478, -10800: 478, -10700: 478, -10600: 527, -10500: 471, -10400: 475, -10300: 496, -10200: 505, -10100: 467}

Viendo M_hist[5], correspondiente a secuencia GGAAGTG
Largo: 11002
{-11502: 744, -11400: 756, -11300: 714, -11200: 680, -11100: 677, -11000: 745, -10900: 750, -10800: 731, -10700: 747, -10600: 749, -10500: 755, -10400: 766, -10300: 726, -10200: 761, -10100: 680}

Viendo M_hist[6], correspondiente a secuencia GAAAGTG
Largo: 10639
{-11502: 692, -11400: 724, -11300: 694, -11200: 649, -11100: 700, -11000: 733, -10900: 736, -10800: 715, -10700: 715, -10600: 705, -10500: 740, -10400: 691, -10300: 720, -10200: 744, -10100: 660}

Viendo M_hist[7], correspondiente a secuencia TGAAGTG
Largo: 8989
{-11502: 623, -11400: 599, -11300: 581, -11200: 609, -11100: 602, -11000: 614, -10900: 561, -10800: 648, -10700: 594, -10600: 632, -10500: 580, -10400: 582, -10300: 605, -10200: 562, -10100: 534}

Viendo M_hist[8], correspondiente a secuencia TTAAGTG
Largo: 7365
{-11502: 513, -11400: 488, -11300: 494, -11200: 473, -11100: 483, -11000: 507, -10900: 474, -10800: 514, -10700: 508, -10600: 472, -10500: 476, -10400: 467, -10300: 506, -10200: 488, -10100: 460}

Viendo M_hist[9], correspondiente a secuencia GTAAGTG
Largo: 5883
{-11502: 360, -11400: 400, -11300: 454, -11200: 389, -11100: 405, -11000: 324, -10900: 460, -10800: 378, -10700: 369, -10600: 383, -10500: 401, -10400: 357, -10300: 441, -10200: 414, -10100: 348}

Histograma total
Largo: 81986
{-11502: 5502, -11400: 5484, -11300: 5485, -11200: 5381, -11100: 5434, -11000: 5438, -10900: 5524, -10800: 5503, -10700: 5403, -10600: 5525, -10500: 5462, -10400: 5338, -10300: 5538, -10200: 5556, -10100: 5056}
'''
