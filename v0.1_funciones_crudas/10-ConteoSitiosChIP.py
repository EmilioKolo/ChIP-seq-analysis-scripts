
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
from referencias_funciones import ConsultaSecuencia, IDchr 

mm9 = EnsemblRelease(67, species='mouse');
hg19 = EnsemblRelease(102, species='human');
GRCm38 = EnsemblRelease(102, species='mouse');

'''
Cuenta cuantos de los picos de ChIP-seq de los archivos seleccionados tienen sitios de union

Devuelve los resultados en display (y tarda un rato)
'''

#################################### VARIABLES ####################################

# Variables de archivos fuente

lista_nombres_peaks = ['Anderson2018', 'Dupays2015', 'Dupays2015_s1e14',
                       'Dupays2015_s4e14', 'He2011', 'vandenBoogard2012'];
lista_nombres_peaks_hg19 = ['Anderson2018'];
lista_nombres_peaks_mm9 = ['Dupays2015', 'Dupays2015_s1e14', 'Dupays2015_s4e14',
                           'He2011', 'vandenBoogard2012'];

path_archivos_peaks = '.\\csv con secuencias\\';

# Variables generales

sitios_union_nkx25 = ['GCAAGTG', 'GGAAGTG', 'GAAAGTG', 'ATAAGTG', 'GTAAGTG',
                      'CTAAGTG', 'TCAAGTG', 'TGAAGTG', 'TAAAGTG', 'TTAAGTG']; 

ext_csv = '.csv';
sep_csv = ';';
ext_bed = '.bed';
sep_bed = ';';
curr_path = os.path.dirname(__file__);

# Variables main()

L_nom_arch_main = lista_nombres_peaks;
L_sitios_main = sitios_union_nkx25;


#################################### FUNCIONES ####################################


def abrir_resultados_csv(dir_arch, nom_arch, ext_arch='.csv', sep_arch=';'):
    # Abre el archivo en dir_arch y devuelve todas las columnas como matriz
    # Default extension .csv pero se puede cambiar

    # Inicializo la matriz a devolver
    M_csv = [];

    # Abro el archivo
    with open(dir_arch + ext_arch, 'r') as F:
        print('Archivo ' + str(nom_arch) + str(ext_arch) + ' abierto.')
        for curr_line in F:
            L_curr_line = curr_line.rstrip().split(sep_arch);
            M_csv.append(L_curr_line[:])
    return M_csv


def buscar_en_secuencia(busq, seq):
    # Busca una secuencia "busq" en una secuencia mas larga "seq"
    # Devuelve la posicion en la que se encontro busq y True si se encuentra busq
    # Devuelve len(seq)-len(busq)+1 y False si no se encuentra busq
    pos = 0;
    encontrado = False;
    # Recorro una por una las posiciones de seq hasta encontrar busq
    while pos < (len(seq)-len(busq)+1) and (not encontrado):
        if seq[pos:pos+len(busq)] == busq:
            encontrado = True;
        else:
            pos += 1;
    return pos, encontrado


def buscar_en_seq_ambas_direcciones(busq, seq):
    # Funcion que usa buscar_en_secuencia() para buscar todas las ocurrencias de busq en seq
    # Busca en ambas direcciones y registra varias ocurrencias
    L_pos = [];
    loop_1 = True;
    curr_pos = 0;
    seq_pos = 0;
    # Primero corro buscar_en_secuencia() para busq y seq normalmente
    while loop_1:
        curr_pos, loop_1 = buscar_en_secuencia(busq, seq[seq_pos:]);
        if loop_1:
            L_pos.append(seq_pos+curr_pos);
            seq_pos = seq_pos + curr_pos + 1;
    # Despues vuelvo a correr con la secuencia reversa
    rev_seq = complemento_secuencia(seq, adn=True);
    loop_1 = True;
    curr_pos = 0;
    seq_pos = 0;
    while loop_1:
        curr_pos, loop_1 = buscar_en_secuencia(busq, rev_seq[seq_pos:]);
        if loop_1:
            L_pos.append(-(seq_pos+curr_pos));
            seq_pos = seq_pos + curr_pos + 1;
    return L_pos


def complemento(N,adn=True):
    # Devuelve el complemento de un nucleotido en adn o arn
    dict_adn = {'T':'A','U':'A','A':'T','C':'G','G':'C','N':'N'};
    dict_arn = {'T':'A','U':'A','A':'U','C':'G','G':'C','N':'N'};
    if adn:
        ret = dict_adn[N];
    else:
        ret = dict_arn[N];
    return ret


def complemento_secuencia(seq, adn=True):
    # Devuelve el complemento de una secuencia de adn o arn
    # La secuencia del complemento se devuelve en la misma orientacion (al reves que la referencia)
    ret_seq = '';
    for i in range(len(seq)):
        ret_seq = ret_seq + complemento(seq[-i-1],adn=adn);
    return ret_seq


def pipeline_main_por_arch(dir_arch, nom_arch, sep_arch, L_sitios, ext_arch='.csv', n_verb=1000):
    #

    # Inicializo la matriz con los datos de peaks
    M_SU = abrir_resultados_csv(dir_arch, nom_arch, ext_arch=ext_arch, sep_arch=sep_arch);

    # Inicializo contador de peaks con sitios de union
    cont_SU = 0;
    # Contador en formato diccionario para ver distribucion por sitio
    dict_cont = {};
    for s in L_sitios:
        dict_cont[s] = 0;
    # Display
    largo_M_SU = len(M_SU);
    # Reviso cada peak
    for i in range(largo_M_SU):
        curr_peak = M_SU[i];
        curr_seq = curr_peak[4];
        # Para determinar si el peak tiene por lo menos un sitio de union
        peak_con_su = False;
        # Reviso con cada sitio de union
        for s in range(len(L_sitios)):
            curr_su = L_sitios[s];
            # Busco curr_su en curr_seq
            L_pos_su = buscar_en_seq_ambas_direcciones(curr_su, curr_seq);
            # Si curr_su aparece por lo menos una vez, hago peak_con_su True y anoto en dict_cont
            if len(L_pos_su) > 0:
                peak_con_su = True;
                dict_cont[curr_su] = dict_cont[curr_su] + 1;
        # Habiendo revisado todos los sitios de union, veo si peak_con_su es True
        if peak_con_su:
            cont_SU = cont_SU + 1;
        # Display
        if (i+1)%n_verb==0:
            print('Analizado peak numero ' + str(i+1) + ' de ' + str(largo_M_SU))
    # Al final, hago print de dict_cont y cont_SU
    print('Peaks con sitios de union encontrados: ' + str(cont_SU) + ' de ' + str(largo_M_SU))
    print('Distribucion de los tipos de sitios encontrados:')
    print(dict_cont)
    # Meto cont_SU en dict_cont para devolver un solo item
    dict_cont['total'] = cont_SU;
    return dict_cont


def _main():
    #

    # Inicializo la matriz que se devuelve
    M_out = [];

    # Corro pipeline_main_por_arch() para cada nombre de archivo en L_nom_arch_main
    for nom_arch_main in L_nom_arch_main:
        print('Pipeline con archivo ' + nom_arch_main + '.csv')
        dir_arch_main = os.path.join(curr_path, path_archivos_peaks + nom_arch_main);
        dict_cont = pipeline_main_por_arch(dir_arch_main, nom_arch_main, sep_csv, L_sitios_main);
        # Copio dict_cont para evitar errores
        M_out.append(dict_cont.copy());
        print()
        print()
    return M_out


def _main_test():
    #

    #
    M_out = [];

    return M_out


#################################### RESULTADOS ###################################

'''
Pipeline con archivo Anderson2018.csv
Peaks con sitios de union encontrados: 3846 de 6871
Distribucion de los tipos de sitios encontrados:
{'GCAAGTG': 362, 'GGAAGTG': 632, 'GAAAGTG': 648, 'ATAAGTG': 437, 'GTAAGTG': 389, 'CTAAGTG': 409, 'TCAAGTG': 872, 'TGAAGTG': 638, 'TAAAGTG': 607, 'TTAAGTG': 566}

Pipeline con archivo Dupays2015.csv
Peaks con sitios de union encontrados: 2022 de 2610
Distribucion de los tipos de sitios encontrados:
{'GCAAGTG': 321, 'GGAAGTG': 441, 'GAAAGTG': 347, 'ATAAGTG': 214, 'GTAAGTG': 228, 'CTAAGTG': 359, 'TCAAGTG': 639, 'TGAAGTG': 537, 'TAAAGTG': 287, 'TTAAGTG': 385}

Pipeline con archivo Dupays2015_s1e14.csv
Peaks con sitios de union encontrados: 3375 de 4418
Distribucion de los tipos de sitios encontrados:
{'GCAAGTG': 537, 'GGAAGTG': 775, 'GAAAGTG': 570, 'ATAAGTG': 347, 'GTAAGTG': 397, 'CTAAGTG': 584, 'TCAAGTG': 994, 'TGAAGTG': 851, 'TAAAGTG': 446, 'TTAAGTG': 628}

Pipeline con archivo Dupays2015_s4e14.csv
Peaks con sitios de union encontrados: 4378 de 5441
Distribucion de los tipos de sitios encontrados:
{'GCAAGTG': 730, 'GGAAGTG': 1144, 'GAAAGTG': 875, 'ATAAGTG': 482, 'GTAAGTG': 573, 'CTAAGTG': 768, 'TCAAGTG': 1245, 'TGAAGTG': 1158, 'TAAAGTG': 689, 'TTAAGTG': 873}

Pipeline con archivo He2011.csv
Peaks con sitios de union encontrados: 11176 de 20573
Distribucion de los tipos de sitios encontrados:
{'GCAAGTG': 1445, 'GGAAGTG': 1693, 'GAAAGTG': 1351, 'ATAAGTG': 914, 'GTAAGTG': 1256, 'CTAAGTG': 1617, 'TCAAGTG': 2507, 'TGAAGTG': 2244, 'TAAAGTG': 1249, 'TTAAGTG': 1637}

Pipeline con archivo vandenBoogard2012.csv
Peaks con sitios de union encontrados: 2242 de 6705
Distribucion de los tipos de sitios encontrados:
{'GCAAGTG': 266, 'GGAAGTG': 368, 'GAAAGTG': 188, 'ATAAGTG': 107, 'GTAAGTG': 196, 'CTAAGTG': 276, 'TCAAGTG': 513, 'TGAAGTG': 416, 'TAAAGTG': 129, 'TTAAGTG': 204}

'''

output_dump = '';

if __name__=='__main__':
    output_dump = _main();
    #output_dump = _main_test();
