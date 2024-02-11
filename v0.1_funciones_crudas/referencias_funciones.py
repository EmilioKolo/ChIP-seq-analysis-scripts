
import os
import math
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


##################################### GENOMAS #####################################

mm9 = EnsemblRelease(67, species='mouse');
hg19 = EnsemblRelease(102, species='human');
GRCm38 = EnsemblRelease(102, species='mouse');


#################################### FUNCIONES ####################################


def abrir_matriz_pbm(nombre_matriz, sep_pbm='\t'):
    # Funcion que abre una matriz PBM (nombre_matriz incluye extension)
    # La matriz PBM tiene que incluir A:/C:/G:/T: en la columna 0
    # Devuelve una matriz con 4 listas correspondientes a A/C/G/T

    # M_out para guardar el output
    M_out = [];
    nuc_a = [];
    nuc_c = [];
    nuc_g = [];
    nuc_t = [];
    # Abro el archivo nombre_matriz
    with open(nombre_matriz, 'r') as F_PBM:
        # Leo cada linea y busco por A:/C:/G:/T: en los primeros dos lugares
        for curr_line in F_PBM:
            if curr_line[0:2]=='A:':
                nuc_a = curr_line.rstrip().split(sep_pbm)[1:];
                for i in range(len(nuc_a)):
                    nuc_a[i] = float(nuc_a[i]);
            elif curr_line[0:2]=='C:':
                nuc_c = curr_line.rstrip().split(sep_pbm)[1:];
                for i in range(len(nuc_c)):
                    nuc_c[i] = float(nuc_c[i]);
            elif curr_line[0:2]=='T:':
                nuc_t = curr_line.rstrip().split(sep_pbm)[1:];
                for i in range(len(nuc_t)):
                    nuc_t[i] = float(nuc_t[i]);
            elif curr_line[0:2]=='G:':
                nuc_g = curr_line.rstrip().split(sep_pbm)[1:];
                for i in range(len(nuc_g)):
                    nuc_g[i] = float(nuc_g[i]);
            else:
                print('No se reconoce:')
                print(curr_line.rstrip())
    M_out = [nuc_a, nuc_c, nuc_g, nuc_t];
    return M_out


def buscar_chr(gen):
    # Recibe un gen en formato Ensembl y devuelve el contig correspondiente en formato chr_n (e.g. chr1, chrX)

    r = 'chr';
    gen_contig = str(gen.contig);
    L_contigs = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16',
                 '17','18','19','20','21','22','23','X','Y','M'];
    if gen_contig.upper() in L_contigs:
        r = r + gen_contig.upper();
    elif gen_contig.upper() == 'MT':
        r = 'chrM';
    else:
        print('Error con contig: ' + gen_contig)
        print(gen)
    
    return r


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


def buscar_en_seq_2dir_unificado(busq, seq):
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
            n = seq_pos+curr_pos;
            L_pos.append(n-len(seq));
            seq_pos = seq_pos + curr_pos + 1;
    # Despues vuelvo a correr con la secuencia reversa
    rev_seq = complemento_secuencia(seq, adn=True);
    loop_1 = True;
    curr_pos = 0;
    seq_pos = 0;
    while loop_1:
        curr_pos, loop_1 = buscar_en_secuencia(busq, rev_seq[seq_pos:]);
        if loop_1:
            n = seq_pos+curr_pos;
            L_pos.append(-n-len(busq));
            seq_pos = seq_pos + curr_pos + 1;
    return L_pos


def buscar_pos0(gen):
    # Recibe un gen en formato Ensembl y devuelve la posicion del +1 y un booleano para la direccion

    if gen.strand == '+':
        pos0 = gen.start;
        forward = True;
    elif gen.strand == '-':
        pos0 = gen.end;
        forward = False;
    else:
        print('ERROR PARSEANDO STRAND DE GEN:')
        print(gen)
        pos0 = gen.start;
        forward = True;

    return pos0, forward


def complemento(N,adn=True):
    # Devuelve el complemento de un nucleotido en adn o arn
    dict_adn = {'T':'A','U':'A','A':'T','C':'G','G':'C','N':'N'};
    dict_arn = {'T':'A','U':'A','A':'U','C':'G','G':'C','N':'N'};
    if adn:
        ret = dict_adn[N];
    else:
        ret = dict_arn[N];
    return(ret)


def complemento_secuencia(seq, adn=True):
    # Devuelve el complemento de una secuencia de adn o arn
    # La secuencia del complemento se devuelve en la misma orientacion (al reves que la referencia)
    ret_seq = '';
    for i in range(len(seq)):
        ret_seq = ret_seq + complemento(seq[-i-1],adn=adn);
    return(ret_seq)


def ConsultaSecuencia(id_chr, seq_start, seq_finish, strand=1, sleep_time=60):
    # Devuelve una secuencia dado un ID de cromosoma (incluye info de especie) y posicion inicial/final
    time.sleep(0.1);
    rec_seq = '';
    try:
        handle = Entrez.efetch(db='nucleotide', id=id_chr, rettype='fasta',
                               strand=strand, seq_start=seq_start, seq_stop=seq_finish);
        record = SeqIO.read(handle, 'fasta');
        handle.close();
        rec_seq = record.seq;
    except:
        print('Exception raised for chr ' + str(id_chr) + ' between positions ' +
              str(seq_start) + ' and ' + str(seq_finish) + '.')
        time.sleep(sleep_time);
        try:
            handle = Entrez.efetch(db='nucleotide', id=id_chr, rettype='fasta',
                                   strand=strand, seq_start=seq_start, seq_stop=seq_finish);
            record = SeqIO.read(handle, 'fasta');
            handle.close();
            rec_seq = record.seq;
        except:
            print('Retry failed. Returning empty string.')
    return rec_seq


def extraer_columnas_csv(dir_archivo, nom_archivo, L_col, ext='.csv', sep_arch=';'):
    # Abre un archivo .csv y devuelve las columnas en L_col como matriz
    # Si L_col es una lista vacia, se devuelven todas las columnas

    # Creo la matriz que se devuelve
    M_out = [];
    # Abro el archivo .csv
    with open(str(dir_archivo) + str(ext), 'r') as F:
        print('Archivo ' + str(nom_archivo) + str(ext) + ' abierto.')
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
            # Guardo L_out en M_out
            M_out.append(L_out[:]);
    return M_out


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
    maxfreq = n.max();
    cutoff_freq = 10 ** (math.ceil(math.log10(maxfreq))-1);
    plt.ylim(ymax=np.ceil(maxfreq / cutoff_freq) * cutoff_freq if maxfreq % cutoff_freq else maxfreq + cutoff_freq);
    plt.xlim(xmin=rango_hist[0], xmax=rango_hist[1]);
    if guardar:
        plt.savefig('Histograma_' + str(nombre_histograma) + '.png');
    else:
        plt.show();
    plt.close();
    return L_numbers, dict_hist


def IDchr(chromosome,genome='hg19'):
    # Devuelve el ID de cromosoma para ConsultaSecuencia dado el numero de cromosoma y el genoma correspondiente
    # Programado a mano, funciona solo con hg19 (humano) y mm9 (raton)
    ret = '';
    b = True;
    if genome.lower() == 'hg19': # GRCh38.p13
        dict_IDchr = {'1':'NC_000001.11', '2':'NC_000002.12', '3':'NC_000003.12', '4':'NC_000004.12',
                      '5':'NC_000005.10', '6':'NC_000006.12', '7':'NC_000007.14', '8':'NC_000008.11',
                      '9':'NC_000009.12', '10':'NC_000010.11', '11':'NC_000011.10', '12':'NC_000012.12',
                      '13':'NC_000013.11', '14':'NC_000014.9', '15':'NC_000015.10', '16':'NC_000016.10',
                      '17':'NC_000017.11', '18':'NC_000018.10', '19':'NC_000019.10', '20':'NC_000020.11',
                      '21':'NC_000021.9', '22':'NC_000022.11', 'X':'NC_000023.11', 'Y':'NC_000024.10',
                      'M':'NC_012920.1', 'MT':'NC_012920.1'};
    elif genome.lower() == 'mm9': # MGSCv37
        dict_IDchr = {'1':'NC_000067.5', '2':'NC_000068.6', '3':'NC_000069.5', '4':'NC_000070.5',
                      '5':'NC_000071.5', '6':'NC_000072.5', '7':'NC_000073.5', '8':'NC_000074.5',
                      '9':'NC_000075.5', '10':'NC_000076.5', '11':'NC_000077.5', '12':'NC_000078.5',
                      '13':'NC_000079.5', '14':'NC_000080.5', '15':'NC_000081.5', '16':'NC_000082.5',
                      '17':'NC_000083.5', '18':'NC_000084.5', '19':'NC_000085.5', 'X':'NC_000086.6',
                      'Y':'NC_000087.6', 'M':'NC_005089.1', 'MT':'NC_005089.1'};
    elif genome.lower() == 'mouse102' or genome.lower() == 'grcm38':
        dict_IDchr = {'1':'CM000994.3', '10':'CM001003.3', '11':'CM001004.3', '12':'CM001005.3',
                      '13':'CM001006.3', '14':'CM001007.3', '15':'CM001008.3', '16':'CM001009.3',
                      '17':'CM001010.3', '18':'CM001011.3', '19':'CM001012.3', '2':'CM000995.3',
                      '3':'CM000996.3', '4':'CM000997.3', '5':'CM000998.3', '6':'CM000999.3',
                      '7':'CM001000.3', '8':'CM001001.3', '9':'CM001002.3', 'MT':'AY172335.1',
                      'X':'CM001013.3', 'Y':'CM001014.3', 'M':'AY172335.1'};
    else:
        print('No se pudo encontrar genoma ' + str(genome))
        b = False;
        dict_IDchr = {};

    if str(chromosome).upper() in dict_IDchr.keys():
        ret = dict_IDchr[str(chromosome).upper()];
    elif b:
        print('No se pudo encontrar cromosoma ' + str(chromosome))
    return ret


def largo_archivo(dir_arch, nom_arch, ext):
    # Devuelve el largo (en filas) de un archivo en la direccion dir_arch
    # nom_arch y ext solo sirven para display
    with open(dir_arch, 'r') as F:
        print('> Archivo ' + str(nom_arch) + str(ext) + ' abierto.')
        l_arch = len(F.read().split('\n'));
    return l_arch


def secuencias_preferidas_pbm(nombre_matriz, pos_ini, pos_end, sep_pbm='\t'):
    # Funcion que abre una matriz PBM y devuelve una lista con las secuencias con mas puntaje
    # La matriz PBM empieza en fila_ini y col_ini es 1 si se incluye A/C/G/T en la columna 0

    # Inicializo la lista de secuencias que devuelve la funcion
    L_seq = [];
    # Cargo la matriz de pesos en M_PBM
    M_PBM = abrir_matriz_pbm(nombre_matriz, sep_pbm=sep_pbm);

    L_seq = M_PBM;
    return L_seq


def _consulta_secuencia_simple(id_chr, seq_start, seq_finish, strand=1):
    # Prueba de consulta de secuencia sin try/except
    handle = Entrez.efetch(db='nucleotide', id=id_chr, rettype='fasta',
                           strand=strand, seq_start=seq_start, seq_stop=seq_finish);
    record = SeqIO.read(handle, 'fasta');
    handle.close();
    rec_seq = record.seq;
    return rec_seq


def _test_secuencia_correcta(seq_buscada, chr_id, pos_ini, pos_end, verbose=False):
    # Revisa si una combinacion chr, pos_ini, pos_end se corresponde a una secuencia

    # Inicializo el booleano para devolver
    b = False;

    # Agarro la secuencia dada por chr_id, pos_ini, pos_end
    seq_en_pos = ConsultaSecuencia(chr_id, pos_ini, pos_end);

    # Si la secuencia encontrada es igual a la buscada, se devuelve True
    if seq_en_pos == seq_buscada or seq_en_pos == complemento_secuencia(seq_buscada):
        if verbose:
            print('OK!')
        b = True;
    # Sino se deja el valor default (False)
    elif verbose:
        print('Secuencia diferente. A continuacion agrego 10 bases antes y despues:')
        print(ConsultaSecuencia(chr_id, pos_ini-10, pos_end+10))
    return b


##################################### TESTEO ######################################

if __name__=='__main__':
    n_pbm = 'Matriz PBM.txt';
    M_out = secuencias_preferidas_pbm(n_pbm, 0, 1, sep_pbm='\t');
    print(M_out)

