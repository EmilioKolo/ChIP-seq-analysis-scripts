
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
Genera una tabla de sitios de union en promotores del genoma y peaks de ChIP-seq

Pipeline que corre todos los programas juntos
    * pipeline_main()

Abre archivos .csv con sitios de union en promotores y les agrega info apropiada
    * devolver_SUG()
    Se puede testear con _prueba_secuencia()

Abre archivos .csv con peaks de ChIP-seq, busca sitios de union y los registra
    * devolver_SUP()

Abre archivos .csv sacados de INSECT2 con sitios de union (formato especifico en mis resultados)
    * devolver_INSECT()
    * No todos los sitios testean bien, no usar sitios si no estan confirmados
'''

#################################### VARIABLES ####################################

# Variables sitios de union genoma en .csv (SUG = Sitios de Union del Genoma)

sitios_union_genoma_mm9 = 'SitiosDeUnionGenoma';
sitios_union_genoma_hg19 = 'SitiosDeUnionGenoma_hg19';

path_sitios_union_genoma = '.\\sitios de union genoma\\';

# Variables peaks ChIP-seq con secuencia en .csv (SUP = Sitios de Union en Peaks)

lista_nombres_peaks = ['Anderson2018', 'Dupays2015', 'Dupays2015_s1e14',
                       'Dupays2015_s4e14', 'He2011', 'vandenBoogard2012'];
lista_nombres_peaks_hg19 = ['Anderson2018'];
lista_nombres_peaks_mm9 = ['Dupays2015', 'Dupays2015_s1e14', 'Dupays2015_s4e14',
                           'He2011', 'vandenBoogard2012'];
L_genomes = [hg19, mm9, mm9, mm9, mm9, mm9];

path_archivos_peak = '.\\csv con secuencias\\';

# Variables resultados INSECT2

nombre_resultados_INSECT = 'NKX2-5';
path_resultados_INSECT = '..\\0-Fuentes\\INSECT2\\';

# Variables generales

sitios_union_nkx25 = ['GCAAGTG', 'GGAAGTG', 'GAAAGTG', 'ATAAGTG', 'GTAAGTG',
                      'CTAAGTG', 'TCAAGTG', 'TGAAGTG', 'TAAAGTG', 'TTAAGTG']; 

ext_csv = '.csv';
sep_csv = ';';
curr_path = os.path.dirname(__file__);

# Variables main()

genoma_main = mm9;
nombre_genoma_main = 'mm9';
L_genoma_main = [hg19, mm9, mm9, mm9, mm9, mm9];
L_nombre_genoma_main = ['hg19', 'mm9', 'mm9', 'mm9', 'mm9', 'mm9'];
genoma_INSECT = hg19;
nombre_genoma_INSECT = 'hg19';
SUG_main = sitios_union_genoma_mm9;
dir_SUG_main = os.path.join(curr_path, path_sitios_union_genoma + SUG_main);
dir_INSECT_main = os.path.join(curr_path, path_resultados_INSECT + nombre_resultados_INSECT);
verbose_SUG1 = False;
SUG_test = False;
SUG_test_n = 20;
verbose_SUP1 = False;
SUP_test = False;
SUP_test_n = 20;
INSECT_test = True;
SU_repetido_test = True;
n_SU_repetidos = 1000;

correr_INSECT = True;
correr_SUG = True;
correr_SUP = True;

verbose_main = True;
n_main_display = 1000;

SUG_main_mm9 = sitios_union_genoma_mm9;
SUG_main_hg19 = sitios_union_genoma_hg19;
dir_SUG_main_mm9 = os.path.join(curr_path, path_sitios_union_genoma + SUG_main_mm9);
dir_SUG_main_hg19 = os.path.join(curr_path, path_sitios_union_genoma + SUG_main_hg19);
nombre_out_main = 'SitiosUnionFinal';
nombre_out_mm9 = nombre_out_main + '_mm9';
nombre_out_hg19 = nombre_out_main + '_hg19';


#################################### FUNCIONES ####################################


def abrir_csv_peaks(dir_arch, nom_arch, genoma_usado, nombre_genoma, sep_arch=';'):
    # Abre el archivo .csv en dir_arch y devuelve todas las columnas como matriz
    # Revisa que haya secuencias para todos los peaks
    M_csv = abrir_resultados_csv(dir_arch, nom_arch, sep_arch=sep_arch);

    # Reviso que todos los peaks tengan secuencias
    for i in range(len(M_csv)):
        curr_peak = M_csv[i];
        if len(curr_peak) < 5:
            print('Peak con largo menor a 5:')
            print(curr_peak)
            # Saco chromosome ID con la funcion IDchr()
            #chr_id = IDchr(curr_peak[1], genome=nombre_genoma);
            # Busco la secuencia correspondiente
            #seq_encontrada = str(ConsultaSecuencia(chr_id, curr_peak[2], curr_peak[3], strand=1, sleep_time=5));
            # Agrego la secuencia al final de la lista
            #M_csv[i].append(seq_encontrada);
            #print('Secuencia encontrada: ')
            #print(seq_encontrada)
        elif len(curr_peak[4]) < (int(curr_peak[3])-int(curr_peak[2])):
            print('Peak con secuencia menor al esperado')
            print(curr_peak)
            # Saco chromosome ID con la funcion IDchr()
            #chr_id = IDchr(curr_peak[1], genome=nombre_genoma);
            # Busco la secuencia correspondiente
            #seq_encontrada = str(ConsultaSecuencia(chr_id, curr_peak[2], curr_peak[3], strand=1, sleep_time=5));
            # Agrego la secuencia al final de la lista
            #M_csv[i][4] = seq_encontrada;
            #print('Secuencia encontrada: ')
            #print(seq_encontrada)
    return M_csv


def abrir_resultados_csv(dir_arch, nom_arch, sep_arch=';'):
    # Abre el archivo .csv en dir_arch y devuelve todas las columnas como matriz

    # Inicializo la matriz a devolver
    M_csv = [];

    # Abro el archivo
    with open(dir_arch + '.csv', 'r') as F:
        print('Archivo ' + str(nom_arch) + '.csv abierto.')
        for curr_line in F:
            L_curr_line = curr_line.rstrip().split(sep_arch);
            M_csv.append(L_curr_line[:]);
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
    return ret


def complemento_secuencia(seq, adn=True):
    # Devuelve el complemento de una secuencia de adn o arn
    # La secuencia del complemento se devuelve en la misma orientacion (al reves que la referencia)
    ret_seq = '';
    for i in range(len(seq)):
        ret_seq = ret_seq + complemento(seq[-i-1],adn=adn);
    return ret_seq 


def devolver_INSECT(dir_INSECT, nom_INSECT, genoma_usado, nombre_genoma, sep_arch=';', hacer_pruebas=False, rango_prueba=50000):
    # Abre archivo .csv con sitios de union encontrados usando INSECT2
    # Los parsea y arregla en formato secuencia_sitio, gene_id, gene_name, chr_n, pos_ini, pos_end, INSECT

    # Inicializo la matriz que se devuelve
    M_out = [];
    
    # Extraigo los datos del csv
    M_csv = abrir_resultados_csv(dir_INSECT, nom_INSECT, sep_arch=sep_arch);

    if hacer_pruebas:
        shuffle(M_csv);

    # Contador para genes no encontrados
    cont_no_encontrado = 0;
    # Recorro cada sitio de union
    for i in range(len(M_csv)):
        curr_sitio = M_csv[i];
        gen_encontrado = False;
        try:
            curr_gen = genoma_usado.gene_by_id(curr_sitio[0]);
            gen_encontrado = True;
        except:
            curr_gen = None;
            gen_encontrado = False;
            cont_no_encontrado = cont_no_encontrado + 1;

            # Display
            if cont_no_encontrado%rango_prueba == 0:
                print('IDs de genes no encontrados: ' + str(cont_no_encontrado))

        if gen_encontrado:
            curr_chr, curr_pos_ini, curr_pos_end = extraer_sitio_gen(curr_gen, curr_sitio[1], curr_sitio[2]);
            # Inicializo la lista que se agrega a M_out para devolver
            sitio_append = [];
            # Agrego toda la info a sitio_append
            sitio_append.append(str(curr_sitio[7]));
            sitio_append.append(str(curr_sitio[0]));
            sitio_append.append(curr_gen.gene_name);
            sitio_append.append(str(curr_chr));
            sitio_append.append(int(curr_pos_ini));
            sitio_append.append(int(curr_pos_end));
            sitio_append.append('INSECT');
            # Agrego el sitio definido en sitio_append a M_out
            M_out.append(sitio_append[:]);
            if hacer_pruebas and (i%rango_prueba == 0):
                print(curr_sitio)
                print(curr_gen)
                print(sitio_append)
                _prueba_secuencia(curr_chr, curr_pos_ini, curr_pos_end, curr_sitio[7], genoma_usado, nombre_genoma, verbose=True);
        # Display
        if (i%rango_prueba == 0):
            print('Sitios medidos: ' + str(i) + ' / ' + str(len(M_csv)))
    print('IDs de genes no encontrados (fin de ejecucion): ' + str(cont_no_encontrado))
    return M_out


def devolver_SUG(dir_SUG, nom_SUG, sep_arch, genoma_usado, verbose=False, verbose_range=5000):
    # Abre archivo .csv de sitios de union en promotores del genoma y los devuelve formateados
    # Los sitios de union se devuelven en formato secuencia_sitio, gen_id, gene_name, chr_n, pos_ini, pos_end

    # Inicializo la matriz que se devuelve
    M_out = [];

    # Agarro la info del .csv
    M_out = abrir_resultados_csv(dir_SUG, nom_SUG, sep_arch=sep_arch);

    if verbose:
        print('Iniciando procesamiento de M_out.')
    
    # Recorro M_out y agrego info a cada sitio de union
    for i in range(len(M_out)):
        curr_sitio = M_out[i];
        curr_sitio_parseado = [str(curr_sitio[1]), str(curr_sitio[0])] + parsear_sitio_csv(curr_sitio, genoma_usado);

        # Por si es necesario separar el cromosoma, pos_ini y pos_end del sitio parseado
        #curr_chr = curr_sitio_parseado[0];
        #curr_pos_ini = curr_sitio_parseado[1];
        #curr_pos_end = curr_sitio_parseado[2];

        M_out[i] = curr_sitio_parseado[:];
        
        if verbose and i%verbose_range == 0:
            print(curr_sitio_parseado)
    
    return M_out


def devolver_SUP(L_nom_peaks, path_peaks, L_sitios_union, L_genoma_usado, L_nombre_genoma, verbose=True, hacer_pruebas=False):
    # Abre varios archivos .csv con peaks de ChIP-seq y devuelve una lista de todos los sitios de union en esos peaks
    # Sitios de union en L_sitios_union
    # Nombres de archivos en L_nom_peaks, ubicados en path_peaks
    # Genomas para chequeos

    # Inicializo la matriz que se devuelve
    M_out = [];

    # Inicializo la lista de referencias al final del sitio de union
    L_ref = [];
    for i in range(len(L_nom_peaks)):
        L_ref.append(False);
    #print('PRUEBA: L_ref inicial: ' + str(L_ref))
    
    # Reviso cada nombre de archivo en L_nom_peaks
    for i in range(len(L_nom_peaks)):
        curr_nom_peak = L_nom_peaks[i];
        if verbose:
            print('Abriendo archivo ' + str(curr_nom_peak) + '.csv')
        curr_path_peak = os.path.join(curr_path, path_peaks + curr_nom_peak);
        # Abro el archivo curr_nom_peak y guardo los sitios de union en una matriz
        M_peaks = abrir_csv_peaks(curr_path_peak, curr_nom_peak, L_genoma_usado[i], L_nombre_genoma[i]);
        if verbose:
            print()
            print('Iniciando busqueda de sitios de union en archivo ' + str(curr_nom_peak) + '.csv')
        #print('PRUEBA: L_ref sin cambiar: ' + str(L_ref))
        # Actualizo L_ref
        if True in L_ref:
            L_ref[i-1] = False;
        L_ref[i] = True;
        #print('PRUEBA: L_ref cambiado: ' + str(L_ref))
        # Agrego sitios de union encontrados en M_peaks a M_out
        M_out = M_out + encontrar_SUP(M_peaks, L_sitios_union, L_genoma_usado[i], L_nombre_genoma[i], L_ref, hacer_pruebas=hacer_pruebas);

    return M_out


def encontrar_SUP(M_peaks_seq, L_SU, genoma_usado, nombre_genoma, L_ref, hacer_pruebas=False, rango_prueba=2000):
    # Busca todos los sitios de union en una matriz de peaks
    # Sitios de union definidos en L_SU
    # Peaks en formato [su_seq, gene_id, gene_name, contig, peak_ini, peak_end, peak_seq]

    # Inicializo la matriz que se devuelve
    M_out = [];

    # Inicializo la lista de referencias al final del sitio de union
    L_ref_append = []
    for j in range(len(L_ref)):
        if L_ref[j]:
            L_ref_append.append('X');
        else:
            L_ref_append.append('');
    # Recorro la matriz de peaks
    for i in range(len(M_peaks_seq)):
        curr_peak = M_peaks_seq[i];
        curr_chr = curr_peak[1];
        curr_peak_ini = int(curr_peak[2]);
        curr_peak_end = int(curr_peak[3]);
        curr_peak_seq = curr_peak[4];
        # Busco cada secuencia de L_SU en la secuencia del peak
        for curr_seq in L_SU:
            # Cargo la lista de sitios de union para curr_seq
            L_sitios = buscar_en_seq_ambas_direcciones(curr_seq, curr_peak_seq);
            # Recorro cada sitio de union y lo registro en M_out
            for sitio_union in L_sitios:
                # Inicializo el sitio de union que se agrega a M_out
                curr_sitio_union = [];
                curr_sitio_union.append(str(curr_seq));
                # Agrego dos espacios vacios a la lista por gene_id y gene_name
                curr_sitio_union.append('');
                curr_sitio_union.append('');
                # Si sitio_union es mayor o igual a 0, busco la posicion sumando a peak_ini
                if sitio_union >= 0:
                    curr_pos_ini = curr_peak_ini + sitio_union;
                    curr_pos_end = curr_pos_ini + (len(curr_seq)-1);
                # Si es menor a 0, busco la posicion restando a peak_end
                else:
                    curr_pos_ini = curr_peak_end + sitio_union - (len(curr_seq)-1);
                    curr_pos_end = curr_pos_ini + (len(curr_seq)-1);
                # Agrego cromosoma, posicion inicial y final del sitio de union
                curr_chr_append = 'chr'+str(curr_chr);
                curr_sitio_union.append(str(curr_chr_append));
                curr_sitio_union.append(int(curr_pos_ini));
                curr_sitio_union.append(int(curr_pos_end));
                # Agrego el paper correspondiente al final del sitio de union
                curr_sitio_union = curr_sitio_union + L_ref_append[:];
                # Agrego el sitio de union procesado a M_out
                M_out.append(curr_sitio_union[:]);
                # Cada rango_prueba sitios pruebo que la posicion sea la correcta
                if hacer_pruebas and (len(M_out)%rango_prueba == 0):
                    print('Sitio de union a testear: ' + str([curr_chr,curr_pos_ini,curr_pos_end,curr_seq]))
                    _prueba_secuencia('chr'+curr_chr, curr_pos_ini, curr_pos_end, curr_seq, genoma_usado, nombre_genoma, verbose=True);
    return M_out


def extraer_sitio_gen(ensembl_gen, ref_ini, ref_end):
    # Usa un elemento ensembl_gen y posiciones de referencia para obtener info
    # Devuelve cromosoma, posicion inicial y final del sitio marcado por la referencia

    # Extraigo chr_n del contig
    chr_n = 'chr' + str(ensembl_gen.contig);

    # Defino pos0 del gen
    pos0, forward = buscar_pos0(ensembl_gen);

    # Calculo pos_ini y pos_end en base a pos0, ref_ini y ref_end
    if forward:
        pos_ini = pos0 + int(ref_ini);
        pos_end = pos0 + int(ref_end);
    else:
        pos_ini = pos0 - int(ref_ini) + 1;
        pos_end = pos0 - int(ref_end) + 1;

    return chr_n, pos_ini, pos_end


def guardar_csv(nombre_out, M_out, ext_out='.csv', sep_out=';'):
    # Guarda una matriz M_out en un archivo nombre_out en formato ext_out
    # Pensado para .csv separado por punto y coma

    # Creo el archivo vacio
    with open(nombre_out + ext_out, 'w') as F_out:
        print('Archivo ' + nombre_out + ext_out + ' creado.')
    # Lo vuelvo a abrir en modo append para guardar los contenidos de M_out
    with open(nombre_out + ext_out, 'a') as F_out:
        # Abro cada sitio de union
        for i in range(len(M_out)):
            curr_line = '';
            curr_SU = M_out[i];
            # Recorro cada componente del sitio de union y lo agrego a curr_line
            for j in range(len(curr_SU)):
                curr_line = curr_line + str(curr_SU[j]) + sep_out;
            # Agrego el fin de linea
            curr_line = curr_line + '\n';
            # Guardo curr_line en F_out
            F_out.write(curr_line);
    return M_out


def parsear_sitio_csv(L_sitio, genome):
    # Recibe un sitio en formato [gen_id, seq, relative_pos] y lo devuelve en formato [gene_name, chr_n, pos_ini, pos_end]

    # Inicializo la variable a devolver
    L_out = ['gene_name', 'chr', -1, -1];

    try:
        # Defino el gen correspondiente al sitio
        curr_gen = genome.gene_by_id(L_sitio[0]);
    except:
        print('Error buscando gen ' + str(L_sitio[0]) + ' en genoma ' + str(genome))
        return L_out

    # Agrego gene_name
    L_out[0] = curr_gen.gene_name;
    # Defino el cromosoma correspondiente
    curr_gen_contig = curr_gen.contig;
    L_contigs = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16',
                 '17','18','19','20','21','22','23','X','Y','M'];
    if str(curr_gen_contig).upper() in L_contigs:
        L_out[1] = 'chr' + str(curr_gen_contig).upper();
    elif str(curr_gen_contig).upper() == 'MT':
        L_out[1] = 'chrM';
    else:
        print('Error con contig: ' + str(curr_gen_contig))
        print(curr_gen)

    # Defino pos0 del gen
    pos0, forward = buscar_pos0(curr_gen);

    # Defino el largo del sitio de union con la secuencia en L_sitio
    largo_sitio = len(L_sitio[1]);
    # Defino pos_ini y pos_end
    if forward:
        pos_ini = pos0+int(L_sitio[2])+1;
        pos_end = pos0+int(L_sitio[2])+largo_sitio;
    else:
        pos_ini = pos0-int(L_sitio[2])-largo_sitio;
        pos_end = pos0-int(L_sitio[2])-1;
    L_out[2] = pos_ini;
    L_out[3] = pos_end;
    
    return L_out


def pipeline_main(dir_SUG, nom_SUG, L_nom_peaks, path_peaks, L_sitios_union, sep_csv, genoma,
                  nombre_genoma, insect2, dir_insect2, nom_insect2, verbose=False, n_verbose=10000):
    #

    # Inicializo la variable que se devuelve
    M_out = [];

    # Abrir archivos .csv en carpeta "sitios de union genoma"
    # SUG = Sitios de Union del Genoma
    M_SUG = devolver_SUG(dir_SUG, nom_SUG, sep_csv, genoma);
    if verbose:
        print('devolver_SUG() terminado. Ejemplos M_SUG: ')
        shuffle(M_SUG);
        for i in M_SUG[:20]:
            print(i)
        print()
    L_genomas = [];
    L_nombres_genomas = [];
    for i in L_nom_peaks:
        L_genomas.append(genoma);
        L_nombres_genomas.append(str(nombre_genoma));
    if verbose:
        print('Iniciando busqueda de sitios de union en peaks.')
        print()

    # Buscar sitios de union en peaks del .csv en carpeta "csv con secuencias"
    # SUP = Sitios de Union en Peaks
    M_SUP = devolver_SUP(L_nom_peaks, path_peaks, L_sitios_union, L_genomas, L_nombres_genomas);
    if verbose:
        print()
        print('devolver_SUP() terminado. Ejemplos M_SUP: ')
        shuffle(M_SUP);
        for i in M_SUP[:20]:
            print(i)
        print()
    if insect2:
        if verbose:
            print('Iniciando busqueda de sitios de union en INSECT2.')
            print()
        # Buscar sitios de union INSECT2
        M_INSECT = devolver_INSECT(dir_insect2, nom_insect2, genoma, nombre_genoma);
        if verbose:
            print()
            print('devolver_INSECT() terminado. Ejemplos M_INSECT: ')
            shuffle(M_INSECT)
            for i in M_INSECT[:20]:
                print(i)
            print()
    else:
        M_INSECT = [];

    # Creo la tabla para guardar al final
    if verbose:
        print('Iniciando generacion de tabla con toda la info.')
        print()
        print('Armando base con M_SUG.')
    # Recorro M_SUG
    for k in range(len(M_SUG)):
        curr_SU = M_SUG[k];
        # Agrego lugar para la info de experimentos ChIP-seq
        for p in range(len(L_nom_peaks)):
            curr_SU = curr_SU + [''];
        # Agrego lugar para la info de INSECT
        if insect2:
            curr_SU = curr_SU + [''];
        # Agrego curr_SU a M_out
        M_out.append(curr_SU[:]);

    # Una vez agregados los sitios de M_SUG, agrego M_SUP y M_INSECT
    if verbose:
        print('Agregando info de M_SUP.')
    # Contador para sitios repetidos
    cont_repetidos = 0;
    # Recorro M_SUP
    for k in range(len(M_SUP)):
        curr_SU = M_SUP[k];
        # Contadores para revisar presencia de curr_SU en M_out
        presente_en_M_out = False;
        j = 0;
        # Reviso si curr_SU ya esta en M_out
        while j<len(M_out) and (not presente_en_M_out):
            # Si curr_SU esta en M_out, agrego la info entre M_out[j][5] y M_out[j][-2]
            bool1 = str(curr_SU[3]) == str(M_out[j][3]);
            bool2 = str(curr_SU[4]) == str(M_out[j][4]);
            bool3 = str(curr_SU[5]) == str(M_out[j][5]);
            if bool1 and bool2 and bool3:
                presente_en_M_out = True;
                cont_repetidos = cont_repetidos + 1;
                if verbose and cont_repetidos%n_verbose == 0:
                    print('curr_SU repetidos encontrados en M_SUP: ' + str(cont_repetidos))
                # Recorro lista_nombres_peaks para definir en donde se encuentra curr_SU
                for h in range(len(L_nom_peaks)):
                    # Registro los valores no vacios
                    if curr_SU[6+h] != '':
                        #####################
                        ### PRUEBA INSECT ###
                        #if curr_SU[6+h] != 'X':
                        #    print('ERROR DE INDICE, 6+h=' + str(6+h))
                        #    print(curr_SU)
                        ### PRUEBA INSECT ###
                        #####################
                        M_out[j][6+h] = str(curr_SU[6+h]);
            j = j + 1;
        # Si no esta presente curr_SU en M_out, lo agrego
        if not presente_en_M_out:
            M_out.append(curr_SU[:]);
        # Display
        if verbose and k%n_verbose == 0:
                print('Revisando ' + str(k) + ' / ' + str(len(M_SUP)))

    # Recorro M_INSECT
    if verbose:
        print('curr_SU repetidos encontrados en M_SUP: ' + str(cont_repetidos))
        print()
        if insect2:
            print('Agregando info de M_INSECT.')
    # Contador para sitios repetidos
    cont_repetidos = 0;
    # Recorro M_INSECT
    for k in range(len(M_INSECT)):
        curr_SU = M_INSECT[k];
        # Contadores para revisar presencia de curr_SU en M_out
        presente_en_M_out = False;
        j = 0;
        # Reviso si curr_SU ya esta en M_out
        while j<len(M_out) and (not presente_en_M_out):
            # Si curr_SU esta en M_out, agrego la info en M_out[j][-1]
            bool1 = str(curr_SU[3]) == str(M_out[j][3]);
            bool2 = str(curr_SU[4]) == str(M_out[j][4]);
            bool3 = str(curr_SU[5]) == str(M_out[j][5]);
            if bool1 and bool2 and bool3:
                presente_en_M_out = True;
                cont_repetidos = cont_repetidos + 1;
                if verbose and cont_repetidos%n_verbose == 0:
                    print('curr_SU repetidos encontrados en M_INSECT: ' + str(cont_repetidos))
                if len(M_out[j]) != len(M_out[0]):
                    print('ERROR DE LARGOS EN M_out. Largo M_out[0]=' + str(len(M_out[0])) + '. Corrigiendo M_out[j]: ')
                    if len(M_out[j]) < len(M_out[0]):
                        M_out[j].append('');
                        print(M_out[j])
                    else:
                        print()
                        print('##########################################')
                        print('ERROR GRAVISIMO: M_out[0] puede estar mal.')
                        print(M_out[0])
                        print('##########################################')
                        print()
                # Agrego 'INSECT' al final de M_out[j]
                M_out[j][-1] = 'INSECT';
            j = j + 1;
        ### No agrego valores de M_INSECT no confirmados por M_SUG o M_SUP
        # Si no esta presente curr_SU en M_out, lo agrego
        #if not presente_en_M_out:
        #    M_out.append(curr_SU[:]);
        # Display
        if verbose and k%n_verbose == 0:
                print('Revisando ' + str(k) + ' / ' + str(len(M_INSECT)))
    if verbose:
        if insect2:
            print('curr_SU repetidos encontrados en M_INSECT: ' + str(cont_repetidos))
            print()
        print('Comparaciones para M_out terminadas. Ejemplos M_out: ')
        shuffle(M_out);
        for i in M_out[:20]:
            print(i)
        print()
    return M_out


def _main():
    # Funcion para probar funciones en ejecucion del archivo

    # Inicializo la variable que se devuelve
    M_out = [];
    M_mm9 = [];
    M_hg19 = [];

    # Uso pipeline_main() para generar M_mm9 y M_hg19
    # Guardo M_mm9 y M_hg19 apenas termino cada uno
    # mm9
    #insect_mm9 = False;
    #M_mm9 = pipeline_main(dir_SUG_main_mm9, SUG_main_mm9, lista_nombres_peaks_mm9, path_archivos_peak,
    #                      sitios_union_nkx25, sep_csv, mm9, 'mm9', insect_mm9, dir_INSECT_main,
    #                      nombre_resultados_INSECT, verbose=verbose_main, n_verbose=n_main_display);
    #M_mm9 = guardar_csv(nombre_out_mm9, M_mm9, ext_out='.csv', sep_out=sep_csv);
    # hg19
    insect_hg19 = True;
    M_hg19 = pipeline_main(dir_SUG_main_hg19, SUG_main_hg19, lista_nombres_peaks_hg19, path_archivos_peak,
                           sitios_union_nkx25, sep_csv, hg19, 'hg19', insect_hg19, dir_INSECT_main,
                           nombre_resultados_INSECT, verbose=verbose_main, n_verbose=n_main_display);
    M_hg19 = guardar_csv(nombre_out_hg19, M_hg19, ext_out='.csv', sep_out=sep_csv);
    
    # Devuelvo M_mm9, M_hg19
    M_out = (M_mm9, M_hg19);
    return M_out


'''
devolver_SUG() terminado. Ejemplos M_SUG: 
['ATAAGTG', 'ENSG00000198331', 'HYLS1', 'chr11', 125882476, 125882482]
['TTAAGTG', 'ENSG00000261821', 'AC090826.1', 'chr15', 74365115, 74365121]
['TCAAGTG', 'ENSG00000224167', 'LINC01357', 'chr1', 112848634, 112848640]
['TAAAGTG', 'ENSG00000201613', 'RNU6-200P', 'chr6', 126591089, 126591095]
['CTAAGTG', 'ENSG00000131386', 'GALNT15', 'chr3', 16173531, 16173537]
['CTAAGTG', 'ENSG00000261265', 'AC025271.2', 'chr15', 57529456, 57529462]
['TGAAGTG', 'ENSG00000165309', 'ARMC3', 'chr10', 22927117, 22927123]
['ATAAGTG', 'ENSG00000260338', 'LINC01570', 'chr16', 5617320, 5617326]
['GAAAGTG', 'ENSG00000247853', 'AC006064.2', 'chr12', 6577681, 6577687]
['ATAAGTG', 'ENSG00000287376', 'AL359538.3', 'chr13', 24543020, 24543026]
['TAAAGTG', 'ENSG00000213657', 'RPL31P44', 'chr10', 52390692, 52390698]
['GGAAGTG', 'ENSG00000283093', 'CENPVL2', 'chrX', 51684191, 51684197]
['CTAAGTG', 'ENSG00000266545', 'MIR5701-1', 'chr15', 20938780, 20938786]
['TAAAGTG', 'ENSG00000259983', 'AC018845.2', 'chr16', 47013116, 47013122]
['TCAAGTG', 'ENSG00000213280', 'AC090114.1', 'chr7', 128570924, 128570930]
['GGAAGTG', 'ENSG00000104972', 'LILRB1', 'chr19', 54617135, 54617141]
['GGAAGTG', 'ENSG00000168014', 'C2CD3', 'chr11', 74171875, 74171881]
['ATAAGTG', 'ENSG00000234737', 'KRT18P15', 'chr3', 32258402, 32258408]
['TGAAGTG', 'ENSG00000242326', 'MCUR1P2', 'chr3', 166182168, 166182174]
['GGAAGTG', 'ENSG00000270977', 'AC015849.5', 'chr17', 35913138, 35913144]

devolver_SUP() terminado. Ejemplos M_SUP: 
['TTAAGTG', '', '', 'chr3', 64239753, 64239759, 'X']
['GCAAGTG', '', '', 'chr16', 75871875, 75871881, 'X']
['GAAAGTG', '', '', 'chr2', 142280658, 142280664, 'X']
['TCAAGTG', '', '', 'chr6', 9436018, 9436024, 'X']
['TCAAGTG', '', '', 'chr6', 14623521, 14623527, 'X']
['TCAAGTG', '', '', 'chr8', 11807395, 11807401, 'X']
['TCAAGTG', '', '', 'chr9', 18662765, 18662771, 'X']
['GCAAGTG', '', '', 'chr5', 19624868, 19624874, 'X']
['CTAAGTG', '', '', 'chr4', 14869752, 14869758, 'X']
['GTAAGTG', '', '', 'chr1', 26420925, 26420931, 'X']
['GAAAGTG', '', '', 'chr3', 21479794, 21479800, 'X']
['GAAAGTG', '', '', 'chr2', 85999824, 85999830, 'X']
['GTAAGTG', '', '', 'chr9', 11417863, 11417869, 'X']
['GAAAGTG', '', '', 'chr1', 214713707, 214713713, 'X']
['CTAAGTG', '', '', 'chr1', 54712073, 54712079, 'X']
['TGAAGTG', '', '', 'chr3', 104095282, 104095288, 'X']
['TAAAGTG', '', '', 'chr2', 50040060, 50040066, 'X']
['TTAAGTG', '', '', 'chr8', 131034355, 131034361, 'X']
['TAAAGTG', '', '', 'chr1', 11908328, 11908334, 'X']
['ATAAGTG', '', '', 'chr12', 105062125, 105062131, 'X']

devolver_INSECT() terminado. Ejemplos M_INSECT: 
['TCAAGTG', 'ENSG00000250064', 'AC097480.1', 'chr4', 28435078, 28435084, 'INSECT']
['GGAAGTG', 'ENSG00000243167', 'RPS10P28', 'chr19', 42668527, 42668533, 'INSECT']
['TTAATTG', 'ENSG00000242810', 'MRPL42P6', 'chr3', 151475759, 151475753, 'INSECT']
['GGAAGTG', 'ENSG00000234422', 'UBE2WP1', 'chr1', 96419001, 96419007, 'INSECT']
['TTAATTG', 'ENSG00000235552', 'RPL6P27', 'chr18', 6463938, 6463932, 'INSECT']
['CCAAGTG', 'ENSG00000272788', 'AP000864.1', 'chr18', 8800022, 8800028, 'INSECT']
['TCAAGTG', 'ENSG00000207051', 'SNORA27', 'chr13', 27255448, 27255454, 'INSECT']
['TGAAGTG', 'ENSG00000264354', 'MIR3134', 'chr3', 15698196, 15698190, 'INSECT']
['TTAAGTG', 'ENSG00000253771', 'TPTE2P1', 'chr13', 24969078, 24969072, 'INSECT']
['TGAAGTG', 'ENSG00000239490', 'RPS4XP18', 'chr18', 22586487, 22586493, 'INSECT']
['GGAAGTG', 'ENSG00000223564', 'CYP4F32P', 'chr2', 94762118, 94762112, 'INSECT']
['CTAAGTG', 'ENSG00000244314', 'RN7SL36P', 'chr3', 195148329, 195148323, 'INSECT']
['TAAAGTG', 'ENSG00000260104', 'AC068338.1', 'chr15', 75250731, 75250725, 'INSECT']
['CTAAGTG', 'ENSG00000105717', 'PBX4', 'chr19', 19617701, 19617695, 'INSECT']
['TTAATTG', 'ENSG00000066933', 'MYO9A', 'chr15', 72119331, 72119325, 'INSECT']
['GGAAGTG', 'ENSG00000249360', 'TOMM40P3', 'chr5', 33502334, 33502340, 'INSECT']
['TCAAGTG', 'ENSG00000125735', 'TNFSF14', 'chr19', 6672570, 6672564, 'INSECT']
['TGAAGTG', 'ENSG00000243483', 'AC078785.3', 'chr3', 113041797, 113041791, 'INSECT']
['TTAATTG', 'ENSG00000241621', 'GOLGA2P6', 'chr10', 30374890, 30374884, 'INSECT']
['GTAAGTG', 'ENSG00000121964', 'GTDC1', 'chr2', 144334078, 144334072, 'INSECT']

############################################################### RESULTADOS INSECT ###############################################################

['ENSG00000270308', '-552', '-546', '9,035', '0,939894468', '+', 'NKX25', 'GCAAGTG']
Gene(gene_id='ENSG00000270308', gene_name='TEX101P1', biotype='processed_pseudogene', contig='X', start=99656997, end=99657381, strand='+', genome='GRCh38')
['ENSG00000270308', 'GCAAGTG', 'chrX', 99656445, 99656451, 'INSECT']
Secuencia encontrada: "GCAAGTG"

['ENSG00000170743', '-871', '-865', '7,681', '0,893233409', '+', 'NKX25', 'GGAAGTG']
Gene(gene_id='ENSG00000170743', gene_name='SYT9', biotype='protein_coding', contig='11', start=7238778, end=7469043, strand='+', genome='GRCh38')
['ENSG00000170743', 'GGAAGTG', 'chr11', 7237907, 7237913, 'INSECT']
Secuencia encontrada: "GGAAGTG"

['ENSG00000223576', '-1406', '-1400', '9,425', '0,953334507', '+', 'NKX25', 'TGAAGTG']
Gene(gene_id='ENSG00000223576', gene_name='AL355001.1', biotype='lncRNA', contig='13', start=19841827, end=19843672, strand='+', genome='GRCh38')
['ENSG00000223576', 'TGAAGTG', 'chr13', 19840421, 19840427, 'INSECT']
Secuencia encontrada: "TGAAGTG"

['ENSG00000215000', '-1179', '-1173', '7,506', '0,887202622', '+', 'NKX25', 'TTAAGTA']
Gene(gene_id='ENSG00000215000', gene_name='GAPDHP77', biotype='processed_pseudogene', contig='X', start=109341813, end=109342697, strand='+', genome='GRCh38')
['ENSG00000215000', 'TTAAGTA', 'chrX', 109340634, 109340640, 'INSECT']
Secuencia encontrada: "TTAAGTA"

['ENSG00000159224', '760', '766', '10,779', '0,999995567', '+', 'NKX25', 'TCAAGTG']
Gene(gene_id='ENSG00000159224', gene_name='GIP', biotype='protein_coding', contig='17', start=48958554, end=48968596, strand='-', genome='GRCh38')
['ENSG00000159224', 'TCAAGTG', 'chr17', 48967837, 48967831, 'INSECT']
Complemento encontrado: "TCAAGTG"

['ENSG00000214204', '954', '960', '7,681', '0,893233409', '+', 'NKX25', 'GGAAGTG']
Gene(gene_id='ENSG00000214204', gene_name='HNRNPA1P43', biotype='processed_pseudogene', contig='1', start=115856910, end=115857819, strand='-', genome='GRCh38')
['ENSG00000214204', 'GGAAGTG', 'chr1', 115856866, 115856860, 'INSECT']
Complemento encontrado: "GGAAGTG"

['ENSG00000248569', '708', '714', '7,791', '0,897024189', '+', 'NKX25', 'TCAAGTA']
Gene(gene_id='ENSG00000248569', gene_name='AC026410.3', biotype='processed_pseudogene', contig='5', start=80351021, end=80351956, strand='-', genome='GRCh38')
['ENSG00000248569', 'TCAAGTA', 'chr5', 80351249, 80351243, 'INSECT']
Complemento encontrado: "TCAAGTA"

['ENSG00000015592', '-931', '-925', '8,332', '0,915667936', '+', 'NKX25', 'ATAAGTG']
Gene(gene_id='ENSG00000015592', gene_name='STMN4', biotype='protein_coding', contig='8', start=27235323, end=27258420, strand='-', genome='GRCh38')
['ENSG00000015592', 'ATAAGTG', 'chr8', 27259352, 27259346, 'INSECT']
Complemento encontrado: "ATAAGTG"

['ENSG00000232019', '-811', '-805', '8,275', '0,913703622', '+', 'NKX25', 'TAAAGTG']
Gene(gene_id='ENSG00000232019', gene_name='AC074183.1', biotype='lncRNA', contig='7', start=84939349, end=84940245, strand='-', genome='GRCh38')
['ENSG00000232019', 'TAAAGTG', 'chr7', 84941057, 84941051, 'INSECT']
Complemento encontrado: "TAAAGTG"

['ENSG00000254397', '-1094', '-1088', '9,425', '0,953334507', '+', 'NKX25', 'TGAAGTG']
Gene(gene_id='ENSG00000254397', gene_name='AC132192.1', biotype='lncRNA', contig='11', start=9430356, end=9433486, strand='-', genome='GRCh38')
['ENSG00000254397', 'TGAAGTG', 'chr11', 9434581, 9434575, 'INSECT']
Complemento encontrado: "TGAAGTG"

['ENSG00000266743', '-1138', '-1132', '8,332', '0,915667936', '+', 'NKX25', 'ATAAGTG']
Gene(gene_id='ENSG00000266743', gene_name='AC103808.6', biotype='lncRNA', contig='18', start=76173304, end=76177738, strand='-', genome='GRCh38')
['ENSG00000266743', 'ATAAGTG', 'chr18', 76178877, 76178871, 'INSECT']
Complemento encontrado: "ATAAGTG"

#################################################################### ERRORES ####################################################################

['ENSG00000117114', '510', '516', '10,494', '0,990173999', '+', 'NKX25', 'TTAAGTG']
Gene(gene_id='ENSG00000117114', gene_name='ADGRL2', biotype='protein_coding', contig='1', start=81306147, end=81992436, strand='+', genome='GRCh38')
['ENSG00000117114', 'TTAAGTG', 'chr1', 81306657, 81306663, 'INSECT']
ERROR. Secuencia incorrecta: "TAGGTAC" // "GTACCTA"

['ENSG00000262668', '-1114', '-1108', '10,779', '0,999995567', '+', 'NKX25', 'TCAAGTG']
Gene(gene_id='ENSG00000262668', gene_name='AJ003147.2', biotype='lncRNA', contig='16', start=3188204, end=3224779, strand='+', genome='GRCh38')
['ENSG00000262668', 'TCAAGTG', 'chr16', 3187090, 3187096, 'INSECT']
ERROR. Secuencia incorrecta: "TCCTGGC" // "GCCAGGA"

['ENSG00000196914', '-1013', '-1007', '7,506', '0,887202622', '+', 'NKX25', 'TTAAGTA']
Gene(gene_id='ENSG00000196914', gene_name='ARHGEF12', biotype='protein_coding', contig='11', start=120336413, end=120489937, strand='+', genome='GRCh38')
['ENSG00000196914', 'TTAAGTA', 'chr11', 120335400, 120335406, 'INSECT']
ERROR. Secuencia incorrecta: "GATGTTT" // "AAACATC"

['ENSG00000074181', '-509', '-503', '10,779', '0,999995567', '+', 'NKX25', 'TCAAGTG']
Gene(gene_id='ENSG00000074181', gene_name='NOTCH3', biotype='protein_coding', contig='19', start=15159038, end=15200995, strand='-', genome='GRCh38')
['ENSG00000074181', 'TCAAGTG', 'chr19', 15201505, 15201499, 'INSECT']
ERROR. Secuencia incorrecta: "AGTTGGA" // "TCCAACT"

['ENSG00000172724', '353', '359', '8,751', '0,930107363', '+', 'NKX25', 'GTAAGTG']
Gene(gene_id='ENSG00000172724', gene_name='CCL19', biotype='protein_coding', contig='9', start=34689570, end=34691276, strand='-', genome='GRCh38')
['ENSG00000172724', 'GTAAGTG', 'chr9', 34690924, 34690918, 'INSECT']
ERROR. Secuencia incorrecta: "ACACTTA" // "TAAGTGT"

['ENSG00000060982', '-1107', '-1101', '9,035', '0,939894468', '+', 'NKX25', 'GCAAGTG']
Gene(gene_id='ENSG00000060982', gene_name='BCAT1', biotype='protein_coding', contig='12', start=24810024, end=24949101, strand='-', genome='GRCh38')
['ENSG00000060982', 'GCAAGTG', 'chr12', 24950209, 24950203, 'INSECT']
ERROR. Secuencia incorrecta: "CTTTGTC" // "GACAAAG"
'''


def _main_old():
    # Funcion para probar funciones en ejecucion del archivo

    # Inicializo la variable que se devuelve
    M_out = [];

    if correr_SUG:
        # Abrir archivos .csv en carpeta "sitios de union genoma"
        # SUG = Sitios de Union del Genoma
        M_SUG = devolver_SUG(dir_SUG_main, SUG_main, sep_csv, genoma_main, verbose=verbose_SUG1);

        if SUG_test:
            print('Iniciando pruebas de sitios de union en el genoma.')
            # Uso shuffle para obtener SUG_test_n sitios de union al azar
            shuffle(M_SUG)
            M_SUG_test = M_SUG[:SUG_test_n];
            # Pruebo los SUG_test_n sitios de union seleccionados
            for t in M_SUG_test:
                print(t)
                _prueba_secuencia(t[2], t[3], t[4], t[1], genoma_main, nombre_genoma_main, verbose=True);
                print()
        if verbose_main:
            print('devolver_SUG() terminado. Ejemplos M_SUG: ')
            shuffle(M_SUG);
            for i in M_SUG[:20]:
                print(i)
            print()
    else:
        M_SUG = [];
    if correr_SUP:
        if verbose_main:
            print('Iniciando busqueda de sitios de union en peaks.')
            print()
        # Buscar sitios de union en peaks del .csv en carpeta "csv con secuencias"
        # SUP = Sitios de Union en Peaks
        M_SUP = devolver_SUP(lista_nombres_peaks, path_archivos_peak, sitios_union_nkx25, L_genoma_main, L_nombre_genoma_main, hacer_pruebas=SUP_test);

        if verbose_main:
            print()
            print('devolver_SUP() terminado. Ejemplos M_SUP: ')
            shuffle(M_SUP);
            for i in M_SUP[:20]:
                print(i)
            print()
    else:
        M_SUP = [];
    if correr_INSECT:
        if verbose_main:
            print('Iniciando busqueda de sitios de union en INSECT2.')
            print()

        # Buscar sitios de union INSECT2
        M_INSECT = devolver_INSECT(dir_INSECT_main, nombre_resultados_INSECT, genoma_INSECT, nombre_genoma_INSECT, hacer_pruebas=INSECT_test);

        if verbose_main:
            print()
            print('devolver_INSECT() terminado. Ejemplos M_INSECT: ')
            shuffle(M_INSECT)
            for i in M_INSECT[:20]:
                print(i)
            print()
    else:
        M_INSECT = [];

    # Solo genero la tabla si tengo M_SUG y M_SUP
    if correr_SUG and correr_SUP:
        if verbose_main:
            print('Iniciando generacion de tabla con toda la info.')
            print()
            print('Armando base con M_SUG.')
        # Recorro M_SUG
        for k in range(len(M_SUG)):
            curr_SU = M_SUG[k];
            # Agrego lugar para la info de experimentos ChIP-seq
            for p in range(len(lista_nombres_peaks)):
                curr_SU = curr_SU + [''];
            # Agrego lugar para la info de INSECT
            if correr_INSECT:
                curr_SU = curr_SU + [''];
            # Agrego curr_SU a M_out
            M_out.append(curr_SU[:]);
        if verbose_main:
            print('Agregando info de M_SUP.')
        # Una vez agregados los sitios de M_SUG, agrego M_SUP y M_INSECT
        # Contador para sitios repetidos
        cont_repetidos = 0;
        # Recorro M_SUP
        for k in range(len(M_SUP)):
            curr_SU = M_SUP[k];
            # Contadores para revisar presencia de curr_SU en M_out
            presente_en_M_out = False;
            j = 0;
            # Reviso si curr_SU ya esta en M_out
            while j<len(M_out) and (not presente_en_M_out):
                # Si curr_SU esta en M_out, agrego la info entre M_out[j][5] y M_out[j][-2]
                if (str(curr_SU[2]) == str(M_out[j][2])) and ((str(curr_SU[3]) == str(M_out[j][3])) and (str(curr_SU[4]) == str(M_out[j][4]))):
                    presente_en_M_out = True;
                    cont_repetidos = cont_repetidos + 1;
                    if SU_repetido_test and cont_repetidos%n_SU_repetidos == 0:
                        print('curr_SU repetidos encontrados en M_SUP: ' + str(cont_repetidos))
                    # Recorro lista_nombres_peaks para definir en donde se encuentra curr_SU
                    for h in range(len(lista_nombres_peaks)):
                        # Registro los valores no vacios
                        if curr_SU[4+h] != '':
                            M_out[j][4+h] = str(curr_SU[4+h]);
                j = j + 1;
            # Si no esta presente curr_SU en M_out, lo agrego
            if not presente_en_M_out:
                M_out.append(curr_SU[:]);
            # Display
            if verbose_main and k%n_main_display == 0:
                    print('Revisando ' + str(k) + ' / ' + str(len(M_SUP)))
        if verbose_main:
            print('curr_SU repetidos encontrados en M_SUP: ' + str(cont_repetidos))
            print()
            print('Agregando info de M_INSECT.')
        # Contador para sitios repetidos
        cont_repetidos = 0;
        # Recorro M_INSECT
        for k in range(len(M_INSECT)):
            curr_SU = M_INSECT[k];
            # Contadores para revisar presencia de curr_SU en M_out
            presente_en_M_out = False;
            j = 0;
            # Reviso si curr_SU ya esta en M_out
            while j<len(M_out) and (not presente_en_M_out):
                # Si curr_SU esta en M_out, agrego la info en M_out[j][-1]
                if (str(curr_SU[2]) == str(M_out[j][2])) and ((str(curr_SU[3]) == str(M_out[j][3])) and (str(curr_SU[4]) == str(M_out[j][4]))):
                    presente_en_M_out = True;
                    cont_repetidos = cont_repetidos + 1;
                    if SU_repetido_test and cont_repetidos%n_SU_repetidos == 0:
                        print('curr_SU repetidos encontrados en M_INSECT: ' + str(cont_repetidos))
                    # Agrego 'INSECT' al final de M_out[j]
                    M_out[j][-1] = 'INSECT';
                j = j + 1;
            ### CORRECCION: No agrego valores de M_INSECT no confirmados por M_SUG o M_SUP
            # Si no esta presente curr_SU en M_out, lo agrego
            #if not presente_en_M_out:
            #    M_out.append(curr_SU[:]);
            # Display
            if verbose_main and k%n_main_display == 0:
                    print('Revisando ' + str(k) + ' / ' + str(len(M_INSECT)))
        if verbose_main:
            print('curr_SU repetidos encontrados en M_INSECT: ' + str(cont_repetidos))
            print()
            print('Comparaciones para M_out terminadas. Ejemplos M_out: ')
            shuffle(M_out)
            for i in M_out[:20]:
                print(i)
            print()
        ################ HACER ################
        # Guardar M_out en formato .csv
        # Titulos: 'gene_id', 'seq', 'chr', 'ini', 'end', nombres en lista_nombres_peaks, 'INSECT'

    return M_out


def _main_test():
    #

    #
    M_out = [];

    print('Iniciando busqueda de sitios de union en INSECT2.')
    print()
    # Prueba devolver_INSECT()
    M_out = devolver_INSECT(dir_INSECT_main, nombre_resultados_INSECT, genoma_INSECT, nombre_genoma_INSECT, hacer_pruebas=INSECT_test);

    print()
    print('devolver_INSECT() terminado. Ejemplos M_INSECT: ')
    shuffle(M_out)
    for i in M_out[:20]:
        print(i)
    return M_out


def _prueba_secuencia(chr_n, pos_ini, pos_end, seq_correcta, genoma, nombre_genoma, verbose=False, tiempo_sleep=2):
    # Recibe una posicion en formato chr_n, pos_ini, pos_end y revisa esa secuencia en genoma
    # Devuelve si la secuencia obtenida es igual a seq_correcta (o a su complemento)
    # verbose permite mostrar los resultados de la prueba en la consola

    # Inicializo el booleano que se devuelve
    r = False;
    
    # Saco chromosome ID con la funcion IDchr()
    chr_id = IDchr(chr_n[3:], genome=nombre_genoma);
    # Saco la secuencia correspondiente a pos_ini y pos_end
    seq_encontrada = str(ConsultaSecuencia(chr_id, pos_ini, pos_end, strand=1, sleep_time=tiempo_sleep));

    if seq_encontrada == seq_correcta:
        r = True;
        if verbose:
            print('Secuencia encontrada: "' + str(seq_encontrada) + '"')
    elif complemento_secuencia(seq_encontrada) == seq_correcta:
        r = True;
        if verbose:
            print('Complemento encontrado: "' + str(complemento_secuencia(seq_encontrada)) + '"')
    else:
        print('ERROR. Secuencia incorrecta: "' + seq_encontrada + '" // "' + complemento_secuencia(seq_encontrada) + '"')
        r = False;
    return r


#################################### RESULTADOS ###################################

output_dump = '';

if __name__=='__main__':
    output_dump = _main();
    #output_dump = _main_test();

