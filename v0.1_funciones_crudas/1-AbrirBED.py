
import os
import time
from pyensembl import EnsemblRelease
from Bio import Entrez, SeqIO
Entrez.email = 'ekolomenski@gmail.com';
Entrez.api_key = '9bccbd4949f226a5930868df36c211e8b408';

'''
Agarra archivos .bed con resultados de ChIP-seq

Agrega las secuencias correspondientes a cada pico de ChIP-seq

Guarda todo en un archivo .csv
'''


#################################### VARIABLES ####################################

'''
>Nombres de archivos:
Anderson2018-GSE89457consensus.bed
Dupays2015.bed
Dupays_GSM1087143_nkx_s1e14_peaks.bed
Dupays_GSM1087144_nkx_s4e14_peaks.bed
He2011.bed
vandenBoogard2012.bed
'''

nombre_archivo = 'Dupays2015';
extension = '.bed';
genoma_bed = 'mm9'; # Creo que todos menos anderson son mm9

lista_nombres = ['Anderson2018-GSE89457consensus', 'Dupays2015', 'Dupays_GSM1087143_nkx_s1e14_peaks',
                 'Dupays_GSM1087144_nkx_s4e14_peaks', 'He2011', 'vandenBoogard2012'];
lista_genomas = ['hg19','mm9','mm9','mm9','mm9','mm9'];

path_archivos = '..\\0-Fuentes\\Papers ChIP-seq\\';
curr_path = os.path.dirname(__file__);
file_path = os.path.join(curr_path, str(path_archivos) + str(nombre_archivo) + str(extension));

csv_sep = ';';
bed_sep = '\t';


##################################### GENOMAS #####################################

mm9 = EnsemblRelease(67, species='mouse');
hg19 = EnsemblRelease(102, species='human');


#################################### FUNCIONES ####################################


def ConsultaSecuencia(id_chr, seq_start, seq_finish, strand=1):
    # Devuelve una secuencia dado un ID de cromosoma (incluye info de especie) y posicion inicial/final
    time.sleep(0.1);
    rec_seq = '';
    try:
        handle = Entrez.efetch(db='nucleotide', id=id_chr, rettype='fasta',
                               strand=strand, seq_start=seq_start, seq_stop=seq_finish);
        record = SeqIO.read(handle, 'fasta');
        handle.close();
        rec_seq = record.seq
    except:
        print('Exception raised for chr ' + str(id_chr) + ' between positions ' +
              str(seq_start) + ' and ' + str(seq_finish) + '.')
        time.sleep(60);
        try:
            handle = Entrez.efetch(db='nucleotide', id=id_chr, rettype='fasta',
                                   strand=strand, seq_start=seq_start, seq_stop=seq_finish);
            record = SeqIO.read(handle, 'fasta');
            handle.close();
            rec_seq = record.seq
        except:
            print('Retry failed. Returning empty string.')
    return rec_seq


def largo_archivo(dir_arch, nom_arch, ext):
    # Devuelve el largo (en filas) de un archivo en la direccion dir_arch
    # nom_arch y ext solo sirven para display
    with open(dir_arch, 'r') as F:
        print('> Archivo ' + str(nom_arch) + str(ext) + ' abierto.')
        l_arch = len(F.read().split('\n'));
    return l_arch


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
                      'M':'NC_012920.1'};
    elif genome.lower() == 'mm9': # MGSCv37
        dict_IDchr = {'1':'NC_000067.5', '2':'NC_000068.6', '3':'NC_000069.5', '4':'NC_000070.5',
                      '5':'NC_000071.5', '6':'NC_000072.5', '7':'NC_000073.5', '8':'NC_000074.5',
                      '9':'NC_000075.5', '10':'NC_000076.5', '11':'NC_000077.5', '12':'NC_000078.5',
                      '13':'NC_000079.5', '14':'NC_000080.5', '15':'NC_000081.5', '16':'NC_000082.5',
                      '17':'NC_000083.5', '18':'NC_000084.5', '19':'NC_000085.5', 'X':'NC_000086.6',
                      'Y':'NC_000087.6', 'M':'NC_005089.1'};
    else:
        print('No se pudo encontrar genoma ' + str(chromosome))
        b = False;
        dict_IDchr = {};

    if str(chromosome).upper() in dict_IDchr.keys():
        ret = dict_IDchr[str(chromosome).upper()];
    elif b:
        print('No se pudo encontrar cromosoma ' + str(chromosome))
    return ret

def parsear_bed(dir_arch, nom_arch, ext, sep_arch, genoma, sep_out=';', print_range=25):
    # Abre un archivo .bed y procesa cada pico con parsear_fila_bed()
    
    # Primero mido el largo del archivo
    l_archivo = largo_archivo(dir_arch, nom_arch, ext);
    L_error = [];
    
    # Abro el .bed o .csv para leer por fila
    with open(dir_arch, 'r') as F:
        # Creo el archivo de output vacio o elimino cualquier archivo con el mismo nombre
        with open(str(nom_arch) + '_output.csv', 'w') as F_out:
            print('Archivo ' + str(nom_arch) + '_output.csv creado.')
        # Lo vuelvo a abrir en modo append
        with open(str(nom_arch) + '_output.csv', 'a') as F_out:
            # Contador para display
            k = 0;
            # Leo cada linea del .bed o .csv
            for curr_line in F:
                k += 1;
                # Variable para agarrar errores
                err = []
                # Funcion para todo lo que hago por linea
                r, err = parsear_fila_bed(curr_line.rstrip().split(sep_arch), k, sep_out, genoma);
                # Guardo el error de ser necesario
                if len(err) >= 1:
                    L_error.append(err[:]);
                # Escribo en el archivo
                F_out.write(r + '\n');
                # Display para progreso
                if k%print_range == 0:
                    print('Lineas leidas archivo ' + str(nom_arch) + ': ' + str(k) + '/' + str(l_archivo))

    # Inicializo la variable que se devuelve para incluir errores no resueltos
    ret = [];
    if L_error == []:
        print('Sin errores detectados.')
    else:
        print('Iniciando reintento de errores.')
        for i in range(len(L_error)):
            print(i)
            secuencia_error = ConsultaSecuencia(IDchr(L_error[i][0],genoma), L_error[i][1], L_error[i][2])
            if secuencia_error == '':
                ret.append(L_error[i]);
            else:
                print(secuencia_error)
    return ret

def parsear_fila_bed(L_row, k, sep, genoma):
    # Procesa una fila de un archivo .bed o .csv en formato lista y la devuelve como fila de .csv
    
    # Inicializo las variables que se devuelven al final
    ret = '';
    err = []
    
    # Extraigo cromosoma, posicion inicial y posicion final de L_row
    chromosome = L_row[0][3:];
    pos_ini = min(int(L_row[1]), int(L_row[2]));
    pos_end = max(int(L_row[1]), int(L_row[2]));
    
    # Obtengo el ID del cromosoma para extraer la secuencia
    id_chr = IDchr(chromosome,genome=str(genoma));
    # Extraigo la secuencia con id_chr, pos_ini y pos_end
    if id_chr != '':
        secuencia_posicion = str(ConsultaSecuencia(id_chr, pos_ini, pos_end));
    else:
        secuencia_posicion = '';

    # Registro error si la secuencia sale mal
    if secuencia_posicion == '':
        err = [str(chromosome),int(pos_ini),int(pos_end)];
        
    # Agrego todos los datos a la variable al final
    ret = str(k) + str(sep) + str(chromosome) + str(sep) + str(pos_ini) + str(sep) + str(pos_end) + str(sep) + secuencia_posicion;    

    return ret, err

#################################### RESULTADOS ###################################

for n in range(len(lista_nombres)):
    print('Iniciando parseo de "' + str(lista_nombres[n]) + str(extension) + '".')
    file_path = os.path.join(curr_path, str(path_archivos) + str(lista_nombres[n]) + str(extension));
    M_error = parsear_bed(file_path, lista_nombres[n], extension, bed_sep, lista_genomas[n]);
    print(M_error)
    print('Parseo de "' + str(lista_nombres[n]) + str(extension) + '" finalizado.')

