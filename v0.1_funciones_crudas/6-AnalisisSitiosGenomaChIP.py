
import os
import time
import copy
from random import shuffle
from pyensembl import EnsemblRelease
from Bio import Entrez, SeqIO
Entrez.email = 'ekolomenski@gmail.com';
Entrez.api_key = '9bccbd4949f226a5930868df36c211e8b408';
from referencias_funciones import ConsultaSecuencia, IDchr, extraer_columnas_csv

mm9 = EnsemblRelease(67, species='mouse');
hg19 = EnsemblRelease(102, species='human');
GRCm38 = EnsemblRelease(102, species='mouse');

'''
Coteja los hits de sitios de union en todos los promotores del genoma con los datos de ChIP-seq

Compara los resultados con los sitios de union en regiones alejadas de los promotores

Guarda en un archivo .csv sitios de union ubicados en peaks de ChIP (guarda sitios y peaks)
'''

#################################### VARIABLES ####################################


nombres_archivos_chip_hg19 = ['Anderson2018-GSE89457consensus'];
nombres_archivos_chip_mm9 = ['Dupays2015', 'He2011', 'vandenBoogard2012'];
nombres_archivos_chip_mm9_extra = ['Dupays_GSM1087143_nkx_s1e14_peaks', 'Dupays_GSM1087144_nkx_s4e14_peaks'];
path_chip = '..\\0-Fuentes\\Papers ChIP-seq\\';
ext_chip = '.bed';

nombre_archivo_sitios_promotor_mm9 = 'SitiosDeUnionGenoma';
nombre_archivo_sitios_lejano_mm9 = 'SitiosDeUnionLejanos';
nombre_archivo_sitios_promotor_hg19 = 'SitiosDeUnionGenoma_hg19';
nombre_archivo_sitios_lejano_hg19 = 'SitiosDeUnionLejanos_hg19';
path_sitios = '.\\sitios de union genoma\\';
ext_sitios = '.csv';

curr_path = os.path.dirname(__file__);

# Variables main()
genoma_usado = mm9;
nombre_genoma_usado = 'mm9';
main_nom_arch_chip = nombres_archivos_chip_mm9_extra[1];
main_nom_arch_sitios = nombre_archivo_sitios_promotor_mm9;
m_dist = 0;
print_verbose = False;
final_print = False;

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


def abrir_resultados_csv(dir_arch, nom_arch, sep_arch=';'):
    # Abre el archivo .csv en dir_arch y devuelve todas las columnas como matriz

    # Inicializo la matriz a devolver
    M_csv = [];

    # Abro el archivo
    with open(dir_arch, 'r') as F:
        print('Archivo ' + str(nom_arch) + '.csv abierto.')
        for curr_line in F:
            L_curr_line = curr_line.rstrip().split(sep_arch);
            M_csv.append(L_curr_line[:])
    return M_csv


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


def buscar_sitios_en_chip(nombre_bed, path_bed, nombre_csv, path_csv, genome, genome_name, m=0, cont_timer=10000, verbose=False):
    # Busca sitios de union que se encuentren dentro de picos ChIP-seq
    # Sitios de union en archivo nombre_csv. Picos ChIP-seq en archivo nombre_bed
    # Devuelve dos listas coordinadas, una con los sitios de union en picos ChIP-seq y otra con los picos correspondientes

    # Abro el archivo nombre_bed y guardo el contenido en M_out_bed
    arch_path_chip = os.path.join(curr_path, path_bed + nombre_bed + '.bed');
    M_out_bed = abrir_resultados_bed(arch_path_chip, nombre_bed);
    # Abro el archivo nombre_csv y guardo el contenido en M_out_csv
    arch_path_sitios = os.path.join(curr_path, path_csv + nombre_csv + '.csv');
    M_out_csv = abrir_resultados_csv(arch_path_sitios, nombre_csv);

    # Separador para display
    print()
    
    # Inicializo las variables a devolver
    sitios_en_chip = [];
    picos_con_sitios = [];
    sitios_no_en_chip = [];
    
    # Defino el largo de M_out_csv para optimizar tiempo
    l_csv = len(M_out_csv);
    # Reviso cada uno de los sitios en M_out_csv para ver si estan cerca de algun sitio de M_out_bed
    for s in range(l_csv):
        # Agrego el modificador para sitios lejanos registrados como cercanos
        M_out_csv[s][2] = int(M_out_csv[s][2])+int(m);
        # Paso el sitio al mismo formato que los picos .bed
        curr_sitio = parsear_sitio_csv(M_out_csv[s], genome, genome_name);
        # Busco curr_sitio en todos los picos del .bed, registro hits
        hits_en_bed = buscar_sitio_en_bed(curr_sitio, M_out_bed);
        # Si el sitio aparece en algun pico .bed, lo registro en el output
        if len(hits_en_bed) > 0:
            # Registro los sitios que aparecen en picos de chip
            sitios_en_chip.append(M_out_csv[s]);
            # Registro los picos en los que aparece cada uno de los sitios
            picos_con_sitios.append(hits_en_bed[:]);
        else:
            sitios_no_en_chip.append(M_out_csv[s]);
        # Contador para display
        if (s+1)%cont_timer == 0:
            print('Progreso: ' + str(s+1) + ' / ' + str(l_csv))
            print()
            # Para prueba usar verbose = True
            if verbose:
                print('Prueba de sitio: ' + str(curr_sitio))
                test_seq_sitio = M_out_csv[s][1];
                test_gen = genome.gene_by_id(M_out_csv[s][0]);
                test_gen_contig = str(test_gen.contig);
                _prueba_posicion_sitios(test_seq_sitio, test_gen_contig, curr_sitio[1], curr_sitio[2], genome_name);
                print()
    print('Cantidad de sitios encontrados en picos de ChIP-seq: ' + str(len(sitios_en_chip)))
    print()
    return sitios_en_chip, picos_con_sitios


def buscar_sitio_en_bed(sitio, L_bed):
    # Recibe un sitio y lo busca en una lista de picos de ChIP-seq
    # Sitio y cada pico en formato [chr_n, pos_ini, pos_end]

    # Inicializo la lista que devuelvo
    L_hits = [];
    # Recorro cada uno de los picos de la lista
    for i in range(len(L_bed)):
        curr_peak = L_bed[i];
        curr_peak[1] = int(curr_peak[1]);
        curr_peak[2] = int(curr_peak[2]);
        # Primero reviso que ambas secuencias se encuentren en el mismo cromosoma
        if curr_peak[0] == sitio[0]:
            if secuencia_en_secuencia(sitio[1],sitio[2],curr_peak[1],curr_peak[2]):
                L_hits.append(curr_peak[:]);
                if (sitio[1] < curr_peak[1]) or (sitio[2] > curr_peak[2]):
                    print('Sitio encontrado en un rincon del pico de ChIP')
                    print(str(curr_peak) + ' // ' + str(sitio))
                    print()
    return L_hits


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


def guardar_dos_L_csv(M1, M2, nombre_output, sep=';'):
    # Guarda los contenidos de dos matrices en un csv
    # M1 es una lista de matrices y M2 es una lista de listas

    # Reviso que ambas listas tengan el mismo largo
    if len(M1) != len(M2):
        print('CUIDADO: Matrices de diferente largo. Se copia hasta el largo de la matriz mas corta.')

    # Inicializo una matriz conjunta para devolver
    M_out = [];
    
    # Inicializo el archivo
    with open(nombre_output + '.csv', 'w') as F_out:
        print('Archivo ' + nombre_output + '.csv creado.')
    with open(nombre_output + '.csv', 'a') as F_out:
        # Recorro ambas listas
        for i in range(min(len(M1),len(M2))):
            # La matriz M1 contiene matrices en vez de listas
            curr_M1 = [];
            for hits in M1[i]:
                curr_M1 = curr_M1 + hits;
            ###### Prueba
            if curr_M1 != M1[i][0]:
                print('M1[i][0]: ' + str(M1[i][0]))
                print('curr_M1: ' + str(curr_M1))
            # Defino la fila conjunta
            curr_L = M1[i][0] + M2[i];
            # Inicializo el texto con el id en la primer columna
            r = str(i);
            # Agrego cada dato de la fila conjunta
            for j in range(len(curr_L)):
                r = r + str(sep) + str(curr_L[j]);
            # Cierro el texto de fila con un fin de linea y lo guardo en F_out
            r = r + '\n';
            F_out.write(r);
            # Guardo curr_L en M_out
            M_out.append(curr_L[:]);
    return M_out


def main():
    # Funcion para probar funcion buscar_sitios_en_chip()
    sitios_en_chip, picos_con_sitios = buscar_sitios_en_chip(main_nom_arch_chip, path_chip, main_nom_arch_sitios,
                                                             path_sitios, genoma_usado, nombre_genoma_usado,
                                                             m=m_dist, verbose=print_verbose);

    if final_print:
        # Printeo los sitios encontrados
        for i in range(len(sitios_en_chip)):
            print(sitios_en_chip[i])
            for pico in picos_con_sitios[i]:
                print(pico)
            print()

    nombre_out = main_nom_arch_chip[:6] + '_' + main_nom_arch_sitios[13:];
    M_out = guardar_dos_L_csv(picos_con_sitios, sitios_en_chip, nombre_out);
    return M_out


def parsear_sitio_csv(L_sitio, genome, genome_name):
    # Recibe un sitio en formato [gen_id, seq, relative_pos] y lo devuelve en formato [chr_n, pos_ini, pos_end]

    # Inicializo la variable a devolver
    L_out = ['chr', -1, -1];

    # Defino el gen correspondiente al sitio
    curr_gen = genome.gene_by_id(L_sitio[0]);

    # Defino el cromosoma correspondiente
    curr_gen_contig = curr_gen.contig;
    L_contigs = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16',
                 '17','18','19','20','21','22','23','X','Y','M'];
    if str(curr_gen_contig).upper() in L_contigs:
        L_out[0] = 'chr' + str(curr_gen_contig).upper();
    elif str(curr_gen_contig).upper() == 'MT':
        L_out[0] = 'chrM';
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
    L_out[1] = pos_ini;
    L_out[2] = pos_end;
    
    return L_out


def secuencia_en_secuencia(x1,x2,y1,y2):
    # Revisa que (x1,x2) se encuentre dentro de (y1,y2)
    return ((x2 > y1) and (x1 < y2))


def _prueba_posicion_sitios(seq_sitio, gen_contig, pos_ini, pos_end, genome_name):
    # Funcion para probar si una secuencia esta bien extraida

    # Extraigo la info del cromosoma del gen
    chr_gen = IDchr(str(gen_contig),genome=genome_name);

    # Busco la secuencia correspondiente a chr_gen, pos_ini, pos_end
    seq_pos = ConsultaSecuencia(chr_gen, pos_ini, pos_end);

    if seq_pos == seq_sitio or seq_pos == complemento_secuencia(seq_sitio):
        print('OK!')
    else:
        print('Secuencia extraida: ' + str(seq_pos) + ' // Secuencia buscada: ' + str(seq_sitio))
    return None



#################################### RESULTADOS ###################################


output_dump = '';

if __name__=='__main__':
    output_dump = main();


