
import os
import time
import copy
from random import shuffle
from pyensembl import EnsemblRelease
from Bio import Entrez, SeqIO
Entrez.email = 'ekolomenski@gmail.com';
Entrez.api_key = '9bccbd4949f226a5930868df36c211e8b408';
from referencias_funciones import ConsultaSecuencia, IDchr

mm9 = EnsemblRelease(67, species='mouse');
hg19 = EnsemblRelease(102, species='human');
GRCm38 = EnsemblRelease(102, species='mouse');

'''
Busco sitios de union alrededor de todos los genes de un genoma

Incluye funcion de testeo que se fija que las posiciones registradas sean correctas
'''

#################################### VARIABLES ####################################

sitios_union_confirmados = ['GCAAGTG', 'GGAAGTG', 'TAAGTG', 'TCAAGTG'];
sitios_union_confirmados_v2 = ['GCAAGTG', 'GGAAGTG', 'GAAAGTG', 'ATAAGTG', 'GTAAGTG', 'CTAAGTG',
                               'TCAAGTG', 'TGAAGTG', 'TAAAGTG', 'TTAAGTG'];

nombre_output = 'SitiosDeUnionGenoma_AmbosLados';
rango_distancias = (-1000,1000);

nombre_output_random = 'SitiosDeUnionLejanos';
rango_dist_random = (-11500, -10000);

nombre_output_random_neg = 'SitiosDeUnionLejanos_neg';
rango_dist_random_neg = (-11000, -9000);

nombre_output_random_pos = 'SitiosDeUnionLejanos_pos';
rango_dist_random_pos = (9000, 11000);

genome = mm9;
genome_name = 'mm9';
L_sitios_union = sitios_union_confirmados_v2;

# Para pruebas

nombre_archivo = nombre_output;
nombre_archivo_human = nombre_output + '_hg19';
L_columnas = [];
genoma = mm9;
nombre_genoma = 'mm9';

extension = '.csv';
sep_csv = ';';
path_archivo = '.\\sitios de union genoma\\';
curr_path = os.path.dirname(__file__);

file_path = os.path.join(curr_path, str(path_archivo) + str(nombre_archivo));

#################################### FUNCIONES ####################################


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


def buscar_sitios_cerca_gen(gen_buscado, genome, genome_name, rango_dist, L_sitios_union):
    # Busca, para un gen dado, si un sitio en L_sitio_union aparece a rango_dist del inicio del gen
    # rango_dist es una tupla con una posicion inicial y final respecto a pos0 del gen

    # Lista de sitios a rango_dist del inicio del gen
    M_sitios_cerca = [];

    # Extraigo la info del cromosoma del gen
    chr_gen = IDchr(str(gen_buscado.contig),genome=genome_name);

    # Defino la posicion del +1 en base a la info que incluye el gen
    pos0, forward = buscar_pos0(gen_buscado);

    # Defino booleano por si no se encuentra ningun sitio
    sitio_encontrado = False;

    # Compruebo que se haya encontrado el cromosoma correspondiente antes de seguir
    if chr_gen != '':
        # Defino la secuencia dada por rango_dist solo si se encontro el cromosoma
        seq_promotor = secuencia_respecto_pos0(chr_gen, pos0, rango_dist[0], rango_dist[1], forward);
        # Busco cada sitio de L_sitios_union en la secuencia del promotor
        for sitio in L_sitios_union:
            # Creo lista para anotar las posiciones
            posiciones_sitio = buscar_en_seq_2dir_unificado(sitio, seq_promotor);
            # Modifico posiciones_sitio para que quede segun el +1
            for p in range(len(posiciones_sitio)):
                posiciones_sitio[p] = posiciones_sitio[p] + rango_dist[1];
            if len(posiciones_sitio) > 0:
                sitio_encontrado = True;
                L_sitios_cerca = [str(sitio), posiciones_sitio];
                M_sitios_cerca.append(L_sitios_cerca[:]);
    else:
        print('Gen salteado por estar en cromosoma ' + str(gen_buscado.contig) + ', que no se pudo parsear.')
    return M_sitios_cerca, sitio_encontrado


def buscar_sitios_cerca_genoma(genome, genome_name, rango_dist, L_sitios_union, nom_arch, ext='.csv', sep=';', print_range=50):
    # Busca, para todos los genes de un genoma, si un sitio en L_sitio_union aparece a rango_dist del inicio del gen
    # rango_dist es una tupla con una posicion inicial y final respecto a pos0 del gen
    # Los guarda en un archivo

    # Primero creo y abro el archivo de output
    with open(str(nom_arch) + str(ext), 'w') as F_out:
        print('Archivo ' + str(nom_arch) + str(ext) + ' creado.')
    with open(str(nom_arch) + str(ext), 'a') as F_out:
        # Contadores para display
        k = 0;
        m = 0;
        n_genes = len(genome.genes());
        # Busco cada uno de los genes del genoma
        for gene in genome.genes():
            # Aumento el contador de genes revisados
            k += 1;
            # Busco los sitios dentro del promotor de cada gen
            M_sitios, sitio_encontrado = buscar_sitios_cerca_gen(gene, genome, genome_name, rango_dist, L_sitios_union);
            # Solo guardo info si se encontro por lo menos un sitio
            if sitio_encontrado:
                # Aumento el contador de genes con sitios encontrados
                m += 1;
                # Agarro cada uno de los sitios de union buscados dentro de M_sitios
                for sitios_cerca_secuencia in M_sitios:
                    # Agarro la secuencia del sitio de union
                    curr_sitio = sitios_cerca_secuencia[0];
                    # Agarro la lista de posiciones donde aparece ese sitio de union
                    posiciones_curr_sitio = sitios_cerca_secuencia[1];
                    # Trabajo con cada una de las posiciones donde aparece el sitio de union curr_sitio
                    for posicion_sitio in posiciones_curr_sitio:
                        # Registro id del gen, sitio de union y posicion relativa al +1 en el archivo
                        r = str(gene.gene_id) + str(sep) + str(curr_sitio) + str(sep) + str(posicion_sitio) + '\n';
                        F_out.write(r);
                if m%print_range == 0:
                    print('Genes con sitios encontrados: ' + str(m))
            if k%print_range == 0:
                print('Genes revisados: ' + str(k) + ' de ' + str(n_genes))
    return m


def complemento(N,adn=True):
    # Devuelve el complemento de un nucleotido en adn o arn
    dict_adn = {'T':'A','U':'A','A':'T','C':'G','G':'C','N':'N'};
    dict_arn = {'T':'A','U':'A','A':'U','C':'G','G':'C','N':'N'};

    if not (N in dict_adn.keys()):
        print('Nucleotido "' + str(N) + '" no interpretado. Se devuelve N.')
        ret = 'N';
    elif adn:
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


def main():
    # Funcion que corre las funciones que voy haciendo

    genes_con_sitios = buscar_sitios_cerca_genoma(genome, genome_name, rango_distancias, L_sitios_union, nombre_output);
    print('Genes con sitios encontrados: ' + str(genes_con_sitios))
    print()
    return None


def main_human():
    # Funcion para correr buscar_sitios_cerca_genoma() con el genoma humano

    genes_con_sitios = buscar_sitios_cerca_genoma(hg19, 'hg19', rango_distancias, L_sitios_union, nombre_output + '_hg19');
    return None


def main_random():
    # Corro buscar_sitios_cerca_genoma() para sitios muy lejos de los genes

    genes_con_sitios = buscar_sitios_cerca_genoma(genome, genome_name, rango_dist_random, L_sitios_union, nombre_output_random);
    print('Genes con sitios encontrados: ' + str(genes_con_sitios))
    print()
    return None


def main_random_posneg():
    # Corro buscar_sitios_cerca_genoma() para sitios muy lejos de los genes o dentro de los genes
    genes_con_sitios_neg = buscar_sitios_cerca_genoma(genome, genome_name, rango_dist_random_neg, L_sitios_union, nombre_output_random_neg);
    print('Genes con sitios encontrados en region negativa: ' + str(genes_con_sitios_neg))
    print()
    genes_con_sitios_pos = buscar_sitios_cerca_genoma(genome, genome_name, rango_dist_random_pos, L_sitios_union, nombre_output_random_pos);
    print('Genes con sitios encontrados en region positiva: ' + str(genes_con_sitios_pos))
    print()
    return None


def main_random_human():
    # Corro buscar_sitios_cerca_genoma() para sitios muy lejos de los genes en humano

    genes_con_sitios = buscar_sitios_cerca_genoma(hg19, 'hg19', rango_dist_random, L_sitios_union, nombre_output_random + '_hg19');
    return None


def main_random_human_posneg():
    # Corro buscar_sitios_cerca_genoma() para sitios muy lejos de los genes o dentro de los genes
    genes_con_sitios_neg = buscar_sitios_cerca_genoma(hg19, 'hg19', rango_dist_random_neg, L_sitios_union, nombre_output_random_neg + '_hg19');
    print('Genes con sitios encontrados en region negativa: ' + str(genes_con_sitios_neg))
    print()
    genes_con_sitios_pos = buscar_sitios_cerca_genoma(hg19, 'hg19', rango_dist_random_pos, L_sitios_union, nombre_output_random_pos + '_hg19');
    print('Genes con sitios encontrados en region positiva: ' + str(genes_con_sitios_pos))
    print()
    return None


def main_random_human_test():
    # Funcion para probar los resultados de main_random_human()

    _prueba_sitios_genoma(nombre_output_random + '_hg19', nombre_output_random + '_hg19', L_columnas, extension, sep_csv, hg19, 'hg19');
    return None


def main_test():
    # Funcion para probar los resultados de main()
    
    _prueba_sitios_genoma(nombre_archivo, nombre_archivo, L_columnas, extension, sep_csv, genoma, nombre_genoma);
    return None


def main_test_human():
    # Funcion para probar los resultados de main_human()

    _prueba_sitios_genoma(nombre_archivo_human, nombre_archivo_human, L_columnas, extension, sep_csv, hg19, 'hg19');
    return None


def main_test_random():
    # Funcion para probar los resultados de main_random()
    
    _prueba_sitios_genoma(nombre_output_random, nombre_output_random, L_columnas, extension, sep_csv, genoma, nombre_genoma);
    return None


def secuencia_respecto_pos0(chr_gen, pos0, d_ini, d_end, forward):
    # Devuelve una secuencia respecto al +1 de un gen dado un id de cromosoma y posiciones relativas
    if forward:
        seq_buscada = str(ConsultaSecuencia(chr_gen, pos0+d_ini, pos0+d_end, strand=1));
    else:
        seq_buscada = complemento_secuencia(str(ConsultaSecuencia(chr_gen, pos0-d_end, pos0-d_ini, strand=1)));
    return seq_buscada


def _extraer_columnas_csv(dir_archivo, nom_archivo, L_col, ext='.csv', sep_arch=';'):
    # Abre un archivo .csv y devuelve las columnas en L_col como matriz
    # Si L_col es una lista vacia, se devuelven todas las columnas

    # Creo la matriz que se devuelve
    M_out = [];
    # Abro el archivo .csv
    with open(str(dir_archivo) + str(ext), 'r') as F:
        print('Archivo ' + str(nombre_archivo) + str(ext) + ' abierto.')
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
            M_out.append(L_out[:])
    return M_out


def _prueba_sitios_genoma(dir_arch, nom_arch, L_col, ext, sep_csv, genome, genome_name, modifier=0, n_test=20):
    # Prueba que las secuencias encontradas por buscar_sitios_cerca_genoma() sean correctas
    # n_test define la cantidad de filas que pruebo
    
    # Primero extraigo los resultados del .csv y los guardo en una matriz
    M_out = _extraer_columnas_csv(dir_arch, nom_arch, L_col, ext=ext, sep_arch=sep_csv);

    # Selecciono n_test filas al azar de M_out para confirmar secuencia
    shuffle(M_out)
    M_test = M_out[:n_test];

    # Contador para display
    k = 0;
    # Reviso cada uno de los sitios seleccionados
    for sitio in M_test:
        print('Sitio: ' + str(sitio))
        # Extraigo la secuencia de largo len(sitio[1]) a distancia sitio[2] de pos0 del gen con id sitio[0]
        secuencia_extraida = _secuencia_gen_por_distancia(sitio[0], int(sitio[2])+int(modifier), len(sitio[1]), genome, genome_name);
        if sitio[1] == secuencia_extraida or sitio[1] == complemento_secuencia(secuencia_extraida):
            print('OK!')
        else:
            print('Secuencia extraida: ' + str(secuencia_extraida))
        print()
    return None


def _secuencia_gen_por_distancia(gen_id, distancia, largo, genome, genome_name):
    # Extrae una secuencia de largo a distancia de pos0 del gen con gen_id en genome
    # Devuelve la secuencia como texto

    # Saco el gen del genoma en base a gen_id
    gen_buscado = genome.gene_by_id(gen_id);
    
    # Extraigo la info del cromosoma del gen
    chr_gen = IDchr(str(gen_buscado.contig),genome=genome_name);
    
    # Defino la posicion del +1 en base a la info que incluye el gen
    pos0, forward = buscar_pos0(gen_buscado);

    if forward:
        s = ConsultaSecuencia(chr_gen, pos0+distancia+1, pos0+distancia+largo, strand=1);
    else:
        s = ConsultaSecuencia(chr_gen, pos0-distancia-largo, pos0-distancia-1, strand=1);
    return s


#################################### RESULTADOS ###################################


if __name__=='__main__':
    # mm9
    print('Corriendo main()')
    main();
    print('Probando main_test()')
    main_test();
    ##main_random();
    ##print('Probando main_test_random()')
    ##main_test_random();
    print('Corriendo main_random_posneg()')
    main_random_posneg();
    # hg19
    print('Corriendo main_human()')
    main_human();
    print('Probando main_test_human()')
    main_test_human();
    ##main_random_human();
    ##print('Probando main_random_human_test()')
    ##main_random_human_test();
    print('Corriendo main_random_human_posneg()')
    main_random_human_posneg();


