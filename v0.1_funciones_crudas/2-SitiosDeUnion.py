
import os
import time
import copy
from pyensembl import EnsemblRelease
from Bio import Entrez, SeqIO
Entrez.email = 'ekolomenski@gmail.com';
Entrez.api_key = '9bccbd4949f226a5930868df36c211e8b408';

# Para testing
from referencias_funciones import ConsultaSecuencia, IDchr

mm9 = EnsemblRelease(67, species='mouse');
hg19 = EnsemblRelease(102, species='human');
#

'''
~ Busca sitios de union en las secuencias de los archivos .csv
    ~ Funcion para buscar secuencias
        X Toma secuencias especificas
        - Toma matrices PWM?
            * Se puede hacer usando una lista de secuencias que pueden seleccionarse del PWM
    X Funcion para buscar secuencia complementaria
        X Vuelve a correr funcion de busqueda con complemento
    X Funcion que abre .csv y busca sitios de union en las secuencias de cada fila

~ Busca genes cerca/alrededor de cada sitio de union
'''

#################################### VARIABLES ####################################

'''
>Nombres de archivos:
Anderson2018.csv
Dupays2015.csv
Dupays2015_s1e14.csv
Dupays2015_s4e14.csv
He2011.csv
vandenBoogard2012.csv
'''

nombre_archivo = 'Dupays2015';

lista_nombres = ['Anderson2018', 'Dupays2015', 'Dupays2015_s1e14',
                 'Dupays2015_s4e14', 'He2011', 'vandenBoogard2012'];

L_genomes = [hg19, mm9, mm9, mm9, mm9, mm9];

extension = '.csv';
csv_sep = ';';

sitios_de_union_estricto = ['TTAAGTG'];
sitios_de_union_generalizado1 = ['AAGTG'];
sitios_de_union_generalizado2 = ['TNNAGTG'];
sitios_de_union_generalizado3 = ['GCAAGTG', 'GGAAGTG', 'GAAAGTG', 'ATAAGTG', 'GTAAGTG', 'CTAAGTG',
                                 'TCAAGTG', 'TGAAGTG', 'TAAAGTG', 'TTAAGTG'];
sitios_de_union_PWM = ['']; ### FALTA IMPLEMENTAR

# Variables _main()

sitios_de_union_main = sitios_de_union_generalizado3;

path_archivos = '.\\csv con secuencias\\';
curr_path = os.path.dirname(__file__);


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


def buscar_en_csv(dir_arch, nom_arch, ext, sep_arch, L_seq_busq, busq_name, dist_genes, genome, sep_out=';', seq_pos=3, print_range=500):
    # Abre un archivo .csv con secuencias y busca seq_busq en la posicion seq_pos de cada fila
    # Devuelve lista con los genes que se encontraron cerca de los sitios de union
    # Guarda un archivo con todos los sitios de union y los genes encontrados cerca
    
    # Primero mido el largo del archivo
    l_archivo = largo_archivo(dir_arch, nom_arch, ext);

    # Lista de genes buscados para devolver al final de la funcion
    L_genes = [];
    
    # Abro el archivo para leer por fila
    with open(dir_arch, 'r') as F:
        # Creo el archivo de output vacio o elimino cualquier archivo con el mismo nombre
        with open(str(nom_arch) + '_sitios_de_union_' + str(busq_name) + ext, 'w') as F_out:
            print('Archivo ' + str(nom_arch) + '_sitios_de_union_' + str(busq_name) + ext + ' creado.')
        # Lo abro en modo append
        with open(str(nom_arch) + '_sitios_de_union_' + str(busq_name) + ext, 'a') as F_out:
            # Contador para display
            k = 0;
            # Leo cada linea del .csv
            for curr_line in F:
                k += 1;
                # Transformo cada linea en una lista por columna
                L = curr_line.rstrip().split(sep_arch);
                # Registro los sitios de union
                for seq_busq in L_seq_busq:
                    # Cargo la lista de sitios de union en L_sitios
                    L_sitios = buscar_en_seq_ambas_direcciones(seq_busq, str(L[4]));
                    for sitio_union in L_sitios:
                        # Busco genes alrededor del sitio de union
                        nearby_genes = buscar_genes(L[1],L[2],L[3],dist_genes,genome);
                        for near_gene in nearby_genes:
                            # Paso los primeros items a formato de texto para guardar en .csv de output
                            r = str(L[0]) + str(sep_out) + str(L[1]) + str(sep_out) + str(L[2]) + str(sep_out) + str(L[3]);
                            # Cargo el sitio de union y la secuencia correspondiente al texto
                            r = r + str(sep_out) + str(sitio_union) + str(sep_out);
                            if sitio_union >= 0:
                                r = r + L[4][sitio_union:sitio_union+len(seq_busq)];
                            else:
                                r = r + L[4][sitio_union-len(seq_busq):sitio_union];
                            for info_gene in near_gene:
                                r = r + str(sep_out) + str(info_gene);
                            # Registro el gen en L_genes si no estaba
                            near_gene_name = (near_gene[0],near_gene[1]);
                            if not(near_gene_name in L_genes):
                                L_genes.append(copy.copy(near_gene_name));
                            # Cierro la fila y guardo cada gen en F_out
                            r = r + '\n';
                            F_out.write(r);
                        # Si no habia genes cerca del sitio de union, anoto solo el sitio
                        if len(nearby_genes) == 0:
                            # Paso los primeros items a formato de texto para guardar en .csv de output
                            r = str(L[0]) + str(sep_out) + str(L[1]) + str(sep_out) + str(L[2]) + str(sep_out) + str(L[3]);
                            # Cargo el sitio de union y la secuencia correspondiente al texto
                            r = r + str(sep_out) + str(sitio_union) + str(sep_out);
                            if sitio_union >= 0:
                                r = r + L[4][sitio_union:sitio_union+len(seq_busq)];
                            else:
                                r = r + L[4][sitio_union-len(seq_busq):sitio_union];
                            r = r + '\n';
                            F_out.write(r);
                # Display para progreso
                if k%print_range == 0:
                    print('Lineas leidas archivo ' + str(nom_arch) + ': ' + str(k) + '/' + str(l_archivo))
    return L_genes


def buscar_genes(chromosome, pos_ini, pos_end, rango_dist, genome):
    # Busca genes alrededor de "pos_ini" y "pos_end" en el cromosoma "chromosome" del genoma "genome"
    # Devuelve lista de tuplas con identificador de genes, posiciones iniciales, finales y orientaciones de los genes
    L_ret = [];
    # Busco todos los genes en el rango
    if pos_ini <= pos_end:
        gene_names_raw = genome.genes_at_locus(chromosome, int(pos_ini)-rango_dist, end=int(pos_end)+rango_dist);
    else:
        print('ADVERTENCIA: pos_ini deberia ser menor o igual a pos_end.')
        gene_names_raw = genome.genes_at_locus(chromosome, int(pos_end)-rango_dist, end=int(pos_ini)+rango_dist);
    # Selecciono dentro de los genes del rango
    for near_gene in gene_names_raw:
        L_ret.append((near_gene.gene_id, near_gene.gene_name, near_gene.biotype, near_gene.start, near_gene.end, near_gene.strand))
    ### Ver en "OLD\4-Histogramas distancias\HistogramaDistanciasFIMO.py"
    return L_ret


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


def largo_archivo(dir_arch, nom_arch, ext):
    # Devuelve el largo (en filas) de un archivo en la direccion dir_arch
    # nom_arch y ext solo sirven para display
    with open(dir_arch, 'r') as F:
        print('> Archivo ' + str(nom_arch) + str(ext) + ' abierto.')
        l_arch = len(F.read().split('\n'));
    return l_arch


def _main():
    # 
    M_out = [];

    #
    M_out = '';
    
    return M_out

#################################### RESULTADOS ###################################


output_dump = '';

if __name__=='__main__':
    output_dump = _main();


'''
for n in range(len(lista_nombres)):
    print('>>>Iniciando parseo de "' + str(lista_nombres[n]) + str(extension) + '".')
    file_path = os.path.join(curr_path, str(path_archivos) + str(lista_nombres[n]) + str(extension));
    L_out = buscar_en_csv(file_path, lista_nombres[n], extension, csv_sep, ['TTAAGTG'], 'estricto',
                          10000, L_genomes[n], sep_out=';', seq_pos=4);
    print(L_out)
    print('Parseo de "' + str(lista_nombres[n]) + str(extension) + '" finalizado.')
'''

### Test ###

# gene_names_raw = mm9.genes_at_locus('10', 15000000, end=18000000);
# print(gene_names_raw)
# Gene(gene_id='ENSMUSG00000019854', gene_name='Reps1', biotype='protein_coding', contig='10',
#       start=17775667, end=17844961, strand='+', genome='NCBIM37')

#file_path = os.path.join(curr_path, str(path_archivos) + str(nombre_archivo) + str(extension));

#L_out = buscar_en_csv(file_path, nombre_archivo, extension, csv_sep, ['TTAAGTG'], 'estricto', 10000, mm9, sep_out=';', seq_pos=4);
#print(L_out)

# seq_test = 'AACCTTCCTTAAGGACTGCCTTG';
# busq_test = 'CCTT';

'''
pos_test, encontrado_test = buscar_en_secuencia(busq_test, seq_test);
print(encontrado_test)
print(pos_test)
print(seq_test[pos_test:pos_test+len(busq_test)])
'''

'''
print(seq_test)
print(busq_test)
l_pos = buscar_en_seq_ambas_direcciones(busq_test, seq_test)
print(l_pos)
for i in l_pos:
    if i>=0:
        print(seq_test[i:i+len(busq_test)])
    else:
        print(complemento_secuencia(seq_test[i-len(busq_test):i]))
'''

