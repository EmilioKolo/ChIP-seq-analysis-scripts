
import copy
import os
import time
from pyensembl import EnsemblRelease
from Bio import Entrez, SeqIO
Entrez.email = 'ekolomenski@gmail.com';
Entrez.api_key = '9bccbd4949f226a5930868df36c211e8b408';
from referencias_funciones import ConsultaSecuencia, IDchr, abrir_matriz_pbm

mm9 = EnsemblRelease(67, species='mouse');
hg19 = EnsemblRelease(102, species='human');

'''
Busca genes alrededor de picos de ChIP-seq

Busca picos de ChIP-seq alrededor del inicio de genes dados

Busca sitios de union alrededor del inicio de genes dados
'''

#################################### VARIABLES ####################################


lista_nombres = ['Anderson2018', 'Dupays2015', 'Dupays2015_s1e14',
                 'Dupays2015_s4e14', 'He2011', 'vandenBoogard2012'];

L_genomes = [hg19, mm9, mm9, mm9, mm9, mm9];

extension = '.csv';
csv_sep = ';';

path_archivos = '.\\csv con secuencias\\';
curr_path = os.path.dirname(__file__);

L_mouse_old = ['ENSMUSG00000013936', 'ENSMUSG00000041616', 'ENSMUSG00000029019', 'ENSMUSG00000024803',
           'ENSMUSG00000005583', 'ENSMUSG00000037335', 'ENSMUSG00000037169', 'ENSMUSG00000021604',
           'ENSMUSG00000030046', 'ENSMUSG00000059325', 'ENSMUSG00000024966'];
L_human_old = ['ENSG00000111245', 'ENSG00000175206', 'ENSG00000120937', 'ENSG00000148677',
           'ENSG00000081189', 'ENSG00000113196', 'ENSG00000134323', 'ENSG00000113430',
           'ENSG00000163217', 'ENSG00000171476', 'ENSG00000168439'];

sitio_union_canonico = ['TTAAGTG'];
sitio_union_secundario = ['CATTAATTN','CTTTAATTN','GTTAATTC'];
sitios_union_PWM_largo = ['AAGTACTTAA','AAGCACTTAA','AAGTGCTTAA','AAGTACTCAA'];
sitios_union_PWM_corto = ['AAGTACTT','AAGCACTT','AAGTGCTT','AAGTACTC'];
sitios_union_PWM_largo_corto = ['AAGTACTTAA','AAGCACTTAA','AAGTGCTTAA','AAGTACTCAA','AAGTACTT','AAGCACTT','AAGTGCTT','AAGTACTC'];

sitios_union_confirmados = ['GCAAGTG', 'GGAAGTG', 'TAAGTG', 'TCAAGTG'];

rango_distancias = 40000;
genome = mm9;
genome_name = 'mm9';
L_genes = L_mouse_old;
L_sitios_union = sitios_union_confirmados;


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


def buscar_genes_picos_csv(dir_arch, nom_arch, ext, sep_arch, dist_genes, genome, sep_out=';', print_range=500):
    # Abre un archivo .csv con posiciones y busca genes alrededor de las posiciones
    # Devuelve una lista de los genes encontrados
    # Guarda un archivo con todos los sitios de union y las posiciones de los genes

    # Primero mido el largo del archivo
    l_archivo = largo_archivo(dir_arch, nom_arch, ext);

    # Lista de genes buscados para devolver al final de la funcion
    L_genes = [];

    # Abro el archivo para leer por fila
    with open(dir_arch, 'r') as F:
        # Creo el archivo de output vacio o elimino cualquier archivo con el mismo nombre
        with open(str(nom_arch) + '_genes_cercanos' + ext, 'w') as F_out:
            print('Archivo ' + str(nom_arch) + '_genes_cercanos' + ext + ' creado.')
        # Lo abro en modo append
        with open(str(nom_arch) + '_genes_cercanos' + ext, 'a') as F_out:
            # Contador para display
            k = 0;
            # Leo cada linea del .csv
            for curr_line in F:
                k += 1;
                # Transformo cada linea en una lista por columna
                L = curr_line.rstrip().split(sep_arch);
                # Busco genes alrededor del pico
                nearby_genes = buscar_genes(L[1],L[2],L[3],dist_genes,genome);
                for near_gene in nearby_genes:
                    # Paso los primeros items a formato de texto para guardar en .csv de output
                    r = str(L[0]) + str(sep_out) + str(L[1]) + str(sep_out) + str(L[2]) + str(sep_out) + str(L[3]);
                    for info_gene in near_gene:
                        r = r + str(sep_out) + str(info_gene);
                    # Registro el gen en L_genes si no estaba
                    near_gene_name = (near_gene[0],near_gene[1]);
                    if not(near_gene_name in L_genes):
                        L_genes.append(copy.copy(near_gene_name));
                    # Cierro la fila y guardo cada gen en F_out
                    r = r + '\n';
                    F_out.write(r);
                # Si no habia genes cerca del pico, anoto solo el sitio
                if len(nearby_genes) == 0:
                    # Paso los primeros items a formato de texto para guardar en .csv de output
                    r = str(L[0]) + str(sep_out) + str(L[1]) + str(sep_out) + str(L[2]) + str(sep_out) + str(L[3]);
                    r = r + '\n';
                    F_out.write(r);
                # Display para progreso
                if k%print_range == 0:
                    print('Lineas leidas archivo ' + str(nom_arch) + ': ' + str(k) + '/' + str(l_archivo))
    return L_genes


def buscar_picos(chromosome, pos_ini, pos_end, rango_dist, dir_arch, nom_arch, ext, sep_arch):
    # Busca picos en un .csv que sean cercanos a pos_ini y pos_end
    # Idealmente se usa el extremo inicial del gen, con pos_ini = pos_end
    # Devuelve lista de tuplas con identificador del pico de ChIP-seq

    # Primero mido el largo del archivo
    l_archivo = largo_archivo(dir_arch, nom_arch, ext, verbose=False);

    # Lista de picos encontrados cerca del gen
    L_picos = [];

    # Reviso que pos_ini sea menor o igual a pos_end
    if pos_ini <= pos_end:
        y_ini = int(pos_ini)-rango_dist;
        y_end = int(pos_end)+rango_dist;
    else:
        y_end = int(pos_ini)-rango_dist;
        y_ini = int(pos_end)+rango_dist;

    # Abro el archivo para leer por fila
    with open(dir_arch, 'r') as F:
        # Contador para display
        k = 0;
        # Leo cada linea del .csv
        for curr_line in F:
            k += 1;
            # Transformo cada linea en una lista por columna
            L = curr_line.rstrip().split(sep_arch);

            if chromosome == L[1] and esta_en_rango(int(L[2]),int(L[3]),y_ini,y_end):
                r = (str(L[0]), str(L[1]), str(L[2]), str(L[3]));
                L_picos.append(copy.copy(r));
    return L_picos


def buscar_picos_cerca_de_genes(L_genes_busq, genome, rango_dist, dir_arch, nom_arch, ext, sep_arch, sep_out=';', extension_rango=2):
    # Busca una lista de genes L_genes_busq en el archivo dir_arch

    # Creo una lista para recopilar genes no encontrados
    L_error = [];
    # Creo un archivo de output para guardar los picos
    with open(str(nom_arch) + '_picos_cerca_de_genes' + ext, 'w') as F_out:
        print('Archivo ' + str(nom_arch) + '_picos_cerca_de_genes' + ext + ' creado.')
    # Lo abro en modo append
    with open(str(nom_arch) + '_picos_cerca_de_genes' + ext, 'a') as F_out:
        # Contador para display
        k = 0;
        # Busco picos cercanos a cada gen con buscar_picos()
        for gen_id in L_genes_busq:
            gen_busq = genome.gene_by_id(gen_id);
            chr_gen = gen_busq.contig;
            print('Iniciando busqueda cerca de gen ' + str(gen_busq.gene_name) + ' en cromosoma ' + str(chr_gen) + '.')
            # Uso ini_gen para strand '+' y end_gen para strand '-'
            if gen_busq.strand == '+':
                ini_gen = gen_busq.start;
                end_gen = gen_busq.start;
            elif gen_busq.strand == '-':
                ini_gen = gen_busq.end;
                end_gen = gen_busq.end;
            else:
                print('ERROR PARSEANDO STRAND DE GEN:')
                print(gen_busq)
                ini_gen = gen_busq.start;
                end_gen = gen_busq.end;
            # Uso buscar_picos() para encontrar los picos cerca de cada gen
            L_picos_gen = buscar_picos(chr_gen, ini_gen, end_gen, rango_dist, dir_arch, nom_arch, ext, sep_arch);
            # Chequeo si se encontraron genes en el rango dado
            if len(L_picos_gen) == 0:
                print('No se encontraron picos con el rango dado cerca de gen:')
                print(gen_busq)
                rango_extendido = rango_dist;
                print('Probando extension de rango')
                # Pruebo aumentar el rango hasta llegar a mil millones de bases
                while (len(L_picos_gen) == 0) and rango_extendido < 1000000000:
                    rango_extendido = rango_extendido * extension_rango;
                    L_picos_gen = buscar_picos(chr_gen, ini_gen, end_gen, rango_extendido, dir_arch, nom_arch, ext, sep_arch);
                if len(L_picos_gen) == 0:
                    print('ERROR: No se encontraron picos con el rango extendido al maximo.')
                    L_error.append((gen_busq.gene_id, gen_busq.gene_name));
                else:
                    print('Se encontraron picos con el rango extendido a ' + str(rango_extendido) + '.')
            # Guardo cada pico en F_out
            for pico_cercano in L_picos_gen:
                # Inicializo la variable para cada fila con la info de gen_busq separada por sep_out
                r = str(gen_busq.gene_id) + str(sep_out) + str(gen_busq.gene_name) + str(sep_out);
                r = r + str(gen_busq.biotype) + str(sep_out) + str(gen_busq.contig) + str(sep_out);
                r = r + str(gen_busq.start) + str(sep_out) + str(gen_busq.end) + str(sep_out);
                r = r + str(gen_busq.strand) + str(sep_out);
                # Agrego cada uno de los datos del pico separados por sep_out
                for p in range(len(pico_cercano)):
                    r = r + str(pico_cercano[p]) + str(sep_out);
                # Calculo la distancia del gen al pico
                x1 = min(int(pico_cercano[2]),int(pico_cercano[3]));
                x2 = max(int(pico_cercano[2]),int(pico_cercano[3]));
                dist = distancia_rangos(x1, x2, int(ini_gen), int(end_gen));
                # Agrego la distancia al fin de la linea
                r = r + str(dist);
                # Agrego el fin de linea y lo guardo en F_out
                r = r + '\n';
                F_out.write(r);
    return L_error


def distancia_rangos(x_ini, x_end, y_ini, y_end):
    if esta_en_rango(x_ini, x_end, y_ini, y_end):
        d = 0;
    else:
        d = min(abs(x_ini-y_end),abs(x_end-y_ini));
    return d


def esta_en_rango(x_ini, x_end, y_ini, y_end):
    # Funcion que verifica que los rangos x=(x_ini,x_end) e y=(y_ini,y_end) se solapen
    b1 = x_ini <= y_end;
    b2 = y_ini <= x_end;
    return (b1 and b2)


def buscar_sitios_cerca_gen(gen_buscado, genome, genome_name, rango_dist, L_sitios_union):
    # Busca, para un gen dado, si un sitio en L_sitio_union aparece a rango_dist del inicio del gen

    # Lista de sitios a rango_dist del inicio del gen
    M_sitios_cerca = [];

    # Extraigo la info del cromosoma del gen
    chr_gen = IDchr(str(gen_buscado.contig),genome=genome_name);

    # Determino la direccion del gen y defino la posicion del 0
    if gen_buscado.strand == '+':
        forward = True;
        pos0 = gen_buscado.start;
    elif gen_buscado.strand == '-':
        forward = False;
        pos0 = gen_buscado.end;
    else:
        print('ERROR en strand de gen ' + str(gen_buscado))
        forward = True;
        pos0 = gen_buscado.start;

    # Busco las secuencias rio abajo y rio arriba del gen
    seq_rio_arriba = '';
    seq_rio_abajo = '';
    if forward:
        seq_rio_arriba = ConsultaSecuencia(chr_gen, pos0-rango_dist, pos0, strand=1);
        seq_rio_abajo = ConsultaSecuencia(chr_gen, pos0, pos0+rango_dist, strand=1);
    else:
        seq_rio_arriba = complemento_secuencia(ConsultaSecuencia(chr_gen, pos0, pos0+rango_dist, strand=1));
        seq_rio_abajo = complemento_secuencia(ConsultaSecuencia(chr_gen, pos0-rango_dist, pos0, strand=1));

    # Busco cada uno de los sitios de union
    for sitio_union in L_sitios_union:
        # Inicializo la lista con los datos de cada sitio de union
        L_sitios_cerca = [str(sitio_union)];
        L_sitios_cerca_complemento = [str(complemento_secuencia(sitio_union))];
        # Busco los sitios de union rio arriba y rio abajo con la funcion buscar_en_seq_ambas_direcciones()
        L_sitios_rio_arriba = buscar_en_seq_ambas_direcciones(str(sitio_union), seq_rio_arriba);
        L_sitios_rio_abajo = buscar_en_seq_ambas_direcciones(str(sitio_union), seq_rio_abajo);
        # Reviso los sitios de union rio arriba
        for i in range(len(L_sitios_rio_arriba)):
            # Si la secuencia rio arriba tiene un sitio de union, lo registro como numero negativo
            if L_sitios_rio_arriba[i] >= 0:
                # Si el sitio rio arriba es positivo, se agrega como sitio forward a L_sitios_cerca
                # L_sitios_rio_arriba[i] es la distancia del sitio a pos0-rango_dist
                # Se resta rango_dist a L_sitios_rio_arriba[i] para obtener la distancia del sitio a pos0
                curr_sitio = L_sitios_rio_arriba[i]-rango_dist;
                if forward:
                    ### PRUEBA:
                    #print('PRUEBA: Sitio de union buscado vs region encontrada. Sitio rio arriba positivo, gen forward.')
                    #print('Sitio de union buscado: ' + str(sitio_union))
                    #print('Sitio encontrado: ' + str(ConsultaSecuencia(chr_gen, pos0+curr_sitio, pos0+curr_sitio+len(sitio_union)-1, strand=1)))
                    L_sitios_cerca.append(curr_sitio);
                else:
                    ### PRUEBA:
                    #print('PRUEBA: Sitio de union buscado vs region encontrada. Sitio rio arriba positivo, gen reverse.')
                    #print('Sitio de union buscado: ' + str(sitio_union))
                    #print('Sitio encontrado: ' + str(ConsultaSecuencia(chr_gen, pos0-curr_sitio-len(sitio_union)+1, pos0-curr_sitio, strand=1)))
                    L_sitios_cerca.append(-(curr_sitio+len(sitio_union)-1));
            else:
                # Si el sitio rio arriba es negativo, se agrega como sitio reverse a L_sitios_cerca_complemento
                # L_sitios_rio_arriba[i] es la distancia del final del sitio a pos0
                # Se devuelve L_sitios_rio_arriba[i]-len(sitio_union)+1
                curr_sitio = int(L_sitios_rio_arriba[i])-len(sitio_union)+1;
                if forward:
                    ### PRUEBA:
                    #print('PRUEBA: Sitio de union buscado vs region encontrada. Sitio rio arriba negativo, gen forward.')
                    #print('Sitio de union buscado: ' + str(complemento_secuencia(sitio_union)))
                    #print('Sitio encontrado: ' + str(ConsultaSecuencia(chr_gen, pos0+curr_sitio, pos0+curr_sitio+len(sitio_union)-1, strand=1)))
                    L_sitios_cerca_complemento.append(curr_sitio);
                else:
                    ### PRUEBA:
                    #print('PRUEBA: Sitio de union buscado vs region encontrada. Sitio rio arriba negativo, gen reverse.')
                    #print('Sitio de union buscado: ' + str(complemento_secuencia(sitio_union)))
                    #print('Sitio encontrado: ' + str(ConsultaSecuencia(chr_gen, pos0-curr_sitio-len(sitio_union)+1, pos0-curr_sitio, strand=1)))
                    L_sitios_cerca_complemento.append(-(curr_sitio+len(sitio_union)-1));
            #print(curr_sitio)
        # Reviso los sitios de union rio abajo
        for j in range(len(L_sitios_rio_abajo)):
            # Si la secuencia rio abajo tiene un sitio de union, lo registro como numero positivo
            if L_sitios_rio_abajo[j] >= 0:
                # Si el sitio rio abajo es positivo, se agrega como sitio forward a L_sitios_cerca
                # L_sitios_rio_abajo[j] es la distancia del sitio a pos0
                # Se devuelve L_sitios_rio_abajo[j]
                curr_sitio = int(L_sitios_rio_abajo[j]);
                if forward:
                    ### PRUEBA:
                    #print('PRUEBA: Sitio de union buscado vs region encontrada. Sitio rio abajo positivo, gen forward.')
                    #print('Sitio de union buscado: ' + str(sitio_union))
                    #print('Sitio encontrado: ' + str(ConsultaSecuencia(chr_gen, pos0+curr_sitio, pos0+curr_sitio+len(sitio_union)-1, strand=1)))
                    L_sitios_cerca.append(curr_sitio);
                else:
                    ### PRUEBA:
                    #print('PRUEBA: Sitio de union buscado vs region encontrada. Sitio rio abajo positivo, gen reverse.')
                    #print('Sitio de union buscado: ' + str(sitio_union))
                    #print('Sitio encontrado: ' + str(ConsultaSecuencia(chr_gen, pos0-curr_sitio-len(sitio_union)+1, pos0-curr_sitio, strand=1)))
                    L_sitios_cerca.append(-(curr_sitio+len(sitio_union)-1));
            else:
                # Si el sitio rio abajo es negativo, se agrega como sitio reverse a L_sitios_cerca_complemento
                # L_sitios_rio_abajo[j] es la distancia del sitio a (A QUE) ################## HACER ESTO ##################
                #
                curr_sitio = rango_dist+L_sitios_rio_abajo[j];
                if forward:
                    ### PRUEBA:
                    #print('PRUEBA: Sitio de union buscado vs region encontrada. Sitio rio abajo negativo, gen forward.')
                    #print('Sitio de union buscado: ' + str(complemento_secuencia(sitio_union)))
                    #print('Sitio encontrado: ' + str(ConsultaSecuencia(chr_gen, pos0+curr_sitio-len(sitio_union)+1, pos0+curr_sitio, strand=1)))
                    L_sitios_cerca_complemento.append(curr_sitio-len(sitio_union)+1);
                else:
                    ### PRUEBA:
                    #print('PRUEBA: Sitio de union buscado vs region encontrada. Sitio rio abajo negativo, gen reverse.')
                    #print('Sitio de union buscado: ' + str(complemento_secuencia(sitio_union)))
                    #print('Sitio encontrado: ' + str(ConsultaSecuencia(chr_gen, pos0-curr_sitio, pos0-curr_sitio+len(sitio_union)-1, strand=1)))
                    L_sitios_cerca_complemento.append(-(curr_sitio));
            #print(curr_sitio)
        # Agrego ambas listas a la matriz de sitios de union que devuelve la funcion
        M_sitios_cerca.append(L_sitios_cerca[:]);
        M_sitios_cerca.append(L_sitios_cerca_complemento[:]);
    return M_sitios_cerca


def buscar_sitios_union_cerca_genes(L_genes, genome, genome_name, rango_dist, L_sitios_union, ext_out='.csv', sep_out=';'):
    # Busca, para cada uno de los genes, sitios de union en genoma a una distancia rango_dist del inicio del gen
    # La lista de genes es en formato gene_id

    # Inicializo una lista para devolver con genes para los que no se haya encontrado nada
    L_error = [];

    # Trabajo con cada uno de los genes, de a uno por vez
    for gen_id in L_genes:
        # Agarro el gen del genoma
        gen_busq = genome.gene_by_id(gen_id);
        # Creo el archivo donde guardo el output y elimino cualquier archivo con el mismo nombre
        with open(str(gen_busq.gene_name) + '_sitios_de_union_' + str(genome_name) + str(ext_out), 'w') as F_out:
            print('>Archivo ' + str(gen_busq.gene_name) + '_sitios_de_union_' + str(genome_name) + str(ext_out) + ' creado.')
        print(gen_busq)
        # Lo abro en modo append
        with open(str(gen_busq.gene_name) + '_sitios_de_union_' + str(genome_name) + str(ext_out), 'a') as F_out:
            # Defino el inicio del gen [esto lo uso al final pero se vuelve a calcular en buscar_sitios_cerca_gen()]
            if gen_busq.strand == '+':
                pos0 = gen_busq.start;
            elif gen_busq.strand == '-':
                pos0 = gen_busq.end;
            else:
                print('ERROR en strand de gen ' + str(gen_busq))
                pos0 = gen_busq.start;
            # Busco los sitios de union alrededor del inicio del gen
            L_sitios_cerca_gen = buscar_sitios_cerca_gen(gen_busq, genome, genome_name, rango_dist, L_sitios_union);
            # Guardo cada uno de los sitios encontrados
            for sitio_cerca in L_sitios_cerca_gen:
                #################### PRUEBA ####################
                #print('L_sitios: ' + str(sitio_cerca))
                #for q in sitio_cerca[1:]:
                    #_print_seq_sitio(gen_busq, genome_name, q, len(sitio_cerca[0]));
                ################## FIN PRUEBA ##################
                r = str(gen_busq.gene_id);
                for n in range(len(sitio_cerca)):
                    if n == 0:
                        r = r + str(sep_out) + str(sitio_cerca[n]);
                    else:
                        r = r + str(sep_out) + str(sitio_cerca[n]); # pos0 + sitio_cerca[n] para obtener la posicion en el genoma
                        ### PRUEBA
                        print(sitio_cerca[n])
                        #print(str(ConsultaSecuencia(IDchr(str(gen_busq.contig),genome=genome_name), pos0+int(sitio_cerca[n]), pos0+int(sitio_cerca[n])+len(sitio_cerca[0])-1, strand=1)))
                r = r + '\n';
                F_out.write(r);
        print()
    return L_error


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


def largo_archivo(dir_arch, nom_arch, ext, verbose=True):
    # Devuelve el largo (en filas) de un archivo en la direccion dir_arch
    # nom_arch y ext solo sirven para display
    with open(dir_arch, 'r') as F:
        if verbose:
            print('> Archivo ' + str(nom_arch) + str(ext) + ' abierto.')
        l_arch = len(F.read().split('\n'));
    return l_arch


def _print_seq_sitio(gene, genome_name, posicion_sitio, largo_sitio):
    # Funcion para testeo
    # Printea la secuencia en posicion_sitio respecto del +1 de gene

    print('Sitio en posicion ' + str(posicion_sitio) + ':')

    # Busco el ID del cromosoma
    chr_gen = IDchr(str(gene.contig),genome=genome_name);

    # Determino la direccion del gen y defino la posicion del 0
    if gene.strand == '+':
        pos0 = gene.start;
    elif gene.strand == '-':
        pos0 = gene.end;
    else:
        print('ERROR en strand de gen ' + str(gene))
        pos0 = gene.start;

    # Consulto la secuencia correspondiente al sitio de union y la printeo
    print(str(ConsultaSecuencia(chr_gen, pos0+posicion_sitio, pos0+posicion_sitio+largo_sitio-1, strand=1)))
    return None


def main():
    for n in range(len(lista_nombres)):
        print('>>>Iniciando parseo de "' + str(lista_nombres[n]) + str(extension) + '".')
        file_path = os.path.join(curr_path, str(path_archivos) + str(lista_nombres[n]) + str(extension));
        '''
        L_out1 = buscar_genes_picos_csv(file_path, lista_nombres[n], extension, csv_sep, rango_distancias, L_genomes[n], sep_out=';', print_range=500);
        print(L_out1)
        '''
        '''
        if n == 0:
            L_genes = L_human;
        else:
            L_genes = L_mouse;
        L_out2 = buscar_picos_cerca_de_genes(L_genes, L_genomes[n], rango_distancias, file_path, lista_nombres[n], extension, ';', sep_out=';');
        print(L_out2)
        '''
        print('Parseo de "' + str(lista_nombres[n]) + str(extension) + '" finalizado.')
        print()
    return ''


#################################### RESULTADOS ###################################


if __name__=='__main__':
    #main()
    L_out = buscar_sitios_union_cerca_genes(L_genes, genome, genome_name, rango_distancias, L_sitios_union);
    #print(L_out)
