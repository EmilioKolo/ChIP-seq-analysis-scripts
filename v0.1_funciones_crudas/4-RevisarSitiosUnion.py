
import os
import time
import copy
from pyensembl import EnsemblRelease
from Bio import Entrez, SeqIO
Entrez.email = 'ekolomenski@gmail.com';
Entrez.api_key = '9bccbd4949f226a5930868df36c211e8b408';
from referencias_funciones import ConsultaSecuencia, IDchr

mm9 = EnsemblRelease(67, species='mouse');
hg19 = EnsemblRelease(102, species='human');

GRCm38 = EnsemblRelease(102, species='mouse');

'''
Funciones para revisar sitios de union cerca de un gen

Toma ID del gen, genoma y posicion relativa del sitio de union
    Devuelve secuencia de la region correspondiente
    Variable opcional de distancia extra
'''

#################################### VARIABLES ####################################


# Chequear posiciones de los sitios de union
L_human = ['ENSG00000175206', 'ENSG00000265107', 'ENSG00000081189', 'ENSG00000171476', 
'ENSG00000141052', 'ENSG00000163485', 'ENSG00000179218', 'ENSG00000183023', 'ENSG00000117298', 
'ENSG00000164093', 'ENSG00000073146', 'ENSG00000141448'];
L_mouse = ['ENSMUSG00000041616', 'ENSMUSG00000057123', 'ENSMUSG00000005583', 'ENSMUSG00000059325', 
'ENSMUSG00000020542', 'ENSMUSG00000042429', 'ENSMUSG00000003814', 'ENSMUSG00000054640', 'ENSMUSG00000057530', 
'ENSMUSG00000028023', 'ENSMUSG00000015365', 'ENSMUSG00000005836'];

L_human_falta = ['ENSG00000265107', 'ENSG00000171476', 'ENSG00000163485'];
L_mouse_falta = ['ENSMUSG00000057123', 'ENSMUSG00000059325', 'ENSMUSG00000042429'];
rangos_falta = [(-1500, -1), (-1500, -1), (-500, -1)];

L_human_dificiles = ['ENSG00000164093', 'ENSG00000073146'];
L_mouse_dificiles = ['ENSMUSG00000028023', 'ENSMUSG00000015365'];

rangos_dificiles = [] ############### BUSCAR ###############

'''
Genes:
>Nppa/ANF (NPPA)
Human: ENSG00000175206; Mouse: ENSMUSG00000041616

>Connexin 40 (GJA5)
Human: ENSG00000265107; Mouse: ENSMUSG00000057123

>Mef2c (MEF2C)
Human: ENSG00000081189; Mouse: ENSMUSG00000005583

>Hop (HOPX)
Human: ENSG00000171476; Mouse: ENSMUSG00000059325

>Myocardin (MYOCD)
Human: ENSG00000141052; Mouse: ENSMUSG00000020542

>A1 adenosine receptor (ADORA1)
Human: ENSG00000163485; Mouse: ENSMUSG00000042429

>Calreticulin (CALR)
Human: ENSG00000179218; Mouse: ENSMUSG00000003814

>Sodium-calcium exchanger 1 (SLC8A1)
Human: ENSG00000183023; Mouse: ENSMUSG00000054640

>Endothelin-converting enzyme 1 (ECE1)
Human: ENSG00000117298; Mouse: ENSMUSG00000057530

>Pitx2 (PITX2)
Human: ENSG00000164093; Mouse: ENSMUSG00000028023

>Csm (MOV10L1) (isoforma)
Human: ENSG00000073146; Mouse: ENSMUSG00000015365

>GATA6 (GATA6)
Human: ENSG00000141448; Mouse: ENSMUSG00000005836
'''


L_genes = L_mouse;
genome = mm9;
genome_name = 'mm9';


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


def revisar_secuencia_gen(gen_id, genoma, nombre_genoma, d_ini, d_end):
    # Recibe un id de un gen, lo busca en un genoma y devuelve una secuencia
    # d_ini y d_end son distancias relativas al +1 del gen

    # Extraigo el gen correspondiente del genoma
    curr_gene = genoma.gene_by_id(str(gen_id));

    # Defino la posicion del +1 en base a la info que incluye el gen
    pos0, forward = buscar_pos0(curr_gene);
    
    # Defino el cromosoma de posicion del gen (para busqueda de secuencia
    chr_gen = IDchr(str(curr_gene.contig),genome=nombre_genoma);
    # Busco la secuencia entre pos0 +/- d_ini/d_end
    seq_buscada = secuencia_respecto_pos0(chr_gen, pos0, d_ini, d_end, forward);

    return seq_buscada


def secuencia_respecto_pos0(chr_gen, pos0, d_ini, d_end, forward):
    # Devuelve una secuencia respecto al +1 de un gen dado un id de cromosoma y posiciones relativas
    if forward:
        seq_buscada = str(ConsultaSecuencia(chr_gen, pos0+d_ini, pos0+d_end, strand=1));
    else:
        seq_buscada = str(ConsultaSecuencia(chr_gen, pos0-d_end, pos0-d_ini, strand=1));
    return seq_buscada


def secuencias_sitios(gen_id, genoma, nombre_genoma, L_sitios, largo_sitios):
    # Devuelve lista con secuencias en sitios ubicados a distancia L_sitios del +1 del gen dado

    # Lista de secuencias encontradas para que devuelva la funcion
    L_secuencias = [];

    # Extraigo el gen correspondiente del genoma
    curr_gene = genoma.gene_by_id(str(gen_id));
    print(curr_gene)

    # Defino la posicion del +1 en base a la info que incluye el gen
    pos0, forward = buscar_pos0(curr_gene);

    # Defino el cromosoma de posicion del gen (para busqueda de secuencia
    chr_gen = IDchr(str(curr_gene.contig),genome=nombre_genoma);

    # Corro una funcion por cada sitio
    for sitio in L_sitios:
        print('Sitio: ' + str(sitio))
        r = secuencia_un_sitio(chr_gen, pos0, int(sitio), largo_sitios, forward);
        print(r)
        L_secuencias.append(str(r));
    return L_secuencias


def secuencia_un_sitio(chr_gen, pos0, sitio, largo_sitio, forward):
    # Devuelve la secuencia de un sitio ubicado a distancia sitio de pos0 en chr_gen

    if forward:
        seq_sitio = str(ConsultaSecuencia(chr_gen, pos0+sitio, pos0+sitio+largo_sitio-1, strand=1));
    else:
        seq_sitio = str(ConsultaSecuencia(chr_gen, pos0-sitio, pos0-sitio+largo_sitio, strand=1));
    return seq_sitio


def main1():
    # Primera revision a ojo de secuencias de promotores
    print('ANF')
    print(genome.gene_by_id(str(L_mouse[0])))
    print(revisar_secuencia_gen(L_mouse[0], genome, genome_name, -96, -76))
    print(revisar_secuencia_gen(L_mouse[0], genome, genome_name, -250, -220))
    print()
    print('Connexin 40') # Buscar VARIOS sitios
    print(genome.gene_by_id(str(L_mouse[1])))
    print(revisar_secuencia_gen(L_mouse[1], genome, genome_name, -970, -1030))
    print(revisar_secuencia_gen(L_mouse[1], genome, genome_name, -740, -800))
    print(revisar_secuencia_gen(L_mouse[1], genome, genome_name, -470, -530))
    print(revisar_secuencia_gen(L_mouse[1], genome, genome_name, -130, -190))
    print(revisar_secuencia_gen(L_mouse[1], genome, genome_name, -100, -40))
    print()
    print('Mef2c')
    print(genome.gene_by_id(str(L_mouse[2])))
    print(revisar_secuencia_gen(L_mouse[2], genome, genome_name, 21748, 21760))
    print()
    print('Hop') # Buscar sitio
    print(genome.gene_by_id(str(L_mouse[3])))
    print(revisar_secuencia_gen(L_mouse[3], genome, genome_name, -270, -200))
    print()
    print('Myocardin')
    print(genome.gene_by_id(str(L_mouse[4])))
    print(revisar_secuencia_gen(L_mouse[4], genome, genome_name, -1070, -1060))
    print()
    print('A1 adenosine receptor')
    print(genome.gene_by_id(str(L_mouse[5])))
    print(revisar_secuencia_gen(L_mouse[5], genome, genome_name, -245, -233))
    print()
    print('Calreticulin') # Buscar sitio
    print(genome.gene_by_id(str(L_mouse[6])))
    print(revisar_secuencia_gen(L_mouse[6], genome, genome_name, -1266, -1200))
    print()
    print('Sodium-calcium exchanger 1')
    print(genome.gene_by_id(str(L_mouse[7])))
    print(revisar_secuencia_gen(L_mouse[7], genome, genome_name, -12, -1))
    print()
    print('Endothelin-converting enzyme 1')
    print(genome.gene_by_id(str(L_mouse[8])))
    print(revisar_secuencia_gen(L_mouse[8], genome, genome_name, -1196, -1185))
    print()
    print('Pitx2')
    print(genome.gene_by_id(str(L_mouse[9])))
    print('HAY QUE BUSCAR LAS PRIMERAS 100 PB DEL ASE')
    #print(revisar_secuencia_gen(L_mouse[9], genome, genome_name, ???, ???))
    print()
    print('Csm')
    print(genome.gene_by_id(str(L_mouse[10])))
    print(revisar_secuencia_gen(L_mouse[10], genome, genome_name, -66, -56))
    print()
    print('GATA6')
    print(genome.gene_by_id(str(L_mouse[11])))
    print(revisar_secuencia_gen(L_mouse[11], genome, genome_name, -3037, -3026))
    print()
    return None

def main2():
    # Busqueda de AAGTG en promotor de 'connexin 40', 'Hop' y 'A1 adenosine receptor'
    for i in range(len(L_human_falta)):
        print(mm9.gene_by_id(str(L_mouse_falta[i])))
        print(hg19.gene_by_id(str(L_human_falta[i])))
        rango_i = rangos_falta[i];
        seq_promotor_mouse = revisar_secuencia_gen(L_mouse_falta[i], mm9, 'mm9', rango_i[0], rango_i[1]);
        seq_promotor_human = revisar_secuencia_gen(L_human_falta[i], hg19, 'hg19', rango_i[0], rango_i[1]);
        if mm9.gene_by_id(str(L_mouse_falta[i])).strand == '-':
            seq_promotor_mouse = complemento_secuencia(seq_promotor_mouse);
        if hg19.gene_by_id(str(L_human_falta[i])).strand == '-':
            seq_promotor_human = complemento_secuencia(seq_promotor_human);
        print('Posiciones de la secuencia AAGTG en humano:')
        print(buscar_en_seq_2dir_unificado('AAGTG', seq_promotor_human))
        print('Posiciones de la secuencia AAGTG en raton:')
        print(buscar_en_seq_2dir_unificado('AAGTG', seq_promotor_mouse))
        print()
    return None

def main2_test():
    # Prueba para ver los sitios de union sacados de main2()
    print('Connexin 40 / GJA5')
    L1 = secuencias_sitios(L_mouse_falta[0], mm9, 'mm9', [-1310, -1039, -735, -374, -126, -66, -272, -868], 5);
    L2 = secuencias_sitios(L_human_falta[0], hg19, 'hg19', [-1271, -1101, -47, -627, -873, -1128], 5);
    #print('Mouse:')
    #for j in L1:
        #print(j)
    #print('Human:')
    #for j in L2:
        #print(j)
    print()

    print('Hop / HOPX')
    L3 = secuencias_sitios(L_mouse_falta[1], mm9, 'mm9', [-1314, -1080, -835, -804, -646], 5);
    L4 = secuencias_sitios(L_human_falta[1], hg19, 'hg19', [-1490, -450, -425, -387, -1050, -1146], 5);
    #print('Mouse:')
    #for j in L3:
        #print(j)
    #print('Human:')
    #for j in L4:
        #print(j)
    print()
    
    print('A1 adenosine receptor / ADORA1')
    L5 = secuencias_sitios(L_mouse_falta[2], mm9, 'mm9', [-120], 5);
    L6 = secuencias_sitios(L_human_falta[2], hg19, 'hg19', [-26], 5);
    #print('Mouse:')
    #for j in L5:
        #print(j)
    #print('Human:')
    #for j in L6:
        #print(j)
    print()
    return None

def main_seq_busq():
    # Prueba de buscar_en_seq_2dir_unificado(busq, seq)
    seq_test = 'CCAACGTATGTT'
    L = buscar_en_seq_2dir_unificado('AA', seq_test)
    for i in L:
        print(seq_test[i:i+2])
    return None

#################################### RESULTADOS ###################################


if __name__=='__main__':
    #main1();
    #main2();
    #main2_test();
    #main_seq_busq();
    pass

'''
>>> RESULTADOS main2() para 'connexin 40', 'Hop' y 'A1 adenosine receptor'

[(-1500, -1), (-1500, -1), (-500, -1)]

> Connexin 40 / GJA5
Posiciones de la secuencia AAGTG en humano:
[-1271, -1132, -1101, -877, -631, -47]
Posiciones de la secuencia AAGTG en raton:
[-1310, -1039, -872, -735, -374, -276, -126, -70]

for i in [-1310, -1039, -735, -374, -126, -70, -276, -872]:
    str(ConsultaSecuencia(chr_3_mouse, 96708616+i-2, 96708616+i+4+2, strand=1))

'GCAAGTGCC'	GCAAGTG		-1310
'GAAAGTGGC'	GAAAGTG		-1039
'GAAAGTGCA'	GAAAGTG		-735
'TTAAGTGTT'	TTAAGTG		-374
'GAAAGTGGG'	GAAAGTG		-126
'AACACTTGG'	CACTTGG		-70
'CTCACTTAC'	CACTTAC		-276
'AGCACTTTG'	CACTTTG		-872

for i in [-1271, -1101, -47, -631, -877, -1132]:
    str(ConsultaSecuencia(chr_1_human, 147773362-i-4-2, 147773362-i+2, strand=1))

'GCCACTTGT'	CACTTGT		-1271
'AACACTTAG'	CACTTAG		-1101
'CCCACTTCT'	CACTTCT		-47
'AAAAGTGTA'	AAAAGTG		-631
'GTAAGTGGT'	GTAAGTG		-877
'GTAAGTGTC'	GTAAGTG		-1132


> Hop / HOPX
Posiciones de la secuencia AAGTG en humano:
[-1490, -1150, -1054, -450, -425, -391]
Posiciones de la secuencia AAGTG en raton:
[-1314, -1080, -835, -804, -646]

for i in [-1314, -1080, -835, -804, -646]:
    str(ConsultaSecuencia(chr_5_mouse, 77544186-i-4-2, 77544186-i+2, strand=1))

'TTCACTTTC'	CACTTTC		-1314
'ATCACTTCA'	CACTTCA		-1080
'GACACTTGC'	CACTTGC		-835
'TTCACTTCT'	CACTTCT		-804
'CCCACTTGA'	CACTTGA		-646

for i in [-1490, -450, -425, -391, -1054, -1150]:
    str(ConsultaSecuencia(chr_4_human, 56681899-i-4-2, 56681899-i+2, strand=1))

'GGCACTTGG'	CACTTGG		-1490
'CCCACTTGA'	CACTTGA		-450
'TTCACTTCA'	CACTTCA		-425
'GTAAGTGCC'	GTAAGTG		-391
'AGAAGTGGC'	AGAAGTG		-1054
'CAAAGTGCT'	CAAAGTG		-1150


> A1 adenosine receptor / ADORA1
Posiciones de la secuencia AAGTG en humano:
[-26]
Posiciones de la secuencia AAGTG en raton:
[-125]

str(ConsultaSecuencia(chr_1_mouse, 136132008+125-4-2, 136132008+125+2, strand=1))

'GGAAGTGAT'	GGAAGTG		-125

str(ConsultaSecuencia(chr_1_human, 203090654-26-2, 203090654-26+4+2, strand=1))

'GGAAGTGAC'	GGAAGTG		-26
'''

##### PRUEBAS

gen_gja5_human = hg19.gene_by_id('ENSG00000265107');
chr_1_human = IDchr('1',genome='hg19');
gen_gja5_mouse = mm9.gene_by_id('ENSMUSG00000057123');
chr_3_mouse = IDchr('3',genome='mm9');
gen_gja5_mouse_modern = GRCm38.gene_by_id('ENSMUSG00000057123');
chr_3_mouse_modern = IDchr('3',genome='mouse102');

gen_hopx_human = hg19.gene_by_id('ENSG00000171476');
chr_4_human = IDchr('4',genome='hg19');
gen_hopx_mouse = mm9.gene_by_id('ENSMUSG00000059325');
chr_5_mouse = IDchr('5',genome='mm9');
gen_hopx_mouse_modern = GRCm38.gene_by_id('ENSMUSG00000059325');
chr_5_mouse_modern = IDchr('5',genome='mouse102');

gen_adora1_human = hg19.gene_by_id('ENSG00000163485');
chr_1_human = IDchr('1',genome='hg19');
gen_adora1_mouse = mm9.gene_by_id('ENSMUSG00000042429');
chr_1_mouse = IDchr('1',genome='mm9');
gen_adora1_mouse_modern = GRCm38.gene_by_id('ENSMUSG00000042429');
chr_1_mouse_modern = IDchr('1',genome='mouse102');

gen_pitx2_human = hg19.gene_by_id('ENSG00000164093');
chr_4_human = IDchr('4',genome='hg19');
gen_pitx2_mouse = mm9.gene_by_id('ENSMUSG00000028023');
chr_3_mouse = IDchr('3',genome='mm9');
gen_pitx2_mouse_modern = GRCm38.gene_by_id('ENSMUSG00000028023');
chr_3_mouse_modern = IDchr('3',genome='mouse102');

gen_mov10l1_human = hg19.gene_by_id('ENSG00000073146');
chr_22_human = IDchr('22',genome='hg19');
gen_mov10l1_mouse = mm9.gene_by_id('ENSMUSG00000015365');
chr_15_mouse = IDchr('15',genome='mm9');
gen_mov10l1_mouse_modern = GRCm38.gene_by_id('ENSMUSG00000015365');
chr_15_mouse_modern = IDchr('15',genome='mouse102');

exon_ids_mov10l1_mouse_modern = GRCm38.exon_ids_of_gene_id('ENSMUSG00000015365');
exon_ids_mov10l1_human = hg19.exon_ids_of_gene_id('ENSG00000073146');

exon_15_mov10l1_mouse_modern = GRCm38.exon_by_id('ENSMUSE00000128734');
exon_16_mov10l1_mouse_modern = GRCm38.exon_by_id('ENSMUSE00000311628');
exon_17_mov10l1_mouse_modern = GRCm38.exon_by_id('ENSMUSE00000379757');

'''
Exones importantes MOV10L1 ('ENSMUSG00000015365')
Raton:
exon 15: 89011900 - 89012195 ('ENSMUSE00000128734') # sin AAGTG
exon 16: 89018150 - 89018249 ('ENSMUSE00000311628') # AAGTG en .end -14 (GCAAGTG en -16)
exon 17: 89020261 - 89020439 ('ENSMUSE00000379757') # AAGTG en .start -8 (CTAAGTG en -10)

GRCm38.exon_by_id(exon_id)
'''

# ConsultaSecuencia(chr_gen, pos_ini, pos_end, strand=1)


############ Genes:

### Nppa/ANF (NPPA)

# NKE1: -94 / -78
# NKE2: -247 / ~ -220

#Human: ENSG00000175206
#Gene(gene_id='ENSG00000175206', gene_name='NPPA', biotype='protein_coding', contig='1', start=11845709, end=11848345, strand='-', genome='GRCh38')
gen_nppa_human = hg19.gene_by_id('ENSG00000175206');
chr_1_human = IDchr('1',genome='hg19');
#buscar_en_seq_2dir_unificado('AAGTG',complemento_secuencia(str(ConsultaSecuencia(chr_1_human, gen_nppa_human.end, gen_nppa_human.end+500, strand=1))))
# [-406, -413]
str(ConsultaSecuencia(chr_1_human, gen_nppa_human.end+399, gen_nppa_human.end+412, strand=1)) # -412 a -399; 'GCAAGTGGCAAGTG'

#Mouse: ENSMUSG00000041616
#Gene(gene_id='ENSMUSG00000041616', gene_name='Nppa', biotype='protein_coding', contig='4', start=147374831, end=147376188, strand='+', genome='NCBIM37')
gen_nppa_mouse = mm9.gene_by_id('ENSMUSG00000041616');
chr_4_mouse = IDchr('4',genome='mm9');
#buscar_en_seq_2dir_unificado('AAGTG',str(ConsultaSecuencia(chr_4_mouse, gen_nppa_mouse.start-500, gen_nppa_mouse.start, strand=1)))
# [-240, -66]
str(ConsultaSecuencia(chr_4_mouse, gen_nppa_mouse.start-67, gen_nppa_mouse.start-52, strand=1))   # -67 a -52; 'GCAAGTGACAGAATGG'
str(ConsultaSecuencia(chr_4_mouse, gen_nppa_mouse.start-241, gen_nppa_mouse.start-235, strand=1)) # -241 a -235; 'TGAAGTG'



### Myocardin (MYOCD)

# NKE: Empieza en -1069

#Human: ENSG00000141052
#Gene(gene_id='ENSG00000141052', gene_name='MYOCD', biotype='protein_coding', contig='17', start=12665890, end=12768949, strand='+', genome='GRCh38')
gen_myocd_human = hg19.gene_by_id('ENSG00000141052');
chr_17_human = IDchr('17',genome='hg19');
#buscar_en_seq_2dir_unificado('AAGTG', str(ConsultaSecuencia(chr_17_human, gen_myocd_human.start-5000, gen_myocd_human.start, strand=1)))
# [-4826, -4799, -2631, -1541, -2002, -3342, -4364]
str(ConsultaSecuencia(chr_17_human, gen_myocd_human.start-1540, gen_myocd_human.start-1534, strand=1)) # -1540 a -1534; 'CACTTTG'

#Mouse: ENSMUSG00000020542
#Gene(gene_id='ENSMUSG00000020542', gene_name='Myocd', biotype='protein_coding', contig='11', start=64990063, end=65083491, strand='-', genome='NCBIM37')
gen_myocd_mouse = mm9.gene_by_id('ENSMUSG00000020542');
chr_11_mouse = IDchr('11',genome='mm9');
#buscar_en_seq_2dir_unificado('AAGTG',complemento_secuencia(str(ConsultaSecuencia(chr_11_mouse, gen_myocd_mouse.end, gen_myocd_mouse.end+1500, strand=1))))
# [-1086, -491, -652]
str(ConsultaSecuencia(chr_11_mouse, gen_myocd_mouse.end+1081, gen_myocd_mouse.end+1087, strand=1)) # -1087 a -1081; 'CACTTGA'
str(ConsultaSecuencia(chr_11_mouse, gen_myocd_mouse.end+484, gen_myocd_mouse.end+490, strand=1))   # -490 a -484; 'AGAAGTG'
str(ConsultaSecuencia(chr_11_mouse, gen_myocd_mouse.end+645, gen_myocd_mouse.end+651, strand=1))   # -651 a -645; 'AAAAGTG'



### Calreticulin (CALR)

# CRT2: Entre -1266 y -902 (cerca de -1266)

#Human: ENSG00000179218
#Gene(gene_id='ENSG00000179218', gene_name='CALR', biotype='protein_coding', contig='19', start=12938578, end=12944489, strand='+', genome='GRCh38')
gen_calr_human = hg19.gene_by_id('ENSG00000179218');
chr_19_human = IDchr('19',genome='hg19');
#buscar_en_seq_2dir_unificado('AAGTG',str(ConsultaSecuencia(chr_19_human, gen_calr_human.start-1500, gen_calr_human.start, strand=1)))
# [-1372, -643, -55, -658, -1276]
str(ConsultaSecuencia(chr_19_human, gen_calr_human.start-56, gen_calr_human.start-50, strand=1))     # -56 a -50; 'CAAAGTG'
str(ConsultaSecuencia(chr_19_human, gen_calr_human.start-644, gen_calr_human.start-638, strand=1))   # -644 a -638; 'TTAAGTG'
str(ConsultaSecuencia(chr_19_human, gen_calr_human.start-657, gen_calr_human.start-651, strand=1))   # -657 a -651; 'CACTTCA'
str(ConsultaSecuencia(chr_19_human, gen_calr_human.start-1275, gen_calr_human.start-1269, strand=1)) # -1275 a -1269; 'CACTTTA'
str(ConsultaSecuencia(chr_19_human, gen_calr_human.start-1373, gen_calr_human.start-1367, strand=1)) # -1373 a -1367; 'CCAAGTG'

#Mouse: ENSMUSG00000003814
#Gene(gene_id='ENSMUSG00000003814', gene_name='Calr', biotype='protein_coding', contig='8', start=87365749, end=87370833, strand='-', genome='NCBIM37')
gen_calr_mouse = mm9.gene_by_id('ENSMUSG00000003814');
chr_8_mouse = IDchr('8',genome='mm9');
#buscar_en_seq_2dir_unificado('AAGTG',complemento_secuencia(str(ConsultaSecuencia(chr_8_mouse, gen_calr_mouse.end, gen_calr_mouse.end+1500, strand=1))))
# [-1166, -831]
str(ConsultaSecuencia(chr_8_mouse, gen_calr_mouse.end+826, gen_calr_mouse.end+832, strand=1))   # -832 a -826; 'CACTTCT'
str(ConsultaSecuencia(chr_8_mouse, gen_calr_mouse.end+1161, gen_calr_mouse.end+1167, strand=1)) # -1167 a -1161; 'CACTTGA'



### Sodium-calcium exchanger 1 (SLC8A1)

# NKE-10: -10 / -3

#Human: ENSG00000183023
#Gene(gene_id='ENSG00000183023', gene_name='SLC8A1', biotype='protein_coding', contig='2', start=40097270, end=40611053, strand='-', genome='GRCh38')
gen_slc8a1_human = hg19.gene_by_id('ENSG00000183023');
chr_2_human = IDchr('2',genome='hg19');
#buscar_en_seq_2dir_unificado('AAGTG',complemento_secuencia(str(ConsultaSecuencia(chr_2_human, gen_slc8a1_human.end, gen_slc8a1_human.end+500, strand=1))))
# [-71, -284, -333]
str(ConsultaSecuencia(chr_2_human, gen_slc8a1_human.end+66, gen_slc8a1_human.end+72, strand=1))   # -72 a -66; 'CACTTAT'
str(ConsultaSecuencia(chr_2_human, gen_slc8a1_human.end+277, gen_slc8a1_human.end+283, strand=1)) # -283 a -277; 'GTAAGTG'
str(ConsultaSecuencia(chr_2_human, gen_slc8a1_human.end+326, gen_slc8a1_human.end+332, strand=1)) # -332 a -326; 'CAAAGTG'

#Mouse: ENSMUSG00000054640
#Gene(gene_id='ENSMUSG00000054640', gene_name='Slc8a1', biotype='protein_coding', contig='17', start=81785359, end=82048947, strand='-', genome='NCBIM37')
gen_slc8a1_mouse = mm9.gene_by_id('ENSMUSG00000054640');
chr_17_mouse = IDchr('17',genome='mm9');
#buscar_en_seq_2dir_unificado('AAGTG',complemento_secuencia(str(ConsultaSecuencia(chr_17_mouse, gen_slc8a1_mouse.end, gen_slc8a1_mouse.end+500, strand=1))))
# [-121, -397]
str(ConsultaSecuencia(chr_17_mouse, gen_slc8a1_mouse.end+116, gen_slc8a1_mouse.end+122, strand=1)) # -122 a -116; 'CACTTAT'
str(ConsultaSecuencia(chr_17_mouse, gen_slc8a1_mouse.end+390, gen_slc8a1_mouse.end+396, strand=1)) # -396 a -390; 'ACAAGTG'



### Endothelin-converting enzyme 1 (ECE1)

# pECE-1b: Empieza en -1197

#Human: ENSG00000117298
#Gene(gene_id='ENSG00000117298', gene_name='ECE1', biotype='protein_coding', contig='1', start=21217247, end=21345572, strand='-', genome='GRCh38')
gen_ece1_human = hg19.gene_by_id('ENSG00000117298');
chr_1_human = IDchr('1',genome='hg19');
#buscar_en_seq_2dir_unificado('AAGTG',complemento_secuencia(str(ConsultaSecuencia(chr_1_human, gen_ece1_human.end, gen_ece1_human.end+1500, strand=1))))
# [-1255, -676]
str(ConsultaSecuencia(chr_1_human, gen_ece1_human.end+669, gen_ece1_human.end+675, strand=1))   # -675 a -669; 'GAAAGTG'
str(ConsultaSecuencia(chr_1_human, gen_ece1_human.end+1250, gen_ece1_human.end+1256, strand=1)) # -1256 a -1250; 'CACTTTT'

#Mouse: ENSMUSG00000057530
#Gene(gene_id='ENSMUSG00000057530', gene_name='Ece1', biotype='protein_coding', contig='4', start=137418152, end=137521144, strand='+', genome='NCBIM37')
gen_ece1_mouse = mm9.gene_by_id('ENSMUSG00000057530');
chr_4_mouse = IDchr('4',genome='mm9');
#buscar_en_seq_2dir_unificado('AAGTG',str(ConsultaSecuencia(chr_4_mouse, gen_ece1_mouse.start-1500, gen_ece1_mouse.start, strand=1)))
# [-601, -33]
str(ConsultaSecuencia(chr_4_mouse, gen_ece1_mouse.start-32, gen_ece1_mouse.start-26, strand=1))   # -32 a -26; 'CACTTCT'
str(ConsultaSecuencia(chr_4_mouse, gen_ece1_mouse.start-602, gen_ece1_mouse.start-596, strand=1)) # -602 a -596; 'CTAAGTG'



### GATA6 (GATA6)

# NKE: Empieza en -3035

#Human: ENSG00000141448
#Gene(gene_id='ENSG00000141448', gene_name='GATA6', biotype='protein_coding', contig='18', start=22169589, end=22202528, strand='+', genome='GRCh38')
gen_gata6_human = hg19.gene_by_id('ENSG00000141448');
chr_18_human = IDchr('18',genome='hg19');
#buscar_en_seq_2dir_unificado('AAGTG',str(ConsultaSecuencia(chr_18_human, gen_gata6_human.start-4000, gen_gata6_human.start, strand=1)))
# [-3799, -1755, -1739, -1679, -2012, -2536]
str(ConsultaSecuencia(chr_18_human, gen_gata6_human.start-1680, gen_gata6_human.start-1674, strand=1)) # -1680 a -1674; 'AGAAGTG'
str(ConsultaSecuencia(chr_18_human, gen_gata6_human.start-1756, gen_gata6_human.start-1734, strand=1)) # -1756 a -1734; 'CAAAGTGGAACGTAATAAAAGTG'
#    str(ConsultaSecuencia(chr_18_human, gen_gata6_human.start-1740, gen_gata6_human.start-1734, strand=1)) # -1740 a -1734; 'AAAAGTG'
#    str(ConsultaSecuencia(chr_18_human, gen_gata6_human.start-1756, gen_gata6_human.start-1750, strand=1)) # -1756 a -1750; 'CAAAGTG'
str(ConsultaSecuencia(chr_18_human, gen_gata6_human.start-2011, gen_gata6_human.start-2005, strand=1)) # -2011 a -2005; 'CACTTAA'
str(ConsultaSecuencia(chr_18_human, gen_gata6_human.start-2535, gen_gata6_human.start-2529, strand=1)) # -2535 a -2529; 'CACTTAT'
str(ConsultaSecuencia(chr_18_human, gen_gata6_human.start-3800, gen_gata6_human.start-3794, strand=1)) # -3800 a -3794; 'CGAAGTG'

#Mouse: ENSMUSG00000005836
#Gene(gene_id='ENSMUSG00000005836', gene_name='Gata6', biotype='protein_coding', contig='18', start=11052508, end=11085633, strand='+', genome='NCBIM37')
gen_gata6_mouse = mm9.gene_by_id('ENSMUSG00000005836');
chr_18_mouse = IDchr('18',genome='mm9');
#buscar_en_seq_2dir_unificado('AAGTG',str(ConsultaSecuencia(chr_18_mouse, gen_gata6_mouse.start-5000, gen_gata6_mouse.start, strand=1)))
# [-3053, -1661, -1600, -3278]
str(ConsultaSecuencia(chr_18_mouse, gen_gata6_mouse.start-1601, gen_gata6_mouse.start-1595, strand=1)) # -1601 a -1595; 'AGAAGTG'
str(ConsultaSecuencia(chr_18_mouse, gen_gata6_mouse.start-1662, gen_gata6_mouse.start-1656, strand=1)) # -1662 a -1656; 'AAAAGTG'
str(ConsultaSecuencia(chr_18_mouse, gen_gata6_mouse.start-3054, gen_gata6_mouse.start-3048, strand=1)) # -3054 a -3048; 'GCAAGTG'
str(ConsultaSecuencia(chr_18_mouse, gen_gata6_mouse.start-3277, gen_gata6_mouse.start-3271, strand=1)) # -3277 a -3271; 'CACTTTT'

'''
### Mef2c (MEF2C)

# N: ~ +21750 / +21758

#Human: ENSG00000081189
#Gene(gene_id='ENSG00000081189', gene_name='MEF2C', biotype='protein_coding', contig='5', start=88717117, end=88904257, strand='-', genome='GRCh38')
gen_mef2c_human = hg19.gene_by_id('ENSG00000081189');
chr_5_human = IDchr('5',genome='hg19');

#buscar_en_seq_2dir_unificado('AAGTG',complemento_secuencia(str(ConsultaSecuencia(chr_5_human, gen_mef2c_human.end, gen_mef2c_human.end+5000, strand=1))))
# [-4755, -3052, -2381, -2142, -1856, -962, -775, -1083, -1524, -2344, -2365, -2406, -2433, -3954, -4545, -4968]

#buscar_en_seq_2dir_unificado('AAGTG',complemento_secuencia(str(ConsultaSecuencia(chr_5_human, gen_mef2c_human.end-22000, gen_mef2c_human.end-20000, strand=1))))
# [-1672, -157, -811, -917, -1110, -1178] (restar 20000 al usarlos en gen.end)

#buscar_en_seq_2dir_unificado('AAGTG',complemento_secuencia(str(ConsultaSecuencia(chr_5_human, gen_mef2c_human.end-22000, gen_mef2c_human.end-21500, strand=1))))
# [-157] (restar 21500 al usarlos en gen.end)


#Mouse: ENSMUSG00000005583
#Gene(gene_id='ENSMUSG00000005583', gene_name='Mef2c', biotype='protein_coding', contig='13', start=83643033, end=83806684, strand='+', genome='NCBIM37')
gen_mef2c_mouse = mm9.gene_by_id('ENSMUSG00000005583');
chr_13_mouse = IDchr('13',genome='mm9');

#buscar_en_seq_2dir_unificado('AAGTG',str(ConsultaSecuencia(chr_13_mouse, gen_mef2c_mouse.start-5000, gen_mef2c_mouse.start, strand=1)))
# [-3376, -2745, -2184, -2109, -1939, -1304, -942, -734, -432, -234, -757, -1062, -1340, -1478, -1503, -3077, -4379, -4487]

#buscar_en_seq_2dir_unificado('AAGTG',str(ConsultaSecuencia(chr_13_mouse, gen_mef2c_mouse.start+20000, gen_mef2c_mouse.start+22000, strand=1)))
# [-1939, -1018, -358, -152, -78, -199, -1007, -1113, -1157, -1392, -1468] (sumar 22000 al usarlos en gen.start)

#buscar_en_seq_2dir_unificado('AAGTG',str(ConsultaSecuencia(chr_13_mouse, gen_mef2c_mouse.start+21500, gen_mef2c_mouse.start+22000, strand=1)))
# [-358, -152, -78, -199] (sumar 22000 al usarlos en gen.start)

'''
