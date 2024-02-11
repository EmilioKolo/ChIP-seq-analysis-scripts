
ib3 = False; 
http_proxy_ib3 = 'http://proxy.fcen.uba.ar:8080'; 
https_proxy_ib3 = 'https://proxy.fcen.uba.ar:8080'; 
ftp_proxy_ib3 = 'ftp://proxy.fcen.uba.ar:8080'; 

# Generales
import os
import time
import copy
from random import shuffle
from pyensembl import EnsemblRelease

# Analisis de secuencias y genomas
from Bio import Entrez, SeqIO, motifs, Seq
Entrez.email = 'ekolomenski@gmail.com'; 
Entrez.api_key = '9bccbd4949f226a5930868df36c211e8b408'; 
import biomart

# Importo analisischip
import sys
path_analisischip_casa = 'C:\\Users\\Admin\\Documents\\Scripts\\Git\\AnalisisSecuenciasChIP\\'; 
path_analisischip_ib3 = 'C:\\Users\\emili\\Documents\\Archivos nube\\Git\\analisissecuenciaschip\\'; 
if ib3:
    path_analisischip_main = path_analisischip_ib3; 
else:
    path_analisischip_main = path_analisischip_casa; 
sys.path.append(path_analisischip_main); 
from analisischip.seq import seq_data, seq_handler

# Genomas de referencia
hg19 = EnsemblRelease(75, species='human'); 
hg38 = EnsemblRelease(102, species='human'); 
mm9 = EnsemblRelease(54, species='mouse'); 
mm10 = EnsemblRelease(102, species='mouse'); 


'''
Funcion que busca sitios de union alrededor del +1 de un gen

gen_id_gen_name() busca el id de un nombre de gen en un genoma dado

seq_alrededor_ini_gen() busca la secuencia alrededor de un id de gen

buscar_sitios_gen() busca con lista de sitios alrededor de un gen por id

buscar_L_pssm_gen() busca con lista de matrices pssm alrededor de un gen por id

Las matrices generadas por buscar_sitios_gen() y buscar_L_pssm_gen() se pueden sumar

guardar_hits() guarda la matriz generada por buscar_sitios_gen() y buscar_L_pssm_gen() 
'''

#################################### VARIABLES ####################################

# Direcciones de archivos y outputs para ib3 y casa
path_dropbox_casa = 'D:\\Users\\Admin\\Dropbox\\'; 
path_dropbox_ib3 = 'C:\\Users\\emili\\Dropbox\\'; 
path_fasta_casa = 'D:\\Archivos doctorado\\Genomas\\'; 
path_fasta_ib3 = 'X:\\Genomas\\'; 
path_output_dump_casa = 'D:\\Archivos doctorado\\Output_dump\\'; 
path_output_dump_ib3 = 'X:\\Output_dump\\'; 

# Nombres input
nkx25_mouse = 'Nkx2-5'; 
nkx25_human = 'NKX2-5'; 

L_sitios = ['GCAAGTG', 'GGAAGTG', 'GAAAGTG', 'ATAAGTG', 'GTAAGTG', 'CTAAGTG', 'TCAAGTG', 'TGAAGTG', 'TAAAGTG', 'TTAAGTG']; 

# Nombres de archivos de matrices de peso y factores de transcripcion asociados
L_arch_pwm_raton = ['NKX25_MOUSE.H11MO.0.A.pcm', 'TBX20_MOUSE.H11MO.0.C.pcm', 'MEIS1_MOUSE.H11MO.1.A.pcm', 'TGIF1_MOUSE.H11MO.0.A.pcm', 
                    'HAND1_MOUSE.H11MO.0.C.pcm', 'MAF_MOUSE.H11MO.1.A.pcm', 'GATA1_MOUSE.H11MO.1.A.pcm', 'GATA6_MOUSE.H11MO.0.A.pcm', 'GATA4_MOUSE.H11MO.0.A.pcm']; 
L_arch_pwm_humano = ['NKX25_HUMAN.H11MO.0.B.pcm', 'TBX20_HUMAN.H11MO.0.D.pcm', 'MEIS1_HUMAN.H11MO.1.B.pcm', 'TGIF1_HUMAN.H11MO.0.A.pcm', 'HAND1_HUMAN.H11MO.0.D.pcm', 
                     'MAF_HUMAN.H11MO.1.B.pcm', 'GATA1_HUMAN.H11MO.1.A.pcm', 'GATA6_HUMAN.H11MO.0.A.pcm', 'GATA4_HUMAN.H11MO.0.A.pcm']; 
L_nombres_pssm_raton = ['Nkx25', 'Tbx20', 'Meis1', 'Tgif1', 'Hand1', 'Maf', 'Gata1', 'Gata6', 'Gata4']; 
L_nombres_pssm_humano = ['NKX25', 'TBX20', 'MEIS1', 'TGIF1', 'HAND1', 'MAF', 'GATA1', 'GATA6', 'GATA4']; 


### Variables main()

# Path main depende de si estoy en ib3 o casa
if ib3:
    path_dropbox_main = path_dropbox_ib3; 
    path_fasta_main = path_fasta_ib3; 
    path_output_dump_main = path_output_dump_ib3; 
else:
    path_dropbox_main = path_dropbox_casa; 
    path_fasta_main = path_fasta_casa; 
    path_output_dump_main = path_output_dump_casa; 
path_bed_main = path_dropbox_main + 'Doctorado\\3-Genes transactivados rio abajo\\0-Fuentes\\Papers ChIP-seq\\'; 
path_rnaseq_main = path_dropbox_main + 'Doctorado\\3-Genes transactivados rio abajo\\1-Interpretacion de datos ChIP-seq\\'; 
path_in_main = path_output_dump_main + ''; ### COMPLETAR
path_out_main = path_output_dump_main + 'SitiosUnionNKX25\\'; 
# Path que dependen de la especie
path_pwm_human = path_dropbox_main + 'Doctorado\\3-Genes transactivados rio abajo\\1-Interpretacion de datos ChIP-seq\\PWM_human\\'; 
path_pwm_mouse = path_dropbox_main + 'Doctorado\\3-Genes transactivados rio abajo\\1-Interpretacion de datos ChIP-seq\\PWM_mouse\\'; 



#################################### FUNCIONES ####################################


def abrir_L_pssm(L_nom_arch, path_arch='', pseudocounts=0.5):
    # Extrae una lista de matrices de counts y devuelve sus PSSM correspondientes

    # Inicializo la lista que se devuelve
    L_pssm = []; 
    # Recorro L_nom_arch
    for nom_arch in L_nom_arch:
        # Uso abrir_pssm() con solo_pssm=True para extraer los PSSM
        L_pssm.append(abrir_pssm(nom_arch, path_arch=path_arch, solo_pssm=True, pseudocounts=pseudocounts)); 
    return L_pssm


def abrir_pssm(nom_arch, path_arch='', solo_pssm=True, pseudocounts=0.5):
    # Extrae la matriz de counts de nom_arch y la pasa a PSSM

    # Extraigo el motif con extraer_motif_pcm()
    m = extraer_motif_pcm(nom_arch, path_arch=path_arch); 
    # Reviso que m no sea lista vacia y genero pwm y pssm
    if m != '':
        pwm = m.counts.normalize(pseudocounts=pseudocounts); 
        pssm = pwm.log_odds(); 
    else:
        pwm = ''; 
        pssm = ''; 
    # Devuelvo solo pssm por default
    if solo_pssm:
        return pssm
    # Si solo_pssm==False, devuelvo m y pwm tambien
    else:
        return m, pwm, pssm


def buscar_sitios_gen(gen_id, L_su, genoma, nom_genoma, path_fasta='', dist=100000):
    # Busca todos los sitios de union de L_su dados alrededor del +1 de un gen a distancia dist y devuelve una lista con los hits

    # Inicializo seq_data
    seq_data_usado = seq_data(nom_genoma, genome_element=genoma, path_fasta=path_fasta); 
    # Agarro el elemento gen con gen_id
    gen_elem = genoma.gene_by_id(gen_id); 
    # Busco el start y chr_n
    gen_start = int(gen_elem.start); 
    chr_n = seq_data_usado._obtener_chr(gen_elem.contig); 
    # Defino rango usado con dist
    rango_usado = [gen_start-dist, gen_start+dist]; 
    # Busco la secuencia con seq_data._consulta_secuencia_fasta()
    seq_alr = str(seq_data_usado._consulta_secuencia_fasta(chr_n, rango_usado[0], rango_usado[1])); 
    # Inicializo la lista de sitios de union
    L_hits = []; 
    # Recorro cada sitio de union
    for curr_sitio in L_su:
        # Busco todos los sitios de union que hayan para curr_sitio
        L_curr_sitio_encontrados = seq_data_usado._buscar_SU_en_seq(curr_sitio, seq_alr); 
        # Recorro cada hit
        for curr_sitio_encontrado in L_curr_sitio_encontrados:
            # Defino pos_ini y pos_end respecto a chr_n
            curr_pos_ini = curr_sitio_encontrado[0]+rango_usado[0]-1; 
            curr_pos_end = curr_sitio_encontrado[1]+rango_usado[0]-1; 
            # Defino la distancia a gen_start
            if curr_pos_end < gen_start:
                # Si el sitio esta antes de gen_start, hago pos_end-start
                dist_pos = curr_pos_end-gen_start; 
            elif curr_pos_ini > gen_start:
                # Si el sitio esta despues de gen_start, hago pos_ini-start+1
                dist_pos = curr_pos_ini-gen_start+1; 
            else:
                # Si el sitio se solapa con el +1, la distancia es 0
                dist_pos = 0; 
            # Creo la fila con los datos ordenados
            SU_encontrado = [curr_pos_ini, curr_pos_end, '', 'nkx25_list', str(curr_sitio), curr_sitio_encontrado[2], int(dist_pos)]; 
            L_hits.append(SU_encontrado[:]); 
    return L_hits, chr_n


def buscar_L_pssm_gen(gen_id, L_pssm, L_nom_pssm, genoma, nom_genoma, path_fasta='', dist=100000, score_lim=0.9):
    # Busca todos los sitios de union de L_pssm alrededor del +1 de un gen a distancia dist y devuelve una lista con los hits

    # Inicializo seq_data
    seq_data_usado = seq_data(nom_genoma, genome_element=genoma, path_fasta=path_fasta); 
    # Agarro el elemento gen con gen_id
    gen_elem = genoma.gene_by_id(gen_id); 
    # Busco el start y chr_n
    gen_start = int(gen_elem.start); 
    chr_n = seq_data_usado._obtener_chr(gen_elem.contig); 
    # Defino rango usado con dist
    rango_usado = [gen_start-dist, gen_start+dist]; 
    # Busco la secuencia con seq_data._consulta_secuencia_fasta()
    seq_alr = str(seq_data_usado._consulta_secuencia_fasta(chr_n, rango_usado[0], rango_usado[1])); 
    # Inicializo la lista de sitios de union
    L_hits = []; 
    # Recorro L_pssm
    for p in range(len(L_pssm)):
        curr_pssm = L_pssm[p]; 
        # Defino score_cutoff para curr_pssm
        score_cutoff = curr_pssm.max * score_lim; 
        # Defino el largo de pssm
        len_pssm = len(curr_pssm.consensus); 
        # Busco todos los sitios de union que hayan para pssm
        L_sitios_curr_pssm = seq_data_usado._buscar_PSSM_en_seq(curr_pssm, seq_alr, score_cutoff=score_cutoff, pos_ini_ref=rango_usado[0]-1, return_forward=True); # pos_ini-1 por como suma _buscar_PSSM_en_seq()
        # Recorro los sitios encontrados
        for sitio_pssm in L_sitios_curr_pssm:
            # Defino posicion inicial/final y score de sitio_pssm
            pos_ini_sitio = sitio_pssm[0]+1; 
            pos_end_sitio = sitio_pssm[0]+len_pssm; 
            score_pssm = sitio_pssm[1]; 
            # Defino la distancia a gen_start
            if pos_end_sitio < gen_start:
                # Si el sitio esta antes de gen_start, hago pos_end-start
                dist_pos = pos_end_sitio-gen_start; 
            elif pos_ini_sitio > gen_start:
                # Si el sitio esta despues de gen_start, hago pos_ini-start+1
                dist_pos = pos_ini_sitio-gen_start+1; 
            else:
                # Si el sitio se solapa con el +1, la distancia es 0
                dist_pos = 0; 
            # Defino la informacion importante en curr_hit
            curr_hit = [int(pos_ini_sitio), int(pos_end_sitio), round(score_pssm, 2), str(L_nom_pssm[p]), str(sitio_pssm[2]), bool(sitio_pssm[3]), int(dist_pos)]; 
            # Agrego curr_hit a L_hits
            L_hits.append(curr_hit[:]); 
    return L_hits


def extraer_motif_pcm(nom_arch, path_arch=''):
    # Abre archivo con motif en formato pfm-four-columns y devuelve el primer motif

    # Defino dir_arch
    if path_arch != '':
        dir_arch = os.path.join(path_arch, nom_arch); 
    else:
        dir_arch = str(nom_arch); 
    # Abro pcm con motifs.parse
    handle = open(dir_arch); 
    record = motifs.parse(handle, 'pfm-four-columns'); 
    handle.close(); 
    # Veo que record tenga un solo motif
    if len(record) == 1:
        m = record[0]; 
    elif len(record) > 1:
        print('WARNING: ' + nom_arch + ' tiene ' + str(len(record)) + ' motifs. Se devuelve solo el primero.')
        m = record[0]; 
        print('Motifs en record:')
        for i in record:
            print(i.name)
            print(i.counts)
    else:
        print('ERROR: No se encontro ningun motif en ' + nom_arch + '. Se devuelve string vacio.')
        m = ''; 
    return m


def gen_id_gen_name(nombre_gen, genoma):
    # Devuelve un id de ensembl asociado a nombre_gen en genoma

    L_id_gen = genoma.gene_ids_of_gene_name(nombre_gen); 
    if len(L_id_gen) > 1:
        print('Mas de un match para ' + str(nombre_gen) + ', se devuelve el primero.')
        print(L_id_gen)
        ret = L_id_gen[0]; 
    elif len(L_id_gen) == 0:
        print('No se encontro ' + str(nombre_gen) + ' en el genoma dado.')
        ret = ''; 
    else:
        ret = L_id_gen[0]; 
    return ret


def guardar_hits(L_hits, L_head, nom_arch, path_arch='', ext='.csv', sep=';'):
    # Guarda L_hits en archivo nom_arch dentro de path_arch con titulos L_headers

    # Defino la direccion del archivo en base a path_arch y nom_arch
    dirarch = os.path.join(path_arch, nom_arch + ext); 
    # Creo el archivo
    with open(dirarch, 'w') as F_out:
        print('Archivo ' + nom_arch + ext + ' creado.')
    # Lo vuelvo a abrir en modo append
    with open(dirarch, 'a') as F_out:
        # Creo el titulo en base a L_head
        if len(L_head) > 0:
            str_head = ''; 
            # Recorro L_head
            for h in L_head:
                str_head = str_head + str(h) + sep; 
            # Elimino la ultima ocurrencia de sep
            str_head = str_head.rstrip(sep); 
            # Agrego end of line y guardo en F_out
            F_out.write(str_head + '\n'); 
        # Recorro L_hits
        for curr_hit in L_hits:
            curr_str = ''; 
            # Recorro L_out
            for i in curr_hit:
                curr_str = curr_str + str(i) + sep; 
            # Elimino la ultima ocurrencia de sep
            curr_str = curr_str.rstrip(sep); 
            # Agrego end of line y guardo en F_out
            F_out.write(curr_str + '\n'); 
    return L_hits


def seq_alrededor_ini_gen(gen_id, genoma, nom_genoma, path_fasta='', dist=100000):
    # Busca la secuencia alrededor del +1 de gen_id en genoma, hasta dist a ambos lados

    # Inicializo seq_data
    seq_data_usado = seq_data(nom_genoma, genome_element=genoma, path_fasta=path_fasta); 
    # Agarro el elemento gen con gen_id
    gen_elem = genoma.gene_by_id(gen_id); 
    # Busco el start y chr_n
    gen_start = int(gen_elem.start); 
    chr_n = seq_data_usado._obtener_chr(gen_elem.contig); 
    # Defino rango usado con dist
    rango_usado = [gen_start-dist, gen_start+dist]; 
    # Busco la secuencia con seq_data._consulta_secuencia_fasta()
    seq_alrededor = str(seq_data_usado._consulta_secuencia_fasta(chr_n, rango_usado[0], rango_usado[1])); 
    return seq_alrededor


def _main():
    # Funcion para probar funciones en ejecucion del archivo

    dist_usada = 100000; 
    default_headers = ['chr_n', 'pos_ini', 'pos_end', 'score', 'nom_tf', 'seq', 'forward', 'dist_pos']; 

    # Human
    print('>Iniciando procesamiento en humano')
    L_pssm_human = abrir_L_pssm(L_arch_pwm_humano, path_arch=path_pwm_human); 
    print('Matrices PSSM recolectadas')
    nkx25_id_human = gen_id_gen_name(nkx25_human, hg19); 
    hits_pssm_human = buscar_L_pssm_gen(nkx25_id_human, L_pssm_human, L_nombres_pssm_humano, hg19, 'hg19', path_fasta=path_fasta_main, dist=dist_usada); 
    print('Hits PSSM encontrados')
    hits_L_human, chr_n_human = buscar_sitios_gen(nkx25_id_human, L_sitios, hg19, 'hg19', path_fasta=path_fasta_main, dist=dist_usada); 
    print('Hits L_sitios encontrados')
    hits_human = hits_pssm_human + hits_L_human; 
    for i in range(len(hits_human)):
        hits_human[i] = [chr_n_human] + hits_human[i]; 
    hits_human = guardar_hits(hits_human, default_headers, 'Sitios_NKX25_human', path_arch=path_out_main); 
    print('Hits guardados')

    # Mouse
    print('>Iniciando procesamiento en raton')
    L_pssm_mouse = abrir_L_pssm(L_arch_pwm_raton, path_arch=path_pwm_mouse); 
    print('Matrices PSSM recolectadas')
    nkx25_id_mouse = gen_id_gen_name(nkx25_mouse, mm9); 
    hits_pssm_mouse = buscar_L_pssm_gen(nkx25_id_mouse, L_pssm_mouse, L_nombres_pssm_raton, mm9, 'mm9', path_fasta=path_fasta_main, dist=dist_usada); 
    print('Hits PSSM encontrados')
    hits_L_mouse, chr_n_mouse = buscar_sitios_gen(nkx25_id_mouse, L_sitios, mm9, 'mm9', path_fasta=path_fasta_main, dist=dist_usada); 
    print('Hits L_sitios encontrados')
    hits_mouse = hits_pssm_mouse + hits_L_mouse; 
    for i in range(len(hits_mouse)):
        hits_mouse[i] = [chr_n_mouse] + hits_mouse[i]; 
    hits_mouse = guardar_hits(hits_mouse, default_headers, 'Sitios_Nkx25_mouse', path_arch=path_out_main); 
    print('Hits guardados')

    return ''


def _probar_hits(L_hits, chr_n, genoma, nom_genoma, path_fasta=''):
    # 

    # Inicializo seq_data
    seq_data_usado = seq_data(nom_genoma, genome_element=genoma, path_fasta=path_fasta); 
    # Recorro L_hits
    for curr_hit in L_hits:
        if curr_hit[2]:
            print('Forward: ' + str(seq_data_usado._consulta_secuencia_fasta(chr_n, curr_hit[0], curr_hit[1])))
        else:
            print('Reverse: ' + str(seq_data_usado._consulta_secuencia_fasta(chr_n, curr_hit[0], curr_hit[1])))
    return



#################################### RESULTADOS ###################################

output_dump = []; 

if __name__=='__main__':
    output_dump.append(_main()); 

