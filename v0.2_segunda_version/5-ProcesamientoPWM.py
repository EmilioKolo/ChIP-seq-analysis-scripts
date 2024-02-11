
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

# Importo analisischip
import sys
path_analisischip_casa = 'C:\\Users\\Admin\\Documents\\Scripts\\Git\\AnalisisSecuenciasChIP\\'; 
path_analisischip_ib3 = 'C:\\Users\\emili\\Documents\\Archivos nube\\Git\\analisissecuenciaschip\\'; 
path_analisischip_main = path_analisischip_ib3; # Cambiar aca
sys.path.append(path_analisischip_main); 
from analisischip.seq import seq_data

# Genomas de referencia
hg19 = EnsemblRelease(75, species='human'); 
hg38 = EnsemblRelease(102, species='human'); 
mm9 = EnsemblRelease(54, species='mouse'); 
mm10 = EnsemblRelease(102, species='mouse'); 


'''
Scripts pensados para buscar Position Weight Matrices (PWMs) en secuencias de ADN usando Bio.motifs

FALTA:
X Buscar ejemplos de PWM/PSSM en MEME y descargarlos
    X Usar resultados de ChIP-seq total
    X Buscar consenso de motif para NKX2-5, GATA1/3/4 y TBX20/MEIS1
X Cargar motifs a elementos Bio.motifs
    * https://biopython-tutorial.readthedocs.io/en/latest/notebooks/14%20-%20Sequence%20motif%20analysis%20using%20Bio.motifs.html#
        * m = motifs.create(instances)
    * MEME: https://biopython-tutorial.readthedocs.io/en/latest/notebooks/14%20-%20Sequence%20motif%20analysis%20using%20Bio.motifs.html#MEME
        * handle = open("meme.dna.oops.txt"); record = motifs.parse(handle, "meme"); handle.close()
            * record es LISTA de elementos motif (ver len(record) y probar los distintos motif)
- Usar motifs para buscar sus ubicaciones en secuencias
- Acoplar scripts a clasificacion de peaks para buscar motifs en peaks de ChIP-seq
'''

#################################### VARIABLES ####################################

path_output_dump_ib3 = 'X:\\Output_dump\\'; 
path_output_dump_casa = 'D:\\Archivos doctorado\\Output_dump\\'; 

path_seq_fasta_casa = path_output_dump_casa + 'FastaMEME\\DistintoLargo\\'; 
path_seq_fasta_ib3 = path_output_dump_ib3 + 'FastaMEME\\DistintoLargo\\'; 
seq_fasta_dupays = 'Dupays\\DupaysTotal.fasta'; 
seq_fasta_anderson = 'Anderson\\AndersonTotal.fasta'; 

path_pwm_casa = 'D:\\Users\\Admin\\Dropbox\\Doctorado\\3-Genes transactivados rio abajo\\1-Interpretacion de datos ChIP-seq\\PWM\\'; 
path_pwm_ib3 = 'C:\\Users\\emili\\Dropbox\\Doctorado\\3-Genes transactivados rio abajo\\1-Interpretacion de datos ChIP-seq\\PWM\\'; 

L_arch_hocomoco_pcm = ['NKX25_MOUSE.H11MO.0.A.pcm', 'TBX20_MOUSE.H11MO.0.C.pcm', 'MEIS1_MOUSE.H11MO.1.A.pcm', 'TGIF1_MOUSE.H11MO.0.A.pcm', 
                       'HAND1_MOUSE.H11MO.0.C.pcm', 'MAF_MOUSE.H11MO.1.A.pcm', 'GATA1_MOUSE.H11MO.1.A.pcm', 'GATA6_MOUSE.H11MO.0.A.pcm', 'GATA4_MOUSE.H11MO.0.A.pcm']; 

L_arch_hocomoco_words = ['NKX25_MOUSE.H11MO.0.A.words', 'TBX20_MOUSE.H11MO.0.C.words', 'MEIS1_MOUSE.H11MO.1.A.words', 'TGIF1_MOUSE.H11MO.0.A.words', 
                         'HAND1_MOUSE.H11MO.0.C.words', 'MAF_MOUSE.H11MO.1.A.words', 'GATA1_MOUSE.H11MO.1.A.words', 'GATA6_MOUSE.H11MO.0.A.words', 'GATA4_MOUSE.H11MO.0.A.words']; 

'''0-Links
Nkx2-5:
https://hocomoco11.autosome.org/motif/NKX25_MOUSE.H11MO.0.A
NKX25_MOUSE.H11MO.0.A.pcm
NKX25_MOUSE.H11MO.0.A.words

Tbx20:
https://hocomoco11.autosome.org/motif/TBX20_MOUSE.H11MO.0.C
TBX20_MOUSE.H11MO.0.C.pcm
TBX20_MOUSE.H11MO.0.C.words

Meis1:
https://hocomoco11.autosome.org/motif/MEIS1_MOUSE.H11MO.1.A
MEIS1_MOUSE.H11MO.1.A.pcm
MEIS1_MOUSE.H11MO.1.A.words

Tgif1:
https://hocomoco11.autosome.org/motif/TGIF1_MOUSE.H11MO.0.A
TGIF1_MOUSE.H11MO.0.A.pcm
TGIF1_MOUSE.H11MO.0.A.words

Hand1:
https://hocomoco11.autosome.org/motif/HAND1_MOUSE.H11MO.0.C
HAND1_MOUSE.H11MO.0.C.pcm
HAND1_MOUSE.H11MO.0.C.words

Maf:
https://hocomoco11.autosome.org/motif/MAF_MOUSE.H11MO.1.A
MAF_MOUSE.H11MO.1.A.pcm
MAF_MOUSE.H11MO.1.A.words

Gata1:
https://hocomoco11.autosome.org/motif/GATA1_MOUSE.H11MO.1.A
GATA1_MOUSE.H11MO.1.A.pcm
GATA1_MOUSE.H11MO.1.A.words

Gata6:
https://hocomoco11.autosome.org/motif/GATA6_MOUSE.H11MO.0.A
GATA6_MOUSE.H11MO.0.A.pcm
GATA6_MOUSE.H11MO.0.A.words

Gata4:
https://hocomoco11.autosome.org/motif/GATA4_MOUSE.H11MO.0.A
GATA4_MOUSE.H11MO.0.A.pcm
GATA4_MOUSE.H11MO.0.A.words
'''

# Variables main()

path_pwm_main = path_pwm_ib3; 
path_seq_fasta_main = path_seq_fasta_ib3; 

#################################### FUNCIONES ####################################


def buscar_motif_en_bed(bed_arch, matriz_su, genome_name, path_bed='', cutoff=3.0):
    # Funcion para buscar hits de una matriz en los peaks de un archivo .bed con resultados de ChIP-seq

    ### FALTA
    # Abrir .bed
    # Conseguir secuencias de peaks
        # Tratar de guardarlas en .csv? para acelerar esto
    # Buscar motif con funciones armadas previamente
        # Ver 3-ClasificacionPeaks.py
            # Funciones pipeline_bed(), procesar_bed() y procesar_peak()
    ###
    return


def buscar_motif_en_seq(seq_ref, pssm, len_m, cutoff=3.0):
    # Devuelve las posiciones en seq_ref donde pssm da score mayor a cutoff
    # Se puede usar seq_ref[pos:pos+len_m] para obtener las secuencias de los hits

    # Inicializo la lista de sitios de union
    L_su = []; 
    # Recorro seq_ref con pssm.search()
    for pos, score in pssm.search(seq_ref, threshold=cutoff):
        # Inicializo la lista que se agrega a L_su
        curr_su = []; 
        # Defino curr_seq en base a pos y len_m
        curr_seq = seq_ref[pos:pos+len_m]; 
        # Agrego los datos relevantes a curr_su
        curr_su.append(pos); 
        curr_su.append(score); 
        curr_su.append(str(curr_seq)); 
        # Agrego curr_su a L_su
        L_su.append(curr_su[:]); 
    return L_su


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
        print('WARNING: record tiene ' + str(len(record)) + ' motifs. Se devuelve solo el primero.')
        m = record[0]; 
        print('Motifs en record:')
        for i in record:
            print(i.name)
            print(i.counts)
    else:
        print('ERROR: No se encontro ningun motif en record. Se devuelve string vacio')
        m = ''; 
    return m


def extraer_motif_pwm_pssm(nom_arch, path_arch='', verbose=False):
    # Extrae el motif de nom_arch con extraer_motif_pcm() y genera pwm y pssm

    ### Display
    if verbose:
        print('### Iniciando procesamiento de archivo ' + nom_arch)
    ###
    # Extraigo el motif con extraer_motif_pcm()
    m = extraer_motif_pcm(nom_arch, path_arch=path_arch); 
    ### Display
    if verbose:
        if m != '':
            print('Motif extraido:')
            print(m.name)
            print(m.consensus)
            print(m.counts)
    ###
    if m != '':
        pwm = m.counts.normalize(pseudocounts=0.5); 
        pssm = pwm.log_odds(); 
    else:
        print('ERROR: Ningun motif enontrado.')
        pwm = ''; 
        pssm = ''; 
    ### Display
    if verbose:
        print('PWM:')
        print(pwm)
        print(pwm.consensus)
        print('PSSM:')
        print(pssm)
        print(pssm.consensus)
    ###
    return m, pwm, pssm


def _main():
    # Funcion para probar funciones en ejecucion del archivo

    # Defino secuencia de testing
    test_seq = Seq.Seq("TACTTAAGTGCTGCATTACAACCCAAGCATTA"); 
    # test_seq = "TACTTAAGTGCTGCATTACAACCCAAGCATTA"; 
    # Corro extraer_motif_pcm con varios motif
    for arch_pcm in L_arch_hocomoco_pcm[:2]: 
        print('### Iniciando procesamiento de archivo ' + arch_pcm)
        m, pwm, pssm = extraer_motif_pwm_pssm(arch_pcm, path_arch=path_pwm_main, verbose=False); 
        # Copiado de tutorial
        for position, score in pssm.search(test_seq, threshold=0.0):
            print('# pssm:')
            #print(pssm)
            print('len(pssm):')
            print(len(pssm))
            print('consensus pssm:')
            print(pssm.consensus)
            print('len(consensus pssm):')
            print(len(pssm.consensus))
            print('# pwm:')
            #print(pwm)
            print('len(pwm):')
            print(len(pwm))
            print('# m')
            #print(m)
            print('len(m):')
            print(len(m))
            print('pssm consensus:')
            print(pssm.consensus)
            print("Position %d: score = %5.3f" % (position, score))
            print(test_seq[position:position+len(m)])
        print()
    return ''


#################################### RESULTADOS ###################################

output_dump = []; 

if __name__=='__main__':
    output_dump.append(_main()); 

