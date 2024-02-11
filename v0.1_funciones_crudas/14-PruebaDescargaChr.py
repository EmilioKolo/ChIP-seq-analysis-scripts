
import os
import time
from Bio import Entrez
from Bio import SeqIO
from pyensembl import EnsemblRelease

#################################### VARIABLES ####################################

Entrez.email = 'ekolomenski@gmail.com';
Entrez.api_key = '9bccbd4949f226a5930868df36c211e8b408';
mm9 = EnsemblRelease(67, species='mouse');
hg19 = EnsemblRelease(102, species='human');
GRCm38 = EnsemblRelease(102, species='mouse');

path_out = 'D:\\Archivos doctorado\\Genomas\\';
#path_out = 'X:\\Genomas\\';


###################################### CLASES #####################################

class seq_data(object):
    '''
    Clase para cargar datos de secuencia en un genoma dado
    Descarga secuencia de cromosomas del genoma dado y registra posiciones
    '''
    def __init__(self, genome):
        # Almaceno datos de secuencia en archivos .fasta al crear la clase
        # Cromosomas seleccionados son todos los numericos mas X, Y y M/MT
        # Devuelve rangos de secuencias para esos cromosomas
        return None


#################################### FUNCIONES ####################################


def consulta_entrez(search, database, return_full=False, return_value='IdList'):
    # Hace una consulta a Entrez con los datos dados
    handle = Entrez.esearch(db=database, term=search, retmax=1000);
    record = Entrez.read(handle);
    if return_full:
        ret = record;
    else:
        ret = record[return_value];
    return ret


def consulta_secuencia_chr(chr_id):
    # Devuelve la secuencia completa de un cromosoma en base al ID de base de datos nucleotide
    nuc_id = consulta_entrez(chr_id, database='nucleotide');
    if len(nuc_id) > 1:
        print('ERROR: Mas de un resultado para el id dado. ' + str(nuc_id))
    handle = Entrez.efetch(db="nucleotide", rettype="gbwithparts", retmode="text", id=nuc_id[0]);
    records = SeqIO.read(handle, 'genbank');
    return records


def consulta_secuencia_chr_retries(chr_id, retries=5):
    # Devuelve la secuencia completa de un cromosoma en base al ID de base de datos nucleotide
    # Repite retries veces si falla la extraccion

    # Chequeo para retries < 1
    if retries < 1:
        retries = 1;
    # Inicializo las variables del ciclo while
    encontrado = False;
    tries = 0;
    # Repito maximo retries veces o hasta que se encuentre
    while not encontrado and tries < retries:
        try:
            seq_chr = consulta_secuencia_chr(chr_id);
            encontrado = True;
        except:
            tries += 1;
            print('Fallo intento ' + str(tries) + '.')
    # Si no se encuentra despues de retries veces, se devuelve texto vacio y mensaje de error
    if not encontrado:
        print('ERROR: No se encontro la secuencia despues de ' + str(tries) + ' intentos. Se devuelve secuencia vacia.')
        seq_chr = '';
    return seq_chr


def leer_fasta(nom_arch, path_arch, ext='fasta'):
    # Abre archivo fasta nom_arch con extension ext en path_arch y devuelve las secuencias e info

    # Inicializo la variable que se devuelve
    L_out = [];

    # Defino dir_arch en base a nom_arch
    dir_arch = os.path.join(path_arch, nom_arch + '.' + ext);

    # Abro el archivo
    arch_fasta = SeqIO.parse(dir_arch, ext);

    # Contador para revisar si hay mas de un registro
    i = 0;
    # Reviso si tiene mas de una secuencia
    for record in arch_fasta:
        i = i + 1;
        L_out.append(record);
    if i == 1:
        L_out = L_out[0];
    else:
        print('Se devuelve mas de un registro')
    return L_out


def guardar_secuencias_dict_chr(dict_chr:dict, path_output, nombre_arch, ext='fasta'):
    # Busca y guarda las secuencias para todos los elementos en dict_chr

    # Repito todo para cada key (chr_n) en dict_chr
    for chr_n in dict_chr.keys():
        print('Iniciando busqueda de secuencia para ' + str(chr_n))
        # Primero extraigo la secuencia de chr_n usando retries
        seq_chr = consulta_secuencia_chr_retries(dict_chr[chr_n], retries=10);
        # Chequeo que se haya encontrado texto
        if len(seq_chr) > 0:
            print('Secuencia para ' + str(chr_n) + ' obtenida.')
            SeqIO.write(seq_chr, os.path.join(path_output, nombre_arch + '_' + chr_n + '.' + ext), ext);
            print('Secuencia para ' + str(chr_n) + ' guardada.')
    return dict_chr


def _main_test():
    # Funcion para probar funciones en ejecucion del archivo

    # Inicializo la variable que se devuelve
    L_out = [];

    # Consulta a mano de cromosomas
    dict_chr_hg19 = {'chr1':'NC_000001.10', 'chr2':'NC_000002.11', 'chr3':'NC_000003.11', 'chr4':'NC_000004.11', 'chr5':'NC_000005.9', 
                     'chr6':'NC_000006.11', 'chr7':'NC_000007.13', 'chr8':'NC_000008.10', 'chr9':'NC_000009.11', 'chr10':'NC_000010.10', 
                     'chr11':'NC_000011.9', 'chr12':'NC_000012.11', 'chr13':'NC_000013.10', 'chr14':'NC_000014.8', 'chr15':'NC_000015.9', 
                     'chr16':'NC_000016.9', 'chr17':'NC_000017.10', 'chr18':'NC_000018.9', 'chr19':'NC_000019.9', 'chr20':'NC_000020.10', 
                     'chr21':'NC_000021.8', 'chr22':'NC_000022.10', 'chrX':'NC_000023.10', 'chrY':'NC_000024.9'};
    dict_chr_hg38 = {'chr1':'NC_000001.11', 'chr2':'NC_000002.12', 'chr3':'NC_000003.12', 'chr4':'NC_000004.12', 'chr5':'NC_000005.10', 
                     'chr6':'NC_000006.12', 'chr7':'NC_000007.14', 'chr8':'NC_000008.11', 'chr9':'NC_000009.12', 'chr10':'NC_000010.11', 
                     'chr11':'NC_000011.10', 'chr12':'NC_000012.12', 'chr13':'NC_000013.11', 'chr14':'NC_000014.9', 'chr15':'NC_000015.10', 
                     'chr16':'NC_000016.10', 'chr17':'NC_000017.11', 'chr18':'NC_000018.10', 'chr19':'NC_000019.10', 'chr20':'NC_000020.11', 
                     'chr21':'NC_000021.9', 'chr22':'NC_000022.11', 'chrX':'NC_000023.11', 'chrY':'NC_000024.10'};
    dict_chr_mm9 = {'chr1':'NC_000067.5', 'chr2':'NC_000068.6', 'chr3':'NC_000069.5', 'chr4':'NC_000070.5', 'chr5':'NC_000071.5', 
                    'chr6':'NC_000072.5', 'chr7':'NC_000073.5', 'chr8':'NC_000074.5', 'chr9':'NC_000075.5', 'chr10':'NC_000076.5', 
                    'chr11':'NC_000077.5', 'chr12':'NC_000078.5', 'chr13':'NC_000079.5', 'chr14':'NC_000080.5', 'chr15':'NC_000081.5', 
                    'chr16':'NC_000082.5', 'chr17':'NC_000083.5', 'chr18':'NC_000084.5', 'chr19':'	NC_000085.5', 
                    'chrX':'NC_000086.6', 'chrY':'NC_000087.6'};
    print('>Empezando descarga de cromosomas para hg19.')
    dict_chr_hg19 = guardar_secuencias_dict_chr(dict_chr_hg19, path_out, 'hg19', ext='fasta');
    print('FIN descarga de cromosomas para hg19.')
    print()
    print('>Empezando descarga de cromosomas para mm9.')
    dict_chr_mm9 = guardar_secuencias_dict_chr(dict_chr_mm9, path_out, 'mm9', ext='fasta');
    print('FIN descarga de cromosomas para mm9.')
    print()
    print('>Empezando descarga de cromosomas para hg38.')
    dict_chr_hg38 = guardar_secuencias_dict_chr(dict_chr_hg38, path_out, 'hg38', ext='fasta');
    print('FIN descarga de cromosomas para hg38.')
    print()

    L_out.append(dict_chr_hg19);
    L_out.append(dict_chr_mm9);
    L_out.append(dict_chr_hg38);

    return L_out


def _main_test2():
    # Funcion para probar funciones en ejecucion del archivo

    # Inicializo la variable que se devuelve
    L_out = [];

    # Pruebo leer_fasta para hg19_chr1.fasta
    #L_out = leer_fasta('hg19_chr1', path_out, ext='fasta');
    #print(L_out)
    #print(len(L_out.seq))

    # Pruebo leer_fasta para todos los cromosomas de todos los genomas
    for i in range(3):
        genome_name = ['hg19', 'mm9', 'hg38'][i];
        print('>Inicializando lectura de genoma ' + genome_name)
        # Defino la lista de chr que comparten los 3 genomas
        L_chr = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 
        'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 
        'chr19', 'chrX', 'chrY'];
        # Agrego chr para genomas humanos
        if genome_name != 'mm9':
            L_chr.append('chr20');
            L_chr.append('chr21');
            L_chr.append('chr22');
        # Recorro la lista de chr_n
        for j in range(len(L_chr)):
            chr_n = L_chr[j];
            record_name = genome_name + '_' + chr_n;
            chr_record = leer_fasta(record_name, path_out, ext='fasta');
            L_out.append(chr_record)
    for i in L_out:
        print(i)
    ##### FALTA:
    ## Inicializar seq_data y cargar rangos + bed (ver analisissecuenciaschip para carga de rangos)
    ## Agarrar secuencias en seq_data de los archivos .fasta
    ## Hacer busquedas de sitios de union en los rangos registrados

    return L_out

##################################### PRUEBAS #####################################

output_dump = [];

if __name__=='__main__':
    #output_dump.append(_main_test());
    output_dump.append(_main_test2());


##################################### VIEJOS ######################################


def consulta_databases():
    # Devuelve las bases de datos en Entrez
    handle = Entrez.einfo();
    record = Entrez.read(handle);
    return record['DbList']
#['pubmed', 'protein', 'nuccore', 'ipg', 'nucleotide', 'structure', 'genome', 'annotinfo', 'assembly', 'bioproject', 
# 'biosample', 'blastdbinfo', 'books', 'cdd', 'clinvar', 'gap', 'gapplus', 'grasp', 'dbvar', 'gene', 'gds', 'geoprofiles', 
# 'homologene', 'medgen', 'mesh', 'ncbisearch', 'nlmcatalog', 'omim', 'orgtrack', 'pmc', 'popset', 'proteinclusters', 
# 'pcassay', 'protfam', 'pccompound', 'pcsubstance', 'seqannot', 'snp', 'sra', 'taxonomy', 'biocollections', 'gtr']


def consulta_genoma_completo(search, database='nucleotide', return_full=False, return_value='IdList', property='complete genome'):
    # Hace consulta a Entrez incluyendo property
    handle = Entrez.esearch(db=database, term=search, retmax=1000, property=property);
    record = Entrez.read(handle);
    if return_full:
        ret = record;
    else:
        ret = record[return_value];
    return ret


'''
#############################
# Retrieve NCBI Data Online #
#############################

genomeAccessions = ['NC_000913', 'NC_002695', 'NC_011750', 'NC_011751', 'NC_017634', 'NC_018658']
search           = " ".join(genomeAccessions)
handle           = Entrez.read(Entrez.esearch(db="nucleotide", term=search, retmode="xml"))
genomeIds        = handle['IdList']
records          = Entrez.efetch(db="nucleotide", id=genomeIds, rettype="gb", retmode="text")

###############################
# Generate Genome Fasta files #
###############################

sequences   = []  # store your sequences in a list
headers     = []  # store genome names in a list (db_xref ids)

for i,record in enumerate(records):

    file_out = open("genBankRecord_"+str(i)+".gb", "w")    # store each genomes .gb in separate files
    file_out.write(record.read())
    file_out.close()

    genomeGenbank   = SeqIO.read("genBankRecord"+str(i)+".gb", "genbank")  # parse in the genbank files
    header         = genome.features[0].qualifiers['db_xref'][0]          # name the genome using db_xfred ID
    sequence       = genome.seq.tostring()                                # obtain genome sequence

    headers.append('>'+header)  # store genome name in list                                     
    sequences.append(sequence)  # store sequence in list

    fasta_out = open("genome"+str(i)+".fasta","w")     # store each genomes .fasta in separate files
    fasta_out.write(header)    # >header ... followed by:
    fasta_out.write(sequence)  # sequence ... 
    fasta_out.close()          # close that .fasta file and move on to next genome
records.close()
'''
