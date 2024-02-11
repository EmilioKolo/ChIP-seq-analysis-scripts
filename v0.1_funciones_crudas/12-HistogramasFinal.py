
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
###### VER QUE ESTE referencias_funciones.py
# from referencias_funciones import ConsultaSecuencia, IDchr 

mm9 = EnsemblRelease(67, species='mouse');
hg19 = EnsemblRelease(102, species='human');
GRCm38 = EnsemblRelease(102, species='mouse');


'''
Script que genera histogramas de todo tipo para la tesis

FALTA:
~ Recopilar datos necesarios
    X Todos los sitios de union en el genoma
    ~ Ubicacion de los sitios de union en peaks de ChIP-seq o no
    ? Secuencias (peaks ChIP-seq y -1000/+1000 en todo el genoma)
- Seleccionar sets de datos
    - Todos los sitios (-1000/+1000 en mm9 y hg19)
    - Solo sitios en ChIP-seq (3 papers mm9 y 1 paper hg19)
    - Preparar otros criterios de seleccion
- Generar graficos para los sets de datos
    - Por secuencia de union + todo junto
    - Por paper en ChIP-seq
    ~ Total (lo que ya tengo)
'''

#################################### VARIABLES ####################################

L_sitios_union_nkx25 = ['GCAAGTG', 'GGAAGTG', 'GAAAGTG', 'ATAAGTG', 'GTAAGTG',
                        'CTAAGTG', 'TCAAGTG', 'TGAAGTG', 'TAAAGTG', 'TTAAGTG'];
rango_sitio_genoma = (-1000, 1000);

nombre_archivo_mm9 = 'SitiosUnionFinal_mm9';
nombre_archivo_hg19 = 'SitiosUnionFinal_hg19';
path_archivo ='.\\sitios de union FINAL\\';

curr_path = os.path.dirname(__file__);

# Variables main()

L_sitios_main = L_sitios_union_nkx25;
rango_sitios_main = rango_sitio_genoma;
genoma_main = 'mm9';

nombre_archivo_main = nombre_archivo_mm9;
dir_archivo_main = os.path.join(curr_path, str(path_archivo) + str(nombre_archivo_main));
ext_main = '.csv';
sep_main = ';';


#################################### FUNCIONES ####################################


def abrir_archivo_fuente(dir_arch, nom_arch, ext='.csv', sep=';'):
    # Abre el archivo con los sitios de union, genes y fuentes
    # Devuelve matriz con toda la informacion

    # Inicializo la matriz que se devuelve
    M_out = [];

    # Abro el archivo
    with open(dir_arch + ext, 'r') as F:
        print('Archivo ' + str(nom_arch) + ext + ' abierto')
        # Leo cada linea del archivo
        for curr_line in F:
            # Separo cada linea en una lista con split(sep)
            L_curr_line = curr_line.rstrip().split(sep);
            # Guardo la lista como fila en M_out
            M_out.append(L_curr_line[:]);
    return M_out


def agregar_dist_pos0(M_csv, genome, gene_id_pos=1):
    # Recorre M_csv y agrega la distancia a pos0 del gen al final de cada fila
    # Gen determinado por M_csv[i][gene_id_pos] y buscado en genome

    # Recorro M_csv
    for i in range(len(M_csv)):
        curr_line = M_csv[i];
        # Agarro el gen de ensembl con gene_id
        curr_gene = genome.gene_by_id(curr_line[gene_id_pos]);
        # Defino pos0 y forward para curr_gene
        pos0, forward = buscar_pos0(curr_gene);
        # Defino la distancia en base a pos0, forward y pos_ini/pos_end del sitio
        dist_sitio = calcular_distancia_sitio(pos0, forward, int(curr_line[4]), int(curr_line[5]));
        # Agrego dist_sitio a curr_line
        M_csv[i].append(int(dist_sitio));
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


def calcular_distancia_sitio(pos0, forward, pos_ini_sitio, pos_end_sitio):
    # Funcion que calcula la distancia de un sitio de union a pos0 de un gen

    # Tres casos primero: sitio sobre 0, sitio antes/despues del 0
    if pos0 >= pos_ini_sitio and pos0 <= pos_end_sitio:
        # Si el sitio esta sobre 0, se devuelve 0
        r = 0;
    elif pos0 < pos_ini_sitio:
        # Si el sitio esta "despues" de pos0, se usa pos_ini
        if forward:
            r = pos_ini_sitio - pos0;
        else:
            r = pos0 - pos_ini_sitio;
    elif pos0 > pos_end_sitio:
        # Si el sitio esta "antes" de pos0, se usa pos_end
        if forward:
            r = pos_end_sitio - pos0;
        else:
            r = pos0 - pos_end_sitio;
    return r


def extraer_sitios_genoma(dir_arch, nom_arch, gene_id, gene_id_pos=1, ext='.csv', sep=';'):
    # Abre el archivo con los sitios de union, genes y fuentes
    # Devuelve matriz solo con los sitios ubicados cerca de genes
    # gene_id='ENSMUSG' para raton y gene_id='ENSG' para humano

    # Inicializo la matriz que se devuelve
    M_out = [];
    
    # Extraigo la matriz del archivo
    M_csv = abrir_archivo_fuente(dir_arch, nom_arch, ext=ext, sep=sep);

    # Recorro M_csv
    for i in range(len(M_csv)):
        curr_line = M_csv[i];
        # Reviso si curr_line[gene_id_pos] empieza con gene_id
        if curr_line[gene_id_pos][:len(gene_id)] == gene_id:
            M_out.append(curr_line[:]);
        elif len(curr_line[gene_id_pos])>0:
            print('ERROR. Gen encontrado no identificado: ' + str(curr_line[gene_id_pos]))
    return M_out


def pipeline_histogramas(dir_arch, nom_arch, genome_name, gene_id_pos=1, ext='.csv', sep=';'):
    # Funcion que corre todas las funciones en orden para generar los distintos histogramas

    ### Display
    print('Inicializando pipeline de creacion de histogramas.')
    
    # Inicializo la matriz que se devuelve
    M_out = [];

    # Defino gene_id en base a genome_name
    if genome_name == 'mm9':
        gene_id = 'ENSMUSG';
        genome_ensembl = mm9;
    elif genome_name == 'hg19':
        gene_id = 'ENSG';
        genome_ensembl = hg19;
    else:
        print('ERROR. genome_name no interpretado.')
        return M_out
    # Extraigo M_csv de dir_arch con los sitios ubicados cerca de genes
    M_csv = extraer_sitios_genoma(dir_arch, nom_arch, gene_id, gene_id_pos=gene_id_pos, ext=ext, sep=sep);

    # Agrego distancia a pos0 del sitio de union registrado
    M_csv = agregar_dist_pos0(M_csv, genome_ensembl, gene_id_pos=gene_id_pos);
    
    ### Display
    print('Datos listos para empezar generacion de histogramas.')

    M_out = M_csv;
    ########################## HACER ESTO ##########################
    ### FALTA:
    ## Terminar de procesar datos (ver mas arriba)
    ## Hacer cada histograma (desguasar este punto)
    return M_out


def _main():
    # Funcion para probar funciones en ejecucion del archivo

    # Inicializo M_out
    M_out = [];

    # Corro el pipeline
    M_out = pipeline_histogramas(dir_archivo_main, nombre_archivo_main, genoma_main,
                                 gene_id_pos=1, ext=ext_main, sep=sep_main);
    return M_out


#################################### RESULTADOS ###################################

output_dump = [];

if __name__=='__main__':
    output_dump.append(_main());

