
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
from referencias_funciones import ConsultaSecuencia, IDchr 

mm9 = EnsemblRelease(67, species='mouse');
hg19 = EnsemblRelease(102, species='human');
GRCm38 = EnsemblRelease(102, species='mouse');

'''
Genera distintas figuras en base a los resultados de 9-TablaFinalSitiosUnion.py

Abre los resultados con abrir_resultados_csv()

FALTA:
- Definir que figuras hacer
    - Venn?
    - Histogramas?
    - Otras?
'''

#################################### VARIABLES ####################################

# Variables de archivos fuente

sitios_union_total = 'SitiosUnionFinal';
ext_SUT = '.csv';
path_SUT = '.\\sitios de union FINAL\\';

# Variables para los graficos



# Variables generales

curr_path = os.path.dirname(__file__);

# Variables main

ext_main = ext_SUT;
nom_main = sitios_union_total;
dir_main = os.path.join(curr_path, path_SUT + nom_main);


#################################### FUNCIONES ####################################


def abrir_resultados_csv(dir_arch, nom_arch, sep_arch=';'):
    # Abre el archivo .csv en dir_arch y devuelve todas las columnas como matriz

    # Inicializo la matriz a devolver
    M_csv = [];

    # Abro el archivo
    with open(dir_arch + '.csv', 'r') as F:
        print('Archivo ' + str(nom_arch) + '.csv abierto.')
        for curr_line in F:
            L_curr_line = curr_line.rstrip().split(sep_arch);
            M_csv.append(L_curr_line[:])
    return M_csv


def _main():
    # Funcion para usar funciones en ejecucion del archivo

    # Inicializo la variable que se devuelve
    M_out = [];

    #
    M_out = abrir_resultados_csv(dir_main, nom_main, sep_arch=';')

    return M_out


def _main_test():
    # Funcion para probar funciones sueltas

    # Inicializo la variable que se devuelve
    M_out = [];

    return M_out


#################################### RESULTADOS ###################################

output_dump = '';

if __name__=='__main__':
    output_dump = _main();
    #output_dump = _main_test();



