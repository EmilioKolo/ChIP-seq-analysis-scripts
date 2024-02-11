
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
Prepara secuencia alrededor de sitios de union para buscar otros sitios en MEME

Abre .csv con sitios de union en peaks de ChIP-seq
    * extraer_SUP()

Procesa los sitios y agrega rangos
    * preparar_rangos_FIMO()

Busca las secuencias delimitadas por los sitios y las pasa a formato .fasta para FIMO
    * guardar_fasta()

Todo junto corre en crear_fasta_SUP() con los resultados de extraer_SUP()

Corre para todos los rangos en L_rangos

FALTA:
- Funcion para seleccionar secuencias de union dentro del .fasta
    * Permite buscar particularidades de cada secuencia de union

? Usar seleccionar_sitios() para agarrar solo los sitios en Dupays2015
    * Resuelto (no me acuerdo si use seleccionar_sitios())
'''

#################################### VARIABLES ####################################

# Variables de archivos fuente

sitios_union_peaks = 'SitiosUnionFinal';
ext_SUP = '.csv';
path_SUP = '.\\sitios de union FINAL\\';

# Variables para output

L_rangos = [100, 200];
out_fasta = 'SecuenciasSitios';

# Variables generales

curr_path = os.path.dirname(__file__);

# Variables main

genome_main = 'mm9';
ext_main = ext_SUP;
nom_main = sitios_union_peaks + '_' + genome_main;
path_main = path_SUP;
dir_main = os.path.join(curr_path, path_main + nom_main);
L_rango_main = L_rangos;
n_prt_main = 200;
out_fasta_main = out_fasta;

# Variables de testeo
n_prt_main_test = 1000;
rango_test_main = 500;


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


def buscar_secuencias_FIMO(M_SU_FIMO, verbose_range=100, tiempo_sleep=300):
    # Busca las secuencias de las regiones delimitadas por M_SU_FIMO y las devuelve en formato .fasta para FIMO
    # Devuelve lista de strings, cada uno correspondiendo a una entrada de .fasta

    # Inicializo la matriz que se devuelve
    M_out = [];
    # Variable de display
    len_M = len(M_SU_FIMO);
    # Recorro M_SU_FIMO
    for i in range(len(M_SU_FIMO)):
        curr_SU = M_SU_FIMO[i];
        # Saco chromosome ID con la funcion IDchr()
        chr_id = IDchr(curr_SU[0][3:], genome=curr_SU[-1]);
        # Busco la secuencia correspondiente al rango en curr_SU
        curr_seq = str(ConsultaSecuencia(chr_id, curr_SU[1], curr_SU[2], strand=1, sleep_time=tiempo_sleep));
        # Escribo el texto en formato .fasta para FIMO
        r = '>ID' + str(i+1) + ',' + curr_SU[3];
        r = r + '\n';
        r = r + curr_seq;
        r = r + '\n';
        # Agrego el texto a M_out
        M_out.append(str(r));
        # Display
        if (i+1)%verbose_range == 0:
            print('Progreso de busqueda de secuencias: ' + str(i+1) + ' / ' + str(len_M))
            print(r)
    return M_out


def crear_fasta_SUP(M_SUP, rango_sitio, nom_out, genome, verbose_range=100):
    # Funcion que aglomera preparar_rangos_FIMO(), buscar_secuencias_FIMO() y guardar_fasta()

    # Inicializo la variable que se devuelve
    M_out = [];
    nombre_archivo_out = str(nom_out) + '_rango' + str(rango_sitio);

    # Paso M_SUP a formato para buscar secuencias y agrego el rango
    M_SU_FIMO = preparar_rangos_FIMO(M_SUP, rango_sitio, genome);
    print('M_SU_FIMO lista para buscar secuencias.')
    # Busco las secuencias y las devuelvo en formato FIMO
    M_str_FIMO = buscar_secuencias_FIMO(M_SU_FIMO, verbose_range=verbose_range);
    print('Secuencias de M_SU_FIMO buscadas y guardadas en M_str_FIMO.')
    # Guardo las secuencias en un archivo .fasta y devuelvo M_str_FIMO
    M_out = guardar_fasta(M_str_FIMO, nombre_archivo_out);
    print('Secuencias de M_str_FIMO guardadas en archivo "' + nombre_archivo_out + '.fasta"')
    print()
    return M_out


def definir_SU_humano_raton(curr_SU):
    # Funcion para determinar si un sitio de union pertenece a humano o raton
    # Pensada para sitios extraidos con extraer_SUP()
    # Hay que actualizarla si cambian los datos

    # Inicializo la variable que se devuelve
    r = '';
    # Reviso que el sitio de union este en Anderson
    if curr_SU[5] == 'X':
        # Si esta en mas lugares que Anderson, es un error
        if 'X' in curr_SU[6:]:
            r = 'ERROR Anderson + Raton';
        else:
            r = 'hg19';
    # Reviso que el sitio de union este en alguna de las demas posiciones
    elif 'X' in curr_SU[6:]:
        r = 'mm9';
    # Si el sitio no esta en ninguna fuente, tiro error
    else:
        r = 'ERROR sin X';
    return r


def extraer_SUP(dir_arch, nom_arch, sep_arch=';'):
    # Usa abrir_resultados_csv() para seleccionar sitios de union de peaks en dir_arch
    # dir_arch es la tabla final de sitios de union

    # Inicializo la matriz que se devuelve
    M_out = [];

    # Extraigo todo del .csv en dir_arch
    M_csv = abrir_resultados_csv(dir_arch, nom_arch, sep_arch=sep_arch);

    # Recorro M_csv
    for i in range(len(M_csv)):
        curr_SU = M_csv[i];
        # Reviso si el sitio de union es sacado de un peak de ChIP
        chip_peak = False;
        # Solo me quedo con sitios de union que tengan X (marcando que son peaks de ChIP-seq)
        for j in curr_SU:
            if j == 'X':
                chip_peak = True;
        # Si es sacado de peak de ChIP-seq, lo guardo en M_out
        if chip_peak:
            M_out.append(curr_SU[:]);
    return M_out


def guardar_fasta(M_fasta, nom_out):
    # Crea un archivo .fasta de nombre nom_out y guarda los contenidos de M_fasta
    # M_fasta tiene que ser una lista de strings

    # Creo el archivo nom_out y borro si existe uno anterior
    with open(nom_out + '.fasta', 'w') as F_out:
        print('Archivo "' + nom_out + '.fasta" creado')
    # Vuelvo a abrir el archivo en formato append para guardar la info en M_fasta
    with open(nom_out + '.fasta', 'a') as F_out:
        # Escribo cada una de las lineas de M_fasta en F_out
        ### SI QUIERO FILTRAR POR TIPO DE ARCHIVO SE HACE ACA
        for r in M_fasta:
            F_out.write(r);
    return M_fasta


def preparar_rangos_FIMO(M_SUP, rango_seq, curr_SU_origin):
    # Prepara los sitios de union sacados con extraer_SUP() para buscar secuencias
    # Devuelve chr_n, pos_ini, pos_end, seq_SU, gen cercano (si hay) y genoma
    ### Genoma dado en curr_SU_origin
    # Defino genoma a mano (primer X es hg19 y el resto son mm9)

    # Inicializo la matriz que se devuelve
    M_out = [];
    # Recorro M_SUP
    for i in range(len(M_SUP)):
        curr_SU = M_SUP[i];
        # Inicializo el sitio de union a devolver
        curr_SU_out = [];
        # Agrego la info que esta disponible de una
        curr_SU_out.append(str(curr_SU[3]));
        curr_SU_out.append(str(int(curr_SU[4])-rango_seq));
        curr_SU_out.append(str(int(curr_SU[5])+rango_seq));
        curr_SU_out.append(str(curr_SU[0]));
        curr_SU_out.append(str(curr_SU[1]));
        ### VIEJO
        # Defino si el sitio de union corresponde a humano o raton
        #curr_SU_origin = definir_SU_humano_raton(curr_SU);
        ###
        if curr_SU_origin == 'mm9':
            curr_SU_out.append('mm9');
        elif curr_SU_origin == 'hg19':
            curr_SU_out.append('hg19');
        else:
            print('ERROR: Sitio de union con origen mal determinado.')
            print(curr_SU)
            print('Origen: ' + str(curr_SU_origin))
            curr_SU_out.append(str(curr_SU_origin));
        # Agrego curr_SU_out a M_out
        M_out.append(curr_SU_out[:]);
    return M_out


def seleccionar_sitios(M_SU, id_busq=0, str_busq='X', hacer_seleccion=False):
    # Selecciona, dentro de M_SU, los sitios de union con str_busq en posicion id_busq
    # Si hacer_seleccion es False, se devuelve M_SU como se recibio
    if hacer_seleccion:
        M_out = [];
        # Revisar cada uno de los sitios de union
        for i in range(len(M_SU)):
            curr_SU = M_SU[i];
            # Si el criterio de seleccion se cumple, se agrega curr_SU a M_out
            if curr_SU[id_busq] == str_busq:
                M_out.append(curr_SU[:]);
    else:
        M_out = M_SU;
    return M_out


def seleccionar_SUP(M_SUP, id_selec, mark_selec='X'):
    # Selecciona los sitios de M_SUP con mark_selec en la posicion id_selec

    # Inicializo la matriz que se devuelve
    M_out = [];

    # Recorro M_SUP
    for i in range(len(M_SUP)):
        curr_SU = M_SUP[i];
        # Selecciono los SU con mark_selec en posicion id_selec
        if curr_SU[id_selec] == mark_selec:
            M_out.append(curr_SU[:]);
    return M_out


def _main():
    # Funcion para usar funciones en ejecucion del archivo

    # Inicializo la variable que se devuelve
    M_out = [];

    # Extraigo los sitios de union de ChIP-seq
    M_SUP = extraer_SUP(dir_main, nom_main, sep_arch=';');
    print('Extraccion SUP terminada.')
    #print(_range_print(M_SUP, n_prt_main_test, shuffle_M=True));

    # Corro crear_fasta_SUP() para todos los rangos seleccionados
    for curr_range in L_rango_main:
        print('Iniciando formateo de .fasta para rango=' + str(curr_range))
        M_str_FIMO = crear_fasta_SUP(M_SUP, curr_range, out_fasta_main, verbose_range=n_prt_main);

    ### OLD
    #M_SU_FIMO = preparar_rangos_FIMO(M_SUP, rango_test_main);
    #print('M_SU_FIMO lista para buscar secuencias.')
    #print(_range_print(M_SU_FIMO, n_prt_main_test, shuffle_M=True))
    #M_str_FIMO = buscar_secuencias_FIMO(M_SU_FIMO, verbose_range=n_prt_main);
    #print('Secuencias de M_SU_FIMO buscadas y guardadas en M_str_FIMO.')
    #print(_range_print(M_str_FIMO, n_prt_main_test, shuffle_M=True))
    #M_out = guardar_fasta(M_str_FIMO, out_fasta_main + 'rango' + str(rango_test_main));
    ###

    M_out = M_str_FIMO;
    
    return M_out


def _main_Dupays():
    # Corro el pipeline de _main() pero solo agarro los resultados de Dupays
    
    # Inicializo la variable que se devuelve
    M_out = [];

    # Extraigo los sitios de union de ChIP-seq
    M_SUP = extraer_SUP(dir_main, nom_main, sep_arch=';');
    print('Extraccion SUP terminada.')
    print()
    print('Seleccionando sitios de union en peaks de Dupays 2015')
    id_dupays = 6; ### REVISAR ESTO ANTES DE CORRER
    M_SUP = seleccionar_SUP(M_SUP, id_dupays, mark_selec='X');
    print('Seleccion de sitios Dupays terminada.')
    print()

    # Corro crear_fasta_SUP() para todos los rangos seleccionados
    for curr_range in L_rango_main:
        print('Iniciando formateo de .fasta para rango=' + str(curr_range))
        M_str_FIMO = crear_fasta_SUP(M_SUP, curr_range, out_fasta_main + '_Dupays', genome_main, verbose_range=n_prt_main);
    return M_out


def _main_test():
    # Funcion para probar funciones sueltas

    # Inicializo la variable que se devuelve
    M_out = [];

    return M_out


def _range_print(M_prt, prt_range, shuffle_M=False):
    # Hace print de elementos de la matriz M_prt a intervalos de prt_range
    # shuffle_M permite ir variando el print
    # Devuelve el largo de M_prt
    if shuffle_M:
        shuffle(M_prt)
    for i in range(len(M_prt)):
        if i%prt_range == 0:
            print(M_prt[i])
    return len(M_prt)

#################################### RESULTADOS ###################################

output_dump = [];

if __name__=='__main__':
    #output_dump.append(_main());
    #output_dump.append(_main_test());
    output_dump.append(_main_Dupays());


