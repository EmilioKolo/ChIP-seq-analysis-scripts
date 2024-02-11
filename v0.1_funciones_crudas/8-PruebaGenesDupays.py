
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
from referencias_funciones import ConsultaSecuencia, IDchr, abrir_matriz_pbm

mm9 = EnsemblRelease(67, species='mouse');
hg19 = EnsemblRelease(102, species='human');
GRCm38 = EnsemblRelease(102, species='mouse');


'''
Pruebas con los 28 genes agarrados al azar en el paper Dupays (L_genes_total y L_genes_id)

Funcion para buscar Ensembl IDs de listas de nombres de genes: L_genes_to_id()

Funcion para buscar distancias a picos mas cercanos en varios csv: analisis_peaks_L_genes()
'''

#################################### VARIABLES ####################################


lista_dupays = ['Dupays2015_s1e14', 'Dupays2015_s4e14'];

L_genes_total = ['3110057O12Rik', 'Actc1', 'Cacna2d1', 'Frem1', 'Gja5', 'H2afx',
                   'Homer1', 'Lrtm1', 'Lyz1', 'Mtbp', 'Mtss1', 'Myl7', 'Pde4b', 'Pde4d',
                   'Pgam2', 'Phka2', 'Ppm1k', 'Rcan1', 'Smpx', 'Tcap', 'Tmem97', 'Tnni2',
                   'Txnip', 'Wbp5', 'Ppp1r17', 'Hist2h4', 'Ppargc1a', 'Tecrl'];

L_genes_dificiles = ['Gsbs', 'Hist2H4', 'Pparga1a', 'Srd5a2l2'];
L_genes_dificiles = ['Ppp1r17', 'Hist2h4', 'Ppargc1a', 'Tecrl'];
# Gsbs == Ppp1r17
# Hist2H4 == Hist2h4
# Pparga1a == Ppargc1a ??
# Srd5a2l2 == Tecrl

L_genes_id = ['ENSMUSG00000037818', 'ENSMUSG00000068614', 'ENSMUSG00000040118', 'ENSMUSG00000059049',
              'ENSMUSG00000057123', 'ENSMUSG00000049932', 'ENSMUSG00000007617', 'ENSMUSG00000045776',
              'ENSMUSG00000069515', 'ENSMUSG00000022369', 'ENSMUSG00000022353', 'ENSMUSG00000020469',
              'ENSMUSG00000028525', 'ENSMUSG00000021699', 'ENSMUSG00000074661', 'ENSMUSG00000020475',
              'ENSMUSG00000031295', 'ENSMUSG00000037826', 'ENSMUSG00000022951', 'ENSMUSG00000041476',
              'ENSMUSG00000007877', 'ENSMUSG00000037278', 'ENSMUSG00000031097', 'ENSMUSG00000038393',
              'ENSMUSG00000042712', 'ENSMUSG00000002930', 'ENSMUSG00000091405', 'ENSMUSG00000029167',
              'ENSMUSG00000049537', 'ENSMUSG00000022803'];

L_genes_id_duplicados = ['ENSMUSG00000021699', 'ENSMUSG00000074661'];
L_genes_id_dificiles = ['ENSMUSG00000002930', 'ENSMUSG00000091405', 'ENSMUSG00000029167', 'ENSMUSG00000049537']

extension = '.csv';
csv_sep = ';';
col_csv = [1,2,3];

path_csv = '.\\csv con secuencias\\';
curr_path = os.path.dirname(__file__);

# Variables main()

L_archivos_main = lista_dupays;
L_genes_nombres_main = L_genes_total;
genoma_usado = mm9;
L_genes_id_main = L_genes_id;

L_csv_main = [1,2,3];

# Variables testing

nombre_csv_test = lista_dupays[1];
dir_csv_test = os.path.join(curr_path, path_csv + nombre_csv_test + '.csv');
L_csv_test = [1,2,3];


#################################### FUNCIONES ####################################


def abrir_csv_extraer_columnas(path_csv, nombre_csv, L_cols, sep_arch=';'):
    # Funcion que abre archivo .csv y devuelve una matriz con columnas seleccionadas

    # Inicializo la matriz que se devuelve
    M_out = [];
    # Abro el archivo .csv
    with open(path_csv, 'r') as F:
        print('Archivo ' + str(nombre_csv) + ' abierto.')
        # Recorro cada linea del archivo
        for curr_line in F:
            # Paso cada linea a formato de lista
            L = curr_line.rstrip().split(sep_arch);
            # Creo una lista para registrar las columnas seleccionadas
            L_out = [];
            # Si hay numeros en L_col, paso solo esas columnas
            if len(L_cols) > 0:
                for col in L_cols:
                        L_out.append(L[col]);
            # Si L_col es una lista vacia, se devuelven todas las columnas
            else:
                L_out = L[:];
            # Agrego L_out a M_out
            M_out.append(L_out[:]);
    return M_out


def analisis_peaks_L_genes(L_nombres_arch, path_arch, L_col_csv, L_genes, genoma, max_dist=100000000, verbose=False):
    # Recibe una lista de archivos .csv con peaks y una lista de genes
    # Busca 

    # Inicializo la matriz que se devuelve
    L_dist_peak = [];

    # Prueba fuerza bruta: Listas de matrices para archivos csv y matriz de distancias a peaks
    M_M_csv = [];

    # Abro cada .csv y los guardo en M_M_csv
    for nombre_curr_csv in L_nombres_arch:
        dir_curr_csv = os.path.join(curr_path, path_arch + nombre_curr_csv + '.csv');
        curr_M_csv = abrir_csv_extraer_columnas(dir_curr_csv, nombre_curr_csv, L_col_csv);
        M_M_csv.append(curr_M_csv[:]);

    # Reviso cada gen de a uno a la vez
    for i in range(len(L_genes)):
        # Anoto el Ensembl id
        curr_gene = L_genes[i];
        # Agarro el elemento gen de Ensembl
        curr_ensembl_gene = genoma.gene_by_id(curr_gene);
        print('Iniciando analisis de gen ' + str(curr_ensembl_gene.gene_name))
        # Inicializo la distancia como el maximo aceptado y registro si se encuentra una menor
        curr_dist = max_dist;
        dist_encontrada = False;
        # Recorro cada matriz de cada .csv
        for j in range(len(M_M_csv)):
            if verbose:
                print('Iniciando analisis de .csv numero ' + str(j+1))
            curr_M_csv = M_M_csv[j];
            # Busco la distancia al peak mas cercano para cada .csv
            d_peak_cercano = peak_cercanos_un_gen(curr_ensembl_gene, curr_M_csv, return_min=True, verbose=0);
            # Si la distancia encontrada es menor que la registrada, se registra y se marca que se encontro por lo menos una distancia
            if int(d_peak_cercano) < curr_dist:
                curr_dist = int(d_peak_cercano);
                dist_encontrada = True;
        # Registro el valor de distancia minima encontrada
        if dist_encontrada:
            L_dist_peak.append(curr_dist);
        # Si no se encontro ninguna distancia, se registra -1
        else:
            print('No se encontro peak a menos de ' + str(max_dist) + ' pares de bases.')
            L_dist_peak.append(-1);
    # Display
    if verbose:
        print()
        print('DISTANCIAS MINIMAS DE LOS GENES AL PEAK MAS CERCANO EN ARCHIVOS:')
        for f in L_nombres_arch:
            print(f + '.csv')
        print()
        for k in range(len(L_genes)):
            curr_gene_display = genoma.gene_by_id(L_genes[k]);
            print('Distancia minima al gen ' + str(curr_gene_display.gene_name) + ': ' + str(L_dist_peak[k])) 
    return L_dist_peak


def L_genes_to_id(L_genes, genoma, verbose=False):
    # Funcion que recibe una lista de nombres de genes y devuelve Ensembl id
    # Usa genoma.gene_ids_of_gene_name()

    # Inicializo matriz que se devuelve
    M_out = [];

    # Recorro cada gen de L_genes
    for i in range(len(L_genes)):
        curr_gene_name = L_genes[i];
        # Busco los Ensembl IDs asociados al nombre en L_genes
        genes_encontrados = genoma.gene_ids_of_gene_name(curr_gene_name);
        # Display
        if verbose:
            print('Gen: ' + str(curr_gene_name))
            print(genes_encontrados)
        # Recorro cada uno de los Ensembl IDs encontrados
        for j in range(len(genes_encontrados)):
            # Display
            if verbose:
                print(genoma.gene_by_id(genes_encontrados[j]))
            # Registro cada Ensembl ID encontrado en M_out
            M_out.append(genes_encontrados[j]);
    return M_out


def main():
    # Funcion para probar funciones en ejecucion del archivo

    M_out = [];

    # Funcion para obtener Ensembl IDs de la lista de genes
    #M_out = L_genes_to_id(L_genes_nombres_main, genoma_usado, verbose=True);

    # Funcion para buscar columnas dentro de un .csv
    #M_csv = abrir_csv_extraer_columnas(dir_csv_test, nombre_csv_test, L_csv_test);

    # Funcion para devolver los peaks mas cercanos a cada gen dentro de un .csv
    #M_out = peak_mas_cercano_L_genes(M_csv, L_genes_id_main, genoma_usado);

    # Funcion para buscar la distancia minima de varios genes a peaks de varios .csv
    M_out = analisis_peaks_L_genes(L_archivos_main, path_csv, L_csv_main, L_genes_id_main, genoma_usado, verbose=True);

    return M_out


def peak_cercanos_un_gen(gen_busq, M_csv, return_min=False, verbose=1):
    # Funcion que devuelve los peaks mas cercanos a gen_busq en M_csv
    # Devuelve la distancia al peak mas cercano a cada lado del gen
    # Devuelve la distancia al mismo peak dos veces si esta por encima o adentro del gen (dist=0)
    # M_csv con filas en formato [contig, pos_ini, pos_end], todo como str()
    # gen_busq es un elemento gen de Ensembl

    if verbose>0:
        print('Inicializando busqueda de peaks cerca de gen ' + str(gen_busq.gene_name))

    # Inicializo la lista/matriz que se devuelve
    L_out = [];

    # Defino contig, pos_ini y pos_end de gen_busq
    contig_gen = str(gen_busq.contig);
    pos_ini = str(gen_busq.start);
    pos_end = str(gen_busq.end);

    # Defino booleanos y peaks mas cercanos a cada lado / dentro del gen
    dist_peak_neg = -1; # Numero positivo que indica distancia a peak mas cercano en 5' del gen
    dist_peak_pos = -1; # Numero positivo que indica distancia a peak mas cercano en 3' del gen
    peak_neg = []; # Peak mas cercano a 5' del gen, inicializa vacio
    peak_pos = []; # Peak mas cercano a 3' del gen, inicializa vacio

    # Contador de chequeo para contig
    cont_contig = 0;
    
    # Recorro cada peak de M_csv
    for i in range(len(M_csv)):
        curr_peak = M_csv[i];
        # Trabajo solo con peaks en el mismo contig que el gen
        if str(curr_peak[0]) == str(contig_gen):
            cont_contig = cont_contig + 1;
            # Defino si el peak esta a 5' o 3' del gen
            # Si esta a 5', curr_peak[2] es menor que pos_ini del gen
            if (int(curr_peak[2]) < int(pos_ini)):
                # Solo registro el peak si no hay otro a distancia 0 o menor
                if (dist_peak_neg > 0) or (len(peak_neg) == 0):
                    dist_curr_peak = int(pos_ini) - int(curr_peak[2]);
                    # Veo si ya hay un peak 5' registrado
                    if len(peak_neg) == 0:
                        # Si no hay ningun peak registrado, registro el peak y la distancia
                        peak_neg = curr_peak[:];
                        dist_peak_neg = max(int(dist_curr_peak), 0);
                    # Si ya hay un peak registrado, reviso que dist_curr_peak sea menor
                    elif dist_curr_peak < dist_peak_neg:
                        peak_neg = curr_peak[:];
                        dist_peak_neg = max(int(dist_curr_peak), 0);
            # Si esta a 3', curr_peak[1] es mayor que pos_end del gen
            elif (int(curr_peak[1]) > int(pos_end)):
                # Solo registro el peak si no hay otro a distancia 0 o menor
                if (dist_peak_pos > 0) or (len(peak_pos) == 0):
                    dist_curr_peak = int(curr_peak[1]) - int(pos_end);
                    # Veo si ya hay un peak 3' registrado
                    if len(peak_pos) == 0:
                        # Si no hay ningun peak registrado, registro el peak y la distancia
                        peak_pos = curr_peak[:];
                        dist_peak_pos = max(int(dist_curr_peak), 0);
                    # Si ya hay un peak registrado, reviso que dist_curr_peak sea menor
                    elif dist_curr_peak < dist_peak_pos:
                        peak_pos = curr_peak[:];
                        dist_peak_neg = max(int(dist_curr_peak), 0);
            # El peak esta dentro del gen, por lo menos una de las distancias es 0
            # Si curr_peak[1] es menor que pos_ini del gen y el peak no pasa sobre todo el gen, es 5'
            elif (int(curr_peak[1]) < int(pos_ini)) and (int(curr_peak[2]) < int(pos_end)):
                dist_curr_peak = 0;
                # Si el peak anterior estaba a distancia mayor a 0, lo sobreescribo
                if (dist_peak_neg > 0) or (len(peak_neg) == 0):
                    peak_neg = curr_peak[:];
                    dist_peak_neg = dist_curr_peak;
                # Si hay mas de un peak a distancia 0, agrego otro mas como lista a peak_neg
                else:
                    if verbose>2:
                        print('Mas de un peak dentro del gen, revisar.')
                    peak_neg.append(curr_peak[:]);
            # Si curr_peak[2] es mayor que pos_end del gen y el peak no pasa sobre todo el gen, es 3'
            elif (int(curr_peak[2]) > int(pos_end)) and (int(curr_peak[1]) > int(pos_ini)):
                dist_curr_peak = 0;
                # Si el peak anterior estaba a distancia mayor a 0, lo sobreescribo
                if (dist_peak_pos > 0) or (len(peak_pos) == 0):
                    peak_pos = curr_peak[:];
                    dist_peak_pos = dist_curr_peak;
                # Si hay mas de un peak a distancia 0, agrego otro mas como lista a peak_pos
                else:
                    if verbose>2:
                        print('Mas de un peak dentro del gen, revisar.')
                    peak_pos.append(curr_peak[:]);
            # Si todo lo demas falla, el peak esta entero dentro del gen o el gen esta entero dentro del peak
            else:
                dist_curr_peak = 0;
                # Si el peak a 5' anterior estaba a distancia mayor a 0, lo sobreescribo
                if (dist_peak_neg > 0) or (len(peak_neg) == 0):
                    peak_neg = curr_peak[:];
                    dist_peak_neg = dist_curr_peak;
                # Si hay mas de un peak a distancia 0, agrego otro mas como lista a peak_neg
                else:
                    if verbose>2:
                        print('Mas de un peak dentro del gen, revisar.')
                    peak_neg.append(curr_peak[:]);
                # Si el peak 3' anterior estaba a distancia mayor a 0, lo sobreescribo
                if (dist_peak_pos > 0) or (len(peak_pos) == 0):
                    peak_pos = curr_peak[:];
                    dist_peak_pos = dist_curr_peak;
                # Si hay mas de un peak a distancia 0, agrego otro mas como lista a peak_pos
                else:
                    if verbose>2:
                        print('Mas de un peak dentro del gen, revisar.')
                    peak_pos.append(curr_peak[:]);
    # Reviso que aunque sea un contig se haya encontrado
    if cont_contig == 0:
        print('Error con el contig, no se encontro ningun peak en contig ' + str(contig_gen))
    # Devuelvo ambos peaks aunque sean iguales, se chequea despues
    if verbose>1:
        print('Gen: ' + str(contig_gen) + ' : ' + str(pos_ini) + ' - ' + str(pos_end))
        print(str(peak_neg) + ' d=' + str(dist_peak_neg))
        print(str(peak_pos) + ' d=' + str(dist_peak_pos))
    L_out = (int(dist_peak_neg), int(dist_peak_pos));
    # return_min hace que devuelva un solo numero
    if return_min:
        L_out = min(int(dist_peak_neg), int(dist_peak_pos))
    return L_out


def peak_mas_cercano_L_genes(M_csv, L_genes, genoma):
    # Devuelve el peak mas cercano de M_csv a cada gen de L_genes, a cada lado del gen
    # M_csv con filas en formato [contig, pos_ini, pos_end], todo como str()
    # L_genes es una lista de Ensembl IDs de genes

    # Inicializo la matriz que se devuelve
    M_out = [];
    # Recorro cada uno de los genes
    for curr_gene in L_genes:
        # Corro una funcion que devuelve los peaks mas cercanos a cada gen de L_gene
        L_peaks_cercanos = peak_cercanos_un_gen(genoma.gene_by_id(curr_gene), M_csv);
        # Agrego L_peaks_cercanos a M_out
        M_out.append(L_peaks_cercanos[:]);
    return M_out


#################################### RESULTADOS ###################################

output_dump = '';

if __name__=='__main__':
    output_dump = main();

