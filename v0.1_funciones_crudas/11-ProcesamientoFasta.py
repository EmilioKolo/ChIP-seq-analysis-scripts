
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
Abre los archivos .fasta creados por 10-SitiosUnionChIP.py y los modifica

Funcion leer_fasta() para abrir los archivos y devolverlos en formato de matriz

Funcion eliminar_sitios_M() para eliminar el sitio de union en el medio
    * Funcion eliminar_sitio() determina por que se reemplaza el sitio
        * Inicialmente se reemplaza por str('N' * largo_sitio)

Funcion guardar_M_fasta() guarda la matriz en un archivo fasta nuevo

FALTA:
- Funcion para crear .fasta en base a .csv con picos ChIP-seq
    * BUSCAR EN OLD
'''

#################################### VARIABLES ####################################

nombre_fasta100 = 'SecuenciasSitios_Dupays_rango100';
nombre_fasta200 = 'SecuenciasSitios_Dupays_rango200';

path_fasta = '.\\fasta para MEME\\';

curr_path = os.path.dirname(__file__);

dir_fasta100 = os.path.join(curr_path, path_fasta + nombre_fasta100);
dir_fasta200 = os.path.join(curr_path, path_fasta + nombre_fasta200);

# Variables main()

nombre_main = nombre_fasta100;
dir_main = dir_fasta100;
largo_sitio_main = 7;
nombre_out_main = nombre_main + '_sin_sitio';


#################################### FUNCIONES ####################################


def eliminar_sitio(str_fasta, len_sitio):
    # Revisa str_fasta y busca la region de largo len_sitio en el medio
    # Cambia la region buscada por seq_reemplazo (definida en funcion)
    # Probablemente no anda si len(str_fasta)-len_sitio no es par

    # Reviso si len(str_fasta)-len_sitio es par
    if (len(str_fasta)-len_sitio)%2 != 0:
        print('CUIDADO: Secuencia str_fasta y len_sitio no compatibles')
    # Determino posicion inicial y final del sitio
    pos_ini = int((len(str_fasta) - len_sitio)/2);
    pos_end = pos_ini + len_sitio;

    # Determino la secuencia que se va a usar para reemplazar al sitio
    seq_reemplazo = 'N'*len_sitio;
    
    # En base a eso, construyo la secuencia que se devuelve
    # Agrego la secuencia antes del sitio
    r = str_fasta[:pos_ini];
    # Agrego la secuencia que reemplaza al sitio
    r = r + str(seq_reemplazo);
    # Agrego la secuencia despues del sitio
    r = r + str_fasta[pos_end:];
    return r


def eliminar_sitios_M(M_fasta, len_sitio):
    # Recorre M_fasta y elimina un sitio de largo len_sitio ubicado en el medio de la secuencia
    # Devuelve una matriz similar con el sitio eliminado (o reemplazado)
    # Ver eliminar_sitio() para determinar que se hace con el sitio eliminado

    # Inicializo la matriz que se devuelve
    M_out = [];

    # Recorro M_fasta
    for i in range(len(M_fasta)):
        curr_fasta = M_fasta[i];
        # Inicializo el fasta sin sitio
        curr_fasta_out = [];
        # Agrego curr_fasta[0] a curr_fasta_out
        curr_fasta_out.append(curr_fasta[0]);
        # Si curr_fasta es de largo 2 se guarda directamente
        if len(curr_fasta) <= 2:
            seq_sin_sitio = eliminar_sitio(curr_fasta[1], len_sitio);
        # Si curr_fasta es de largo mayor a 2 junto todo en curr_fasta[1:]
        else:
            print('L_fasta de largo mayor a 2 encontrado.')
            print(curr_fasta)
            seq_juntada = '';
            for j in range(1,len(curr_fasta)):
                seq_juntada = seq_juntada + str(curr_fasta[j]);
            seq_sin_sitio = eliminar_sitio(seq_juntada, len_sitio);
        # Agrego seq_sin_sitio a curr_fasta_out
        curr_fasta_out.append(seq_sin_sitio);
        # Guardo curr_fasta_out en M_out
        M_out.append(curr_fasta_out[:]);
    return M_out


def guardar_M_fasta(M_fasta, nombre_out):
    # Guardo M_fasta en un archivo .fasta

    # Creo el archivo vacio (borra el anterior) y lo abro en modo append
    with open(nombre_out + '.fasta', 'w') as F_out:
        print('Archivo ' + nombre_out + '.fasta creado.')
    with open(nombre_out + '.fasta', 'a') as F_out:
        # Recorro M_fasta
        for i in range(len(M_fasta)):
            curr_fasta = M_fasta[i];
            curr_fasta_str = '';
            # Agrego cada elemento de curr_fasta como una oracion
            for j in curr_fasta:
                curr_fasta_str = curr_fasta_str + j + '\n';
            # Guardo la version str del fasta en F_out
            F_out.write(curr_fasta_str);
    return M_fasta


def leer_fasta(dir_arch, nom_arch):
    # Abre archivo fasta y registra cada entrada empezando por '>'

    # Inicializo la variable que se devuelve
    M_out = [];

    # Abro el archivo
    with open(dir_arch + '.fasta', 'r') as F:
        print('Archivo ' + str(nom_arch) + '.fasta abierto.')
        curr_fasta = [];
        # Recorro las lineas del archivo
        for curr_line in F:
            # Si inicia con '>' se guarda el fasta anterior y se inicia otro
            if curr_line[0] == '>':
                # Solo se registra el fasta anterior si contiene algo
                if len(curr_fasta) > 0:
                    M_out.append(curr_fasta[:]);
                    curr_fasta = [];
                curr_fasta.append(str(curr_line).rstrip());
            # Si no inicia con '>' se agrega curr_line al fasta
            elif len(curr_line) > 0:
                # Solo se agrega curr_line si tiene largo mayor a 0
                curr_fasta.append(str(curr_line).rstrip());
        # Cierro el fasta guardando curr_fasta
        if len(curr_fasta) > 0:
            M_out.append(curr_fasta[:]);
    return M_out


def _main():
    # Funcion para probar funciones en ejecucion del archivo

    # Inicializo la variable que se devuelve
    M_out = [];

    # Extraigo el fasta
    M_fasta = leer_fasta(dir_main, nombre_main);
    ### Display
    #shuffle(M_fasta);
    #for i in M_fasta[:10]:
    #    print(i)

    # Elimino el sitio de largo largo_sitio_main del medio de todos los fasta
    M_fasta = eliminar_sitios_M(M_fasta, largo_sitio_main);
    ### Display
    #shuffle(M_fasta);
    #for i in M_fasta[:10]:
    #    print(i)
    
    # Guardo M_fasta en un fasta nuevo
    # Devuelvo M_fasta
    M_out = guardar_M_fasta(M_fasta, nombre_out_main);
    return M_out


#################################### RESULTADOS ###################################

output_dump = [];

if __name__=='__main__':
    output_dump.append(_main());

