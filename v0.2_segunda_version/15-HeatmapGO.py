
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
import goatools
from goatools.obo_parser import GODag
from goatools.base import get_godag
from goatools.gosubdag.gosubdag import GoSubDag

# Generacion de graficos
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np

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
Scripts para generar heatmaps en base a listas de GO terms de ShinyGO

Incluye scripts para modificar la matriz que se da al heatmap

Incluye funciones para seleccionar y filtrar GO terms
'''

#################################### VARIABLES ####################################

# Direcciones de archivos y outputs para ib3 y casa
path_dropbox_casa = 'D:\\Users\\Admin\\Dropbox\\'; 
path_dropbox_ib3 = 'C:\\Users\\emili\\Dropbox\\'; 
path_fasta_casa = 'D:\\ArchivosDoctorado\\Genomas\\'; 
path_fasta_ib3 = 'X:\\Genomas\\'; 
path_output_dump_casa = 'D:\\ArchivosDoctorado\\Output_dump\\'; 
path_output_dump_ib3 = 'X:\\Output_dump\\'; 



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
path_in_main = path_rnaseq_main + 'gene_ontology\\shiny_go\\'; 
path_out_main = path_output_dump_main + 'GOterms\\Graficos\\Heatmaps\\'; 
path_go_obo = path_output_dump_main + 'GeneOntology\\'; 
# Path que dependen de la especie
path_pwm_human = path_dropbox_main + 'Doctorado\\3-Genes transactivados rio abajo\\1-Interpretacion de datos ChIP-seq\\PWM_human\\'; 
path_pwm_mouse = path_dropbox_main + 'Doctorado\\3-Genes transactivados rio abajo\\1-Interpretacion de datos ChIP-seq\\PWM_mouse\\'; 



#################################### FUNCIONES ####################################


def _main():
    # Funcion para probar funciones en ejecucion del archivo

    # Pruebas para heatmaps
    #M_test = np.random.random((18, 8)); 
    #crear_heatmap(M_test, nom_out='', path_out=path_out_main); 

    # Pruebas para extraer datos de ShinyGO
    #nom_shinygo_test = 'enrichment_anderson_100kpb_down_genes_seleccionados'; 
    #M_shinygo = abrir_shinygo(nom_shinygo_test, path_shiny=path_in_main); 
    #for i in M_shinygo:
    #    print(i)

    # Pruebas de pipeline_heatmap_shinygo()
    nom_100k = 'prueba_100k'; 
    nom_10k = 'prueba_10k'; 
    nom_1m = 'prueba_1m'; 
    L_nombres_shinygo_100k = ['enrichment_anderson_100kpb_updown_genes_seleccionados', 'enrichment_anderson_100kpb_up_genes_seleccionados', 
                              'enrichment_anderson_100kpb_down_genes_seleccionados', 'enrichment_dupays_100kpb_updown_genes_seleccionados', 
                              'enrichment_dupays_100kpb_up_genes_seleccionados', 'enrichment_dupays_100kpb_down_genes_seleccionados']; 
    L_nombres_shinygo_10k = ['enrichment_anderson_10kpb_updown_genes_seleccionados', 'enrichment_anderson_10kpb_up_genes_seleccionados', 
                             'enrichment_anderson_10kpb_down_genes_seleccionados', 'enrichment_dupays_10kpb_updown_genes_seleccionados', 
                             'enrichment_dupays_10kpb_up_genes_seleccionados', 'enrichment_dupays_10kpb_down_genes_seleccionados']; 
    L_nombres_shinygo_1m = ['enrichment_anderson_1mpb_updown_genes_seleccionados', 'enrichment_anderson_1mpb_up_genes_seleccionados', 
                            'enrichment_anderson_1mpb_down_genes_seleccionados', 'enrichment_dupays_1mpb_updown_genes_seleccionados', 
                            'enrichment_dupays_1mpb_up_genes_seleccionados', 'enrichment_dupays_1mpb_down_genes_seleccionados']; 
    L_eje_y_shinygo = ['anderson_updown', 'anderson_up', 'anderson_down', 'dupays_updown', 'dupays_up', 'dupays_down']; 
    L_orden_x_100k = [
# Cardiac/Heart
'GO:0035051|Cardiocyte differentiation', 
'GO:0003013|Circulatory system proc.', 
'GO:0072359|Circulatory system development', 
'GO:0060372|Reg. of atrial cardiac muscle cell membrane repolarization', 
'GO:0055013|Cardiac muscle cell development', 
'GO:0097746|Blood vessel diameter maintenance', 
'GO:0048738|Cardiac muscle tissue development', 
'GO:0003018|Vascular proc. in circulatory system', 
'GO:0001568|Blood vessel development', 
'GO:0001944|Vasculature development', 
'GO:2000181|Neg. reg. of blood vessel morphogenesis', 
'GO:1901343|Neg. reg. of vasculature development', 
# Muscle
'GO:0003012|Muscle system proc.', 
'GO:0061061|Muscle structure development', 
'GO:0043500|Muscle adaptation', 
'GO:0090257|Reg. of muscle system proc.', 
'GO:0050881|Musculoskeletal movement', 
'GO:0031444|Slow-twitch skeletal muscle fiber contraction', 
'GO:0014888|Striated muscle adaptation', 
'GO:0055001|Muscle cell development', 
'GO:0042692|Muscle cell differentiation', 
# Cell migration
'GO:0030334|Reg. of cell migration', 
'GO:0016477|Cell migration', 
'GO:0048870|Cell motility', 
'GO:0051674|Localization of cell', 
'GO:0040011|Locomotion', 
'GO:0007155|Cell adhesion', 
'GO:0050879|Multicellular organismal movement', 
'GO:0050918|Positive chemotaxis', 
'GO:0042330|Taxis', 
'GO:0006935|Chemotaxis', 
# Nervous system
'GO:0048699|Generation of neurons', 
'GO:0030182|Neuron differentiation', 
'GO:0022008|Neurogenesis', 
'GO:0007416|Synapse assembly', 
'GO:0007268|Chemical synaptic transmission', 
'GO:0099536|Synaptic signaling', 
'GO:0048666|Neuron development', 
'GO:0033564|Anterior/posterior axon guidance', 
'GO:0007411|Axon guidance', 
'GO:0097485|Neuron projection guidance', 
'GO:0061564|Axon development', 
'GO:0007409|Axonogenesis', 
# Embriology
'GO:0009792|Embryo development ending in birth or egg hatching', 
'GO:0043009|Chordate embryonic development', 
'GO:0009790|Embryo development', 
'GO:0048646|Anatomical structure formation involved in morphogenesis', 
'GO:0009887|Animal organ morphogenesis', 
'GO:0035295|Tube development', 
'GO:0035296|Reg. of tube diameter', 
'GO:2000026|Reg. of multicellular organismal development', 
'GO:0051241|Neg. reg. of multicellular organismal proc.', 
# Other
'GO:0061033|Secretion by lung epithelial cell involved in lung growth', 
'GO:0001541|Ovarian follicle development', 
'GO:0010838|Pos. reg. of keratinocyte proliferation', 
'GO:0046660|Female sex differentiation', 
'GO:0006955|Immune response', 
'GO:1990384|Hyaloid vascular plexus regression', 
'GO:0046533|Neg. reg. of photoreceptor cell differentiation', 
'GO:0045653|Neg. reg. of megakaryocyte differentiation', 
'GO:0045652|Reg. of megakaryocyte differentiation', 
# General
'GO:0044057|Reg. of system proc.', 
'GO:0007169|Transmembrane receptor protein tyrosine kinase signaling pathway', 
'GO:0007167|Enzyme linked receptor protein signaling pathway', 
'GO:0007267|Cell-cell signaling', 
'GO:0006335|DNA replication-dependent nucleosome assembly', 
'GO:0050731|Pos. reg. of peptidyl-tyrosine phosphorylation', 
'GO:0043062|Extracellular structure organization', 
'GO:0032989|Cellular component morphogenesis', 
'GO:0042127|Reg. of cell population proliferation', 
'GO:0006336|DNA replication-independent nucleosome assembly', 
'GO:0006334|Nucleosome assembly', 
'GO:0031497|Chromatin assembly', 
'GO:0034728|Nucleosome organization', 
'GO:0006333|Chromatin assembly or disassembly'
]; 
    L_orden_x_10k = [
# Cardiac/Heart
'GO:0035051|Cardiocyte differentiation', 
'GO:0060372|Reg. of atrial cardiac muscle cell membrane repolarization', 
'GO:0055013|Cardiac muscle cell development', 
'GO:0048738|Cardiac muscle tissue development', 
'GO:1903522|Reg. of blood circulation', 
'GO:0008015|Blood circulation', 
'GO:0014898|Cardiac muscle hypertrophy in response to stress', 
'GO:0014887|Cardiac muscle adaptation', 
# Muscle
'GO:0043502|Reg. of muscle adaptation', 
'GO:0043500|Muscle adaptation', 
'GO:0006942|Reg. of striated muscle contraction', 
'GO:0090257|Reg. of muscle system proc.', 
'GO:0006937|Reg. of muscle contraction', 
'GO:0003012|Muscle system proc.', 
'GO:0014745|Neg. reg. of muscle adaptation', 
'GO:0014706|Striated muscle tissue development', 
'GO:0060537|Muscle tissue development', 
'GO:0006936|Muscle contraction', 
'GO:0061061|Muscle structure development', 
'GO:0031444|Slow-twitch skeletal muscle fiber contraction', 
'GO:0014721|Twitch skeletal muscle contraction', 
'GO:0003009|Skeletal muscle contraction', 
'GO:0050881|Musculoskeletal movement', 
'GO:0014888|Striated muscle adaptation', 
'GO:0006941|Striated muscle contraction', 
'GO:0055001|Muscle cell development', 
'GO:0042692|Muscle cell differentiation', 
# Cell migration
'GO:0030334|Reg. of cell migration', 
'GO:0040012|Reg. of locomotion', 
'GO:0016477|Cell migration', 
'GO:0048870|Cell motility', 
'GO:0051674|Localization of cell', 
'GO:0040011|Locomotion', 
# Nervous system
'GO:0050905|Neuromuscular proc.', 
'GO:0051963|Reg. of synapse assembly', 
'GO:0021702|Cerebellar Purkinje cell differentiation', 
'GO:0021694|Cerebellar Purkinje cell layer formation', 
'GO:0021533|Cell differentiation in hindbrain', 
'GO:0021697|Cerebellar cortex formation', 
# Embriology
'GO:0009887|Animal organ morphogenesis', 
# Other
'GO:0052646|Alditol phosphate metabolic proc.', 
'GO:0007565|Female pregnancy', 
'GO:0050678|Reg. of epithelial cell proliferation', 
'GO:0050673|Epithelial cell proliferation', 
'GO:0007169|Transmembrane receptor protein tyrosine kinase signaling pathway', 
'GO:0050879|Multicellular organismal movement', 
'GO:1902259|Reg. of delayed rectifier potassium channel activity', 
'GO:0045653|Neg. reg. of megakaryocyte differentiation', 
'GO:0045652|Reg. of megakaryocyte differentiation', 
# General
'GO:0030198|Extracellular matrix organization', 
'GO:0043062|Extracellular structure organization', 
'GO:0045229|External encapsulating structure organization', 
'GO:0051241|Neg. reg. of multicellular organismal proc.', 
'GO:0044057|Reg. of system proc.', 
'GO:0007167|Enzyme linked receptor protein signaling pathway', 
'GO:0006335|DNA replication-dependent nucleosome assembly', 
'GO:0006334|Nucleosome assembly', 
'GO:0031497|Chromatin assembly', 
'GO:0034728|Nucleosome organization', 
'GO:0006333|Chromatin assembly or disassembly', 
'GO:0006336|DNA replication-independent nucleosome assembly', 
'GO:0006323|DNA packaging', 
'GO:0065004|Protein-DNA complex assembly', 
'GO:0006338|Chromatin remodeling', 
'GO:0071824|Protein-DNA complex subunit organization', 
'GO:0071103|DNA conformation change'
]; 
    L_orden_x_1m = [
# Cardiac/Heart
'GO:0072359|Circulatory system development', 
'GO:0060372|Reg. of atrial cardiac muscle cell membrane repolarization', 
'GO:0097746|Blood vessel diameter maintenance', 
'GO:0003018|Vascular proc. in circulatory system', 
'GO:0003013|Circulatory system proc.', 
# Muscle
'GO:0090257|Reg. of muscle system proc.', 
'GO:0003012|Muscle system proc.', 
'GO:0042692|Muscle cell differentiation', 
'GO:0061061|Muscle structure development', 
# Cell migration
'GO:0016477|Cell migration', 
'GO:0040011|Locomotion', 
'GO:0048870|Cell motility', 
'GO:0051674|Localization of cell', 
'GO:0030334|Reg. of cell migration', 
'GO:0042330|Taxis', 
'GO:0006935|Chemotaxis', 
# Nervous system
'GO:0007409|Axonogenesis', 
'GO:0061564|Axon development', 
'GO:0048667|Cell morphogenesis involved in neuron differentiation', 
'GO:0031175|Neuron projection development', 
'GO:0048666|Neuron development', 
'GO:0030182|Neuron differentiation', 
'GO:0048699|Generation of neurons', 
'GO:0022008|Neurogenesis', 
'GO:0007416|Synapse assembly', 
'GO:0007411|Axon guidance', 
'GO:0097485|Neuron projection guidance', 
'GO:0048812|Neuron projection morphogenesis', 
# Embriology
'GO:0000904|Cell morphogenesis involved in differentiation', 
'GO:2000026|Reg. of multicellular organismal development', 
'GO:0022603|Reg. of anatomical structure morphogenesis', 
'GO:0009887|Animal organ morphogenesis', 
'GO:0032989|Cellular component morphogenesis', 
'GO:0000902|Cell morphogenesis', 
'GO:0048646|Anatomical structure formation involved in morphogenesis', 
'GO:0035296|Reg. of tube diameter', 
'GO:0120039|Plasma membrane bounded cell projection morphogenesis', 
'GO:0035295|Tube development', 
'GO:0048589|Developmental growth', 
'GO:0035239|Tube morphogenesis', 
'GO:0051093|Neg. reg. of developmental proc.', 
'GO:0040007|Growth', 
'GO:0051094|Pos. reg. of developmental proc.', 
'GO:0042127|Reg. of cell population proliferation', 
'GO:0009888|Tissue development', 
'GO:0045596|Neg. reg. of cell differentiation', 
# Other
'GO:0045653|Neg. reg. of megakaryocyte differentiation', 
'GO:0060441|Epithelial tube branching involved in lung morphogenesis', 
'GO:0048880|Sensory system development', 
'GO:0043010|Camera-type eye development', 
'GO:0001654|Eye development', 
'GO:0150063|Visual system development', 
'GO:0016055|Wnt signaling pathway', 
'GO:0198738|Cell-cell signaling by wnt', 
'GO:1905114|Cell surface receptor signaling pathway involved in cell-cell signaling', 
'GO:0007423|Sensory organ development', 
# General
'GO:0044057|Reg. of system proc.', 
'GO:0000165|MAPK cascade', 
'GO:0006335|DNA replication-dependent nucleosome assembly', 
'GO:0007267|Cell-cell signaling'
]; 
    # Defino agrupaciones de GO_terms
    L_grupos_main = []; 
    # Defino distancia usada
    L_dist_usadas = ['1m', '100k', '10k']; 
    dist_usada = L_dist_usadas[0]; 
    correr = False; 
    # Reviso la distancia usada
    if dist_usada=='1m':
        # 1mpb
        L_nom_shinygo_usado = L_nombres_shinygo_1m; 
        orden_x_usado = L_orden_x_1m; 
        nom_usado = nom_1m; 
    elif dist_usada=='100k':
        # 100kpb
        L_nom_shinygo_usado = L_nombres_shinygo_100k; 
        orden_x_usado = L_orden_x_100k; 
        nom_usado = nom_100k; 
    elif dist_usada=='10k':
        # 10kpb
        L_nom_shinygo_usado = L_nombres_shinygo_10k; 
        orden_x_usado = L_orden_x_10k; 
        nom_usado = nom_10k; 
    else:
        print(dist_usada + ' no es valido para este script.')
        correr = False; 
    # Corro si la distancia usada es valida
    if correr:
        ### Corro esto si quiero generar la lista de GO terms
        orden_x_og = copy.deepcopy(orden_x_usado); 
        #orden_x_usado=[]; 
        ###
        M_heatmap, eje_x = pipeline_heatmap_shinygo(L_nom_shinygo_usado, L_eje_y_shinygo, orden_x=orden_x_usado, nom_heatmap=nom_usado, path_shinygo=path_in_main, path_heatmap=path_out_main); 
        ### Solo corre para generar la lista de GO terms
        if False and orden_x_usado==[]:
            # Creo L_eje_x_final con el mismo largo que la lista anterior
            L_eje_x_final = ['' for n in orden_x_og]; 
            # Recorro la lista con GO ID
            for i in range(len(eje_x)):
                curr_goterm = eje_x[i].split('|'); 
                # Defino go_id y go_name del curr_goterm
                curr_go_id = curr_goterm[0]; 
                curr_go_name = curr_goterm[1]; 
                # Veo que curr_go_name este en orden_x_og
                if curr_go_name in orden_x_og:
                    # Uso index() para determinar la posicion en la que lo registro en L_eje_x_final
                    # Veo que el go_name no sobreescriba ningun dato en L_eje_x_final
                    if L_eje_x_final[orden_x_og.index(curr_go_name)] == '':
                        # Si no hay conflicto, redefino
                        L_eje_x_final[orden_x_og.index(curr_go_name)] = eje_x[i]; 
                    else:
                        # Si hay conflicto, tiro warning y lo agrego uno al lado del otro
                        print('WARNING: GO term repetido')
                        L_eje_x_final[orden_x_og.index(curr_go_name)] += '|'+eje_x[i]; 
            print(L_eje_x_final)

    # Pruebo clamp_M() (YA ESTA)
    if False:
        M_test = [[1,2.0,5,10], [], [-9,0], [50, 18.5]]; 
        print(M_test)
        print(clamp_M(M_test, 5))
    # Primer prueba de GOaTools (YA ESTA)
    if False:
        c = 0; 
        godag = GODag(os.path.join(path_go_obo, 'go-basic.obo')); 
        for i in range(len(L_orden_x_1m)):
            curr_goterm = L_orden_x_1m[i].split('|'); 
            curr_go_id = curr_goterm[0]; 
            curr_goname = curr_goterm[1]; 
            if curr_go_id in godag.keys():
                print(godag[curr_go_id])
            else:
                print('GO ID y term ' + curr_go_id + ' || ' + curr_goname + ' no encontrados.')
                c += 1; 
        print()
        print('Cantidad de GO terms no encontrados: ' + str(c))
    # Segunda prueba de GOaTools
    if True:
        godag = GODag(os.path.join(path_go_obo, 'go-basic.obo'), optional_attrs={'relationship'}); 
        pipeline_goatools(L_orden_x_1m, godag, L_grupos_main); 
    # Prueba de GOaTools para buscar GO terms no encontrados
    # (OBSOLETO) GO:0051270 || Reg. of cellular component movement 
    # (OBSOLETO) GO:0022610 || Biological adhesion
    if False:
        godag = GODag(os.path.join(path_go_obo, 'go-basic.obo'), optional_attrs={'relationship'}); 
        for g in godag.values():
            #if "adhesion" in g.name.lower():
            #    print(g)
            if "cellular component" in g.name.lower() and "regulation" in g.name.lower():
                print(g)
    # godag.keys()
    # godag.values()
    # i.name
    # i.item_id
    # i.level
    # i.depth
    return ''


def abrir_csv(nom_arch, path_arch='', ext='.csv', sep=';', ignore_headers=True):
    # Abre un archivo en formato de matriz con filas por lineas y columnas separadas por sep
    # Devuelve una matriz con una lista de filas, cada una con una lista de columnas
    # ignore_headers ignora la primer fila

    # Inicializo la matriz que se devuelve
    M_out = []; 
    # Defino la direccion del archivo con nom_arch y path_arch
    if path_arch=='':
        filepath = nom_arch + ext; 
    else:
        filepath = os.path.join(path_arch, nom_arch + ext); 
    # Booleano para revisar que se ignore la primera fila
    header_ignorado = not ignore_headers; 
    # Abro el archivo en modo reading
    with open(filepath, 'r') as F_open:
        print('Archivo ' + nom_arch + ' abierto.')
        # Recorro cada fila
        for curr_line in F_open.readlines():
            if header_ignorado:
                # Transformo la linea en lista
                L_line = curr_line.rstrip().split(sep=sep); 
                # Cargo L_line en M_out
                M_out.append(L_line[:]); 
            else:
                # Si header_ignorado empieza sin ser True, se ignora la primera linea
                header_ignorado = True; 
    return M_out


def abrir_shinygo(nom_shiny, L_col=[(3,'fold_enrich'),(4,'pathway'),(5,'url'),(6,'genes')], path_shiny='', ext_shiny='.csv', sep_shiny=';'):
    # Funcion que usa abrir_csv() para extraer la informacion considerada relevante de outputs de ShinyGO

    # Uso abrir_csv() para abrir la matriz de shinygo
    M_shinygo = abrir_csv(nom_shiny, path_arch=path_shiny, ext=ext_shiny, sep=sep_shiny, ignore_headers=True); 
    # Inicializo la matriz que se devuelve
    M_out = []; 
    # Recorro M_shinygo
    for i in range(len(M_shinygo)):
        curr_row = M_shinygo[i]; 
        # Inicializo la lista que se agrega a M_out
        curr_L = []; 
        # Recorro L_col para ver las columnas que me quedo
        for k in range(len(L_col)):
            curr_col = L_col[k]; 
            # Agrego el valor de la columna a curr_L
            curr_L.append(curr_row[curr_col[0]]); 
        # Agrego curr_L a M_out
        M_out.append(curr_L[:]); 
    return procesar_M_shinygo(M_out)


def clamp_M(M_in, val_max, val_min=0):
    # Funcion para limitar los valores de una matriz entre val_max y val_min

    # Inicializo la matriz que se devuelve
    M_out = []; 
    ### Valor para testing
    #val_min_not_0 = 0; 
    ###
    # Recorro M_in
    for i in range(len(M_in)): 
        curr_L = M_in[i]; 
        # Inicializo la lista que cargo a M_out
        L_cargada = []; 
        # Recorro cada valor de curr_L
        for j in range(len(curr_L)): 
            # Limito el valor con val_max y val_min
            val_cargado = max(val_min, min(val_max, float(curr_L[j]))); 
            L_cargada.append(float(val_cargado)); 
            ### Para testing
            #if val_cargado>0 and (val_min_not_0==0 or val_cargado<val_min_not_0):
                #val_min_not_0 = float(val_cargado);
            ### 
        # Cargo L_cargada a M_out
        M_out.append(L_cargada[:]); 
    #print('Valor minimo que no sea 0: ' + str(val_min_not_0))
    return M_out


def crear_heatmap(M_datos, L_eje_x, L_eje_y, nom_out='', path_out='.\\', show_anyways=True):
    # Funcion para generar un heatmap usando una matriz de input con numeros

    ### Configuraciones de pyplot
    # Defino dpi
    dpi_usado = 400; 
    # Defino tamaño de la figura
    plt.figure(figsize=(12, 12), dpi=dpi_usado)
    # Defino colores personalizados
    colors_list = ['#ffffff', '#ffffff', '#ffffff', '#ffcc99', '#ff9900', '#ff9900', '#cc3300', '#cc3300', '#990000', '#550000']; 
    cmap = colors.ListedColormap(colors_list)
    # Funcion de pyplot que genera la matriz
    plt.imshow(M_datos,  cmap=cmap, interpolation='nearest'); 
    # Agrego color bar con formato
    cbar = plt.colorbar(ticks=[0, 1, 2, 3, 4, 5])
    cbar.ax.set_yticklabels(['', 'NS', '2', '3', '4', '5'])
    # Titulo del grafico
    plt.title("Heatmap"); 
    # Titulos de eje x/y
    plt.xlabel("GO terms"); 
    plt.ylabel("Experimentos"); 
    # Ticks de eje x/y
    plt.xticks(range(len(L_eje_x)), L_eje_x, rotation=90); 
    plt.yticks(range(len(L_eje_y)), L_eje_y); 
    # Tight layout
    plt.tight_layout(); 

    ### Defino si muestro o guardo el plot
    # Si nom_out es un string vacio, muestro la matriz
    if nom_out=='':
        plt.show(); 
    # Si nom_out no es un string vacio, lo uso de nombre para generar el heatmap
    else:
        plt.savefig(os.path.join(path_out, str(nom_out)+'.png'), dpi=dpi_usado); 
        if show_anyways:
            # Si show_anyways es True, lo muestro aunque se guarde
            plt.show(); 
    return M_datos


def dividir_goterms(L_goterms, sep='|'):
    # Funcion que recibe una lista de GO terms (go_id | go_name) y la devuelve como matriz

    # Inicializo la matriz que se devuelve
    M_out = []; 
    # Recorro L_goterms
    for i in range(len(L_goterms)):
        curr_goterm = L_goterms[i].split(sep); 
        # Agrego [curr_goterm[0], curr_goterm[1]] a M_out
        M_out.append([str(curr_goterm[0]), str(curr_goterm[1])]); 
    return M_out


def pipeline_goatools(L_goterms, godag, L_groups, level_cutoff=0, depth_cutoff=0):
    # Funcion que recibe una lista de GO terms (go_id | go_name)
    # Los trata de agrupar en GO ids dentro de L_groups y devuelve una matriz con la agrupacion (VER FORMATO)
    # Para la busqueda de ancestros, usa level_cutoff y depth_cutoff (si alguno es 0, no corta con ese criterio)

    # Paso L_goterms a L_goterms_usado con dividir_goterms()
    L_goterms_usado = dividir_goterms(L_goterms); 
    # Esto se usa para GoSubDag, que se usa para ancestors
    optional_relationships = set(); 
    ### PRUEBA
    nom_test = 'test_ancestors1m.csv'; 
    with open(path_out_main + nom_test, 'w') as F_out:
        print('Archivo ' + nom_test + ' creado.')
    ###
    # Defino largo de L_goterms_usado en variable
    largo_goterms = len(L_goterms_usado); 
    ### Display
    print('Largo L_goterms_usado: ' + str(largo_goterms))
    with open(path_out_main + nom_test, 'a+') as F_out:
    # Recorro L_goterms_usado
        for i in range(largo_goterms):
            curr_goterm = L_goterms_usado[i]; 
            curr_goid = curr_goterm[0]; 
            curr_goname = curr_goterm[1]; 
            # Uso GoSubDag para ancestors con regulates (copio de GitHub)
            gosubdag_r0 = GoSubDag([curr_goid], godag, relationships=optional_relationships, prt=None); 
            try:
                curr_ancestors = gosubdag_r0.rcntobj.go2ancestors[curr_goid]; 
                # Creo la lista de ancestors 
                L_ancestors = recursive_ancestors_goatools(curr_ancestors, level_cutoff, depth_cutoff, godag, optional_relationships); 
                ### PRUEBA
                F_out.write(stringify_go_term(godag[curr_goid]) + '\n'); 
            except:
                print('ERROR 1: ID ' + str(curr_goid) + ' | ' + str(curr_goname) + ' no encontrado en busqueda de ancestors.')
                L_ancestors = []; 
                ### PRUEBA
                F_out.write(str(curr_goid) + ';' + str(curr_goname) + '\n'); 
            F_out.write('#;#;#;#\n'); 
            for j in L_ancestors:
                F_out.write(stringify_go_term(godag[j]) + '\n'); 
            F_out.write('\n'); 
            ###
            ### Display
            if (i+1)%5 == 0:
                print('Progreso: ' + str(i+1) + '/' + str(largo_goterms))
            ###
    ### FALTA
    # Agrupar por heart, neuronal, development, etc.
    # Ver que las agrupaciones cubran la mayoria de los GO terms
    # Definir como agrupar los numeros de heatmap asociados cuando haya mas de un numero por agrupacion (agarrar max?)
    # Devolver la lista agrupada con los valores de heatmap asociados
    ###


def pipeline_heatmap_shinygo(L_nom_shinygo, L_eje_y, orden_x=[], val_max_heatmap=5, nom_heatmap='', path_shinygo='', ext_shinygo='.csv', sep_shinygo=';', path_heatmap='.\\'):
    # Funcion para generar un heatmap en base a una lista de tablas con resultados shinygo

    # Columnas de GO term
    col_go_name = 1; 
    col_go_id = 2; 
    # Columna de fold_enrich
    col_fe = 0; 
    # Inicializo la matriz del heatmap
    M_heatmap = []; 
    # Inicializo la lista de valores del eje x
    L_eje_x = orden_x[:]; 
    # Recorro cada uno de los nombres de archivos shinygo y valores de eje y
    for i in range(len(L_nom_shinygo)):
        curr_nom_shinygo = L_nom_shinygo[i]; 
        # Abro el archivo curr_nom_shinygo
        curr_M_shinygo = abrir_shinygo(curr_nom_shinygo, path_shiny=path_shinygo, ext_shiny=ext_shinygo, sep_shiny=sep_shinygo); 
        # Inicializo la lista correspondiente al histograma sacado de M_shinygo
        L_heatmap = [0 for k in L_eje_x]; 
        # Recorro curr_M_shinygo
        for j in range(len(curr_M_shinygo)):
            curr_row_shinygo = curr_M_shinygo[j]; 
            # Defino el GO term con col_go
            curr_go_term = curr_row_shinygo[col_go_id] + '|' + curr_row_shinygo[col_go_name]; 
            # Veo si el GO term esta en L_eje_x
            if curr_go_term in L_eje_x:
                # Si esta, agrego fold_enrich a L_heatmap en posicion correspondiente
                L_heatmap[L_eje_x.index(curr_go_term)] = curr_row_shinygo[col_fe]; 
            else:
                # Si no esta, agrego curr_go_term a L_eje_x y fold_enrich a L_heatmap
                L_eje_x.append(str(curr_go_term)); 
                L_heatmap.append(curr_row_shinygo[col_fe]); 
        # Agrego L_heatmap a M_heatmap
        M_heatmap.append(L_heatmap[:]); 
    # Recorro M_heatmap para asegurarme que todos los L_heatmap tengan el largo de L_eje_x
    for i in range(len(M_heatmap)):
        while len(M_heatmap[i]) < len(L_eje_x):
            M_heatmap[i].append(0); 
    # Preproceso M_heatmap con prepro_M_heatmap()
    M_usada, L_eje_x_usado = prepro_M_heatmap(M_heatmap, L_eje_x, val_max_heatmap); 
    # Creo el heatmap
    M_usada = crear_heatmap(M_usada, L_eje_x_usado, L_eje_y, nom_out=nom_heatmap, path_out=path_heatmap); 
    return M_heatmap, L_eje_x


def prepro_M_heatmap(M_heatmap, L_eje_x, val_max):
    # Funcion para preprocesar la matriz usada para generar el heatmap

    # Uso clamp_M() para limitar M_heatmap a val_max
    M_usada = clamp_M(M_heatmap, val_max); 
    ### FALTA
    # Procesar L_eje_x para agrupar GO terms
    # Solo devolver go_name (para el heatmap)
    ###
    # Provisional, para pasar solo el nombre del GO term a L_eje_x_usado
    L_eje_x_usado = []; 
    for i in L_eje_x:
        L_eje_x_usado.append(i.split('|')[1]); 
    return M_usada, L_eje_x_usado


def procesar_M_shinygo(M_shiny):
    # Funcion para procesar la matriz sacada del output de ShinyGO
    # Asume [fold_enrich, go_term_name, go_term_url, lista_genes] con fold_enrich como string y lista_genes separado por " " (espacio)

    # Inicializo la matriz que se devuelve
    M_out = []; 
    # Recorro M_shiny
    for i in range(len(M_shiny)):
        curr_row = M_shiny[i]; 
        # Inicializo la lista que se carga a M_out
        L_out = []; 
        # Paso curr_row[0] a float
        L_out.append(float(curr_row[0])); 
        # Paso curr_row[1] directamente como string
        L_out.append(str(curr_row[1]).rstrip(' ')); 
        # Paso solo la parte despues de GO: de curr_row[2]
        L_out.append(curr_row[2][-10:]); 
        # Paso los genes como lista separada por " "
        L_out.append(curr_row[3].split('  ')); 
        # Agrego L_out a M_out
        M_out.append(L_out[:]); 
    return M_out


def recursive_ancestors_goatools(curr_ancestors, level_cutoff, depth_cutoff, godag, optional_relationships):
    # Funcion recursiva para devolver todos los ancestros hasta level_cutoff y depth_cutoff de una lista de elementos de goatools

    # Inicializo la lista que se devuelve
    L_out = []; 
    # Recorro curr_ancestors
    for i in curr_ancestors:
        try:
            curr_goterm = godag[i]; 
            # Cargo godag[i] en L_out aunque sea level o depth menor a cutoff
            L_out.append(godag[i].item_id); 
            # Veo si godag[i] es level y depth mayor a cutoff
            if curr_goterm.level > level_cutoff and curr_goterm.depth > depth_cutoff:
                # Uso GoSubDag para ancestors con regulates (copio de GitHub)
                gosubdag_r0 = GoSubDag([curr_goterm.item_id], godag, relationships=optional_relationships, prt=None); 
                # Inicializo la lista que va a tener los elementos que no esten en L_out
                L_recursive_ancestors_usados = []; 
                try:
                    # Intento correr el script copiado de GitHub para obtener ancestors
                    L_recursive_ancestors = gosubdag_r0.rcntobj.go2ancestors[curr_goterm.item_id]; 
                    # Veo cada uno de los ancestors que encuentra
                    for j in L_recursive_ancestors:
                        try:
                            if not (godag[j].item_id in L_out):
                                L_recursive_ancestors_usados.append(godag[j].item_id); 
                        except:
                            print('ERROR 3: ID ' + str(j) + ' no encontrado en godag.')
                except:
                    print('WARNING: ID ' + str(curr_goterm.item_id) + ' no encontrado en busqueda de ancestors. Ya agregado a L_out de todas formas.')
                # Hago recursion con la nueva lista de ancestros
                L_out = L_out + recursive_ancestors_goatools(L_recursive_ancestors_usados, level_cutoff, depth_cutoff, godag, optional_relationships); 
        except:
            print('ERROR 2: ID ' + str(i) + ' no encontrado en godag.')
    return L_out


def stringify_go_term(godag_goterm, sep=';'):
    # Funcion que devuelve un string en formato lo mas corto posible para un dado GO term agarrado con goatools
    str_goterm = str(godag_goterm.level) + sep + str(godag_goterm.depth) + sep + str(godag_goterm.item_id) + sep + str(godag_goterm.name); 
    return str_goterm


def unique_go_terms(L_goterms):
    # Funcion para devolver los elementos unicos de L_goterms
    # Uso i.item_id para identificar elementos unicos

    # Inicializo la lista que se devuelve
    L_out = []; 
    # Inicializo una lista de item_id unicos
    L_unique_item_id = []; 
    # Recorro L_goterms
    for i in range(len(L_goterms)):
        curr_goterm = L_goterms[i]; 
        # Veo si curr_goterm.item_id esta en L_out
        if not (str(curr_goterm.item_id) in L_unique_item_id):
            # Si no esta, agrego item_id a L_unique_item_id
            L_unique_item_id.append(str(curr_goterm.item_id)); 
            # Agrego el elemento completo a L_out
            L_out.append(L_goterms[i]); 
    return L_out


#################################### RESULTADOS ###################################

output_dump = []; 

if __name__=='__main__':
    output_dump.append(_main()); 

