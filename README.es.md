# ChIP-seq-analysis-scripts
Compilación de funciones usadas para analizar sitios de unión de factores de transcripción en base a experimentos ChIP-seq y RNA-seq

[![pick](https://img.shields.io/badge/lang-pick-red.svg)](https://github.com/EmilioKolo/ChIP-seq-analysis-scripts/blob/main/README.md)
[![en](https://img.shields.io/badge/lang-EN-green.svg)](https://github.com/EmilioKolo/ChIP-seq-analysis-scripts/blob/main/README.en.md)
___

Las versiones finales de los programas utilizados se encuentran en el archivo pipeline.py de la carpeta v1.0_scripts_simplificados

El archivo pipeline.py contiene 3 funciones principales. Estas incluyen comentarios en español, aunque sin caracteres que contengan tildes.

* pipeline_generador(): 
Genera archivos con picos de ChIP-seq, sitios de unión y genes cerca de los picos/sitios de unión. 
Toma inputs de archivos .bed para ChIP-seq y de .csv para RNA-seq

* pipeline_meme_chip(): 
Genera archivos con secuencias en formato .fasta para mandar a MEME-ChIP. 
Toma como parte de su input los archivos generados por pipeline_generador()

* pipeline_otros_tf(): 
Genera archivos con sitios de unión de varios tf cerca de sitios de unión de pipeline_generador(). 
Toma como parte de su input los archivos generados por pipeline_generador()
