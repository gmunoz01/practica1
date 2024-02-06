# Practica n°1
* Universidad Tecnológica Metropolitana
* Ingeniería Civil en Bioinformática - Universidad de Talca
* Gerardo Muñoz

El presente Github corresponde a los archivos creados durante la práctica profesional I.

## Descripción

Los archivos fueron creados para trabajar la base de datos de Ensembl, una página FTP en donde se pueden encontrar todas las líneas de más de tres mil bacterias. El link que se utilizó desde un comienzo se puede encontrar [aquí](https://ftp.ensemblgenomes.ebi.ac.uk/pub/bacteria/release-57/updated_assemblies.txt).

A partir del enlace anterior, se realizó tres códigos distintos en donde se generaban distintos URL personalizados buscando un orden entre ellos para poder así obtener una descarga directa en base al enlace generado de forma automática.

### Pre-requisitos para Python

_Las liberías utilizadas fueron para el desarrollo del código en Python:_

* [Tqdm](https://tqdm.github.io/) - Barra de progreso en terminal.
* [Numpy](https://numpy.org/doc/stable/) - Herramienta de cálculos estadísticos.
* [BeatifulSoup](https://pypi.org/project/beautifulsoup4/) - Herramienta para visualización de HTML en terminal, útil para separar las columnas sin descargar la página por completo.
* [gzip](https://docs.python.org/3/library/gzip.html) - Herramienta para descomprimir archivos .zip descargados.
* [Matplotlib](https://matplotlib.org/) - Herramienta para la visualización gráfica de datos.
