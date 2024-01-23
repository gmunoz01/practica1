#!/bin/bash

url="https://ftp.ensemblgenomes.ebi.ac.uk/pub/bacteria/release-57/updated_assemblies.txt"

# conexión con curl y filtro con awk de las filas, teniendo el sitio web de donde encontrarlas y el link directo de descarga
download_file=$(curl -s $url | awk 'BEGIN{FS="\t"}($4!="database"){print "https://ftp.ensemblgenomes.ebi.ac.uk/pub/bacteria/release-57/fasta/" $4 "/" $1 "/ncrna/" $1 "." $2 ".ncrna.fa.gz"}' | sed 's/_core[^/]*//' | sed 's/\(.*\/\)\([a-z]\)/\1\U\2/' | sed 's/Ncrna/ncrna/g')
download_site=$(curl -s $url | awk 'BEGIN{FS="\t"}($4!="database"){print $4 "\t" "https://ftp.ensemblgenomes.ebi.ac.uk/pub/bacteria/release-57/fasta/" $4 "/" $1 "/ncrna/"}' | sed 's/_core[^\t]*//' | sed 's/_core[^/]*//')

# guardado de los datos en archivos txt tabulados
printf "%s" "$download_file" > download_files.txt
printf "%s" "$download_site" > download_sites.txt

echo $"FIN"