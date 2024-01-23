#!/bin/bash
for i in {0..128};do

    nombre_carpeta=bacteria_${i}_collection_NCRNA

    url="https://ftp.ensemblgenomes.ebi.ac.uk/pub/bacteria/release-57/fasta/bacteria_${i}_collection/"
    p=$(curl -s $url | grep -Ev "table|head|body|title|html|index|DOCTYPE|h1|th|PARENTDIR" | awk -F'/' '{print $5}' | sed 's/[^a-zA-Z_0-9]//g')

    echo $"COMIENZO PARA bateria_${i}_collection ."
    for x in $p;do
        url2="https://ftp.ensemblgenomes.ebi.ac.uk/pub/bacteria/release-57/fasta/bacteria_${i}_collection/${x}/ncrna/"
        archivo_gz=$(curl -s "$url2" | grep -Eo 'href="[^"]*\.gz"' | sed 's/href="//; s/"//')
        if [ -d "$nombre_carpeta" ]; then
            echo "La carpeta '$nombre_carpeta' ya existe."
            cd $nombre_carpeta

            if [ -f "$archivo_gz" ]; then
                echo "El archivo ${archivo_gz} ya existe."
                cd ..
            else
                # Descarga
                echo $"Descargando "${archivo_gz}""
                wget "https://ftp.ensemblgenomes.ebi.ac.uk/pub/bacteria/release-57/fasta/bacteria_${i}_collection/${x}/ncrna/$archivo_gz"
                gunzip $archivo_gz
                cd ..
            fi

        else
            # Crear la carpeta si no existe
            mkdir "$nombre_carpeta"
            echo "Se ha creado la carpeta '$nombre_carpeta'."
            cd $nombre_carpeta

            if [ -f "$archivo_gz" ]; then
                echo "El archivo ${archivo_gz} ya existe."
                cd ..
            else
                # Descarga
                echo $"Descargando "${archivo_gz}""
                wget "https://ftp.ensemblgenomes.ebi.ac.uk/pub/bacteria/release-57/fasta/bacteria_${i}_collection/${x}/ncrna/$archivo_gz"
                cd ..
            fi
        fi
    done
    echo $"FIN"
done