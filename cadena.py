import urllib.request
from tqdm import tqdm
import numpy as np
import requests
from bs4 import BeautifulSoup
import gzip
import shutil
import os
import csv
import matplotlib.pyplot as plt

def download_db(url_, file_):
    try:
        with tqdm(unit="B", unit_scale=True, desc="Descargando", leave=True) as t:
            urllib.request.urlretrieve(url_, file_, reporthook=lambda blocknum, blocksize, totalsize: t.update(blocksize))
        print(f"\nDescarga completa. Archivo guardado como '{file_}'.")
    except Exception as e:
        print(f"Error al descargar el archivo: {e}")
    
def get_file_size(file_path):
    try:
        B_size = os.path.getsize(file_path)
        KB_size= B_size / 1024
        return KB_size
    except FileNotFoundError:
        print(f"El archivo {file_path} no fue encontrado.")
        return None

def cal_mean_sizes(files):
    all_sizes = []

    for file in files: 
        KB_size = get_file_size(file)
        if KB_size is not None:
            all_sizes.append(KB_size)
        
    if len(all_sizes) > 0:
        mean = sum(all_sizes) / len(all_sizes)
        return mean
    else:
        return None

# configuración de URL y archivo de salida
url = "https://ftp.ensemblgenomes.ebi.ac.uk/pub/bacteria/release-57/updated_assemblies.txt"
file_name = "ensembl_db.txt"

# verificación si existe o no la base de datos descargada
if os.path.exists(file_name):
    print("\nEl FTP ya ha sido descargado y almacenado en su dispositivo.")
else:
    print("El archivo FTP no está descargado. Comenzando descarga...")
    download_db(url, file_name)

with open(file_name, 'r') as file:
    lines = file.readlines()

# arreglos vacíos para posterior separación por columnas
has_variation = []
assembly = []
name = []

print("\nSeparando archivos...")

# separación de columnas del archivo
for line in tqdm(lines, desc="Procesando", unit="linea"):
    
    has_variation_ = line.split('\t')[3]
    assembly_ = line.split('\t')[1]
    name_ = line.split('\t')[0]
    has_variation.append(has_variation_)
    assembly.append(assembly_)
    name.append(name_)

output_file_urls_path = "download_files_urls.txt"
output_site_urls_path = "download_sites_urls.txt"

# creación de archivos, un archivo con todas las URL de descarga directa y otro con los sitios de dicha descarga
with open(output_file_urls_path, "w") as download_file_urls:
    print("\nProcesando enlaces de descarga de archivos.")
    for h, n, a in zip(has_variation, name, assembly):
        # corte de "bacteria_X_collection_core...", se elimina todo lo que venga luego de "_core..."
        h_before_core = h.split('core')[0].strip()
        h_before_core = h_before_core.rstrip('_')

        download_file = f"https://ftp.ensemblgenomes.ebi.ac.uk/pub/bacteria/release-57/fasta/{h_before_core}/{n}/ncrna/{n.capitalize()}.{a}.ncrna.fa.gz"
        download_file_urls.write(download_file + "\n")
    print("Fin, enlaces de descarga de archivos.")

with open(output_site_urls_path, "w") as download_site_urls:
    print("\nProcesando enlaces de los sitios de archivos.")
    for h, n, a in zip(has_variation, name, assembly):
        # corte de "bacteria_X_collection_core...", se elimina todo lo que venga luego de "_core..."
        h_before_core = h.split('core')[0].strip()
        h_before_core = h_before_core.rstrip('_')

        download_site = f"https://ftp.ensemblgenomes.ebi.ac.uk/pub/bacteria/release-57/fasta/{h_before_core}/{n}/ncrna/"
        download_site_urls.write(download_site + "\n")
    print("Fin, enlaces de los sitios de archivos.")

print(f"\nLas URLs de descarga han sido guardadas en {output_file_urls_path}")
print(f"Las URLs del sitio han sido guardadas en {output_site_urls_path}")

# conexión con el FTP
response = requests.get(url)

# arreglos necesarios para cálculos y archivos "momentaneos"
file_sizes = []
files_cal_means_up100 = []
files_cal_means_lw100 = []
files_cal_means_lw100_zip = []
file_size_upper100 = []
file_size_lower100 = []
sizes_before_remove = []

with open('download_sites_urls.txt', 'r') as archivo:
    # Utiliza next para omitir la primera línea
    next(archivo) # salta la prímera linea del archivo

    # contadores necesarios
    count = 1
    count_upper100 = 0
    count_lower100 = 0
    print("\nExtrayendo peso del archivo desde FTP.")

    # bucle de los 20 primeros archivos/líneas (según lo acordado previamente)
    for url in tqdm(archivo, desc="Procesando URLs", unit="URL"):
        if count > 20:
            break  # Sale del bucle después de las primeras 20 líneas

        url = url.strip()  # Elimina espacios en blanco y caracteres de nueva línea
        response = requests.get(url)

        if response.status_code == 200:
            soup = BeautifulSoup(response.text, 'html.parser')
            table_html = soup.find('table', {'id': ''})
            
            # verifica si existe una tabla HTML dado que el código trabaja con la conexión directa del FTP y verifica el contenido real de la tabla HTML
            if table_html:
                # separación de filas y columnas 
                for fila in table_html.find_all('tr'):
                    columnas = fila.find_all(['th', 'td'])
                    if len(columnas) >= 2:
                        segunda_celda = columnas[1]
                        if segunda_celda.text.strip().endswith('.gz'):
                            datos_fila = [columna.text.strip() for columna in columnas]
                            ultimo_valor = datos_fila[3]
                            nombre = datos_fila[1]
                            print(nombre)

                            try:
                                size_k = ''.join(char if char.isdigit() or char == '.' else ' ' for char in ultimo_valor)
                                size = float(size_k.split()[0]) if size_k else 0  # convierte a número el tamaño "char"
                                file_sizes.append(size)
                            except ValueError:
                                print("Error.")

                            final_url_download = f"{url}{nombre}"

                            print(f'\n{count}. URL: {final_url_download} | Peso del archivo ".gz": {ultimo_valor}')

                            download_file_name = f"{count}_{nombre}"
                            download_db(final_url_download, download_file_name)

                            file_gz = f"{download_file_name}"
                            unzip_file = f"{download_file_name}_unzip"

                            current_path = os.getcwd()
                         
                            if not os.path.exists(current_path):
                                os.makedirs(current_path)

                            with gzip.open(file_gz, 'rb') as f_comprimido:
                                with open(unzip_file, 'wb') as f_descomprimido:
                                    shutil.copyfileobj(f_comprimido, f_descomprimido)

                            with open(unzip_file, 'r', encoding='utf-8') as f_descomprimido:
                                file_content = f_descomprimido.read()

                            # contar la cantidad de '>'
                            flag = file_content.count('>')
                            if(flag >= 100):
                                print(f"\nEl archivo ''{unzip_file}'' contiene una cantidad de: {flag} y es útil.")
                                file_size_upper100.append(size)
            
                                if os.path.exists(file_gz):
                                    os.remove(file_gz)   

                                sizefile_upper100 = get_file_size(unzip_file)
                                print(f"Peso del archivo: {sizefile_upper100} KB.") 
                                count_upper100 += 1

                                files_cal_means_up100.append(unzip_file)

                            elif(flag < 100):
                                print(f"\nEl archivo ''{unzip_file}'' contiene una cantidad de: {flag}. Es insuficiente y será eliminado.")
                                file_size_lower100.append(size)
                                sizefile_lower100 = get_file_size(unzip_file)
                                sizes_before_remove.append(sizefile_lower100)
                                print(f"Peso del archivo: {sizefile_lower100} KB.")
                                files_cal_means_lw100.append(unzip_file)
                                files_cal_means_lw100_zip.append(file_gz)
                                count_lower100 += 1
                            else:
                                print("Error.")     

                            count += 1
                            flag2 = cal_mean_sizes(files_cal_means_up100)
                            flag3 = cal_mean_sizes(files_cal_means_lw100)
            else:
                print(f"{count}. URL: {url} | No se encontró la tabla.")
                count += 1
        else:
            print(f"{count}. URL: {url} | Error: {response.status_code}")
            count += 1

print("\n== ESTADÍSTICAS ==")

if flag2 is not None:
    percentage = (count_upper100 / 20) * 100
    print(f"\nEstadísticas para archivos >= 100 '>' contados: \nCantidad de archivos: {count_upper100} | Porcentaje: {percentage}%")
    print(f"Peso promedio de los archivos: {flag2:.2f} KB.")
else:
    print("Error.")

if flag3 is not None:
    percentage = (count_lower100 / 20) * 100
    print(f"\nEstadísticas para archivos < 100 '>' contados: \nCantidad de archivos: {count_lower100} | Porcentaje: {percentage}%")
    print(f"Peso promedio de los archivos: {flag3:.2f} KB.")
else:
    print("Error.")

for unzip_file, zip_file in zip(files_cal_means_lw100, files_cal_means_lw100_zip):
    if os.path.exists(unzip_file):
        os.remove(unzip_file)
    
    if os.path.exists(zip_file):
        os.remove(zip_file)

print("\nFin. Han sido almacenados sólo los archivos que contienen en total '>' 100 o más.")

archivos_y_pesos_upper100 = zip(files_cal_means_up100, file_size_upper100)
archivos_y_pesos_lower100 = zip(files_cal_means_lw100, sizes_before_remove)

archivos_upper100, pesos_upper100 = zip(*archivos_y_pesos_upper100)
archivos_lower100, pesos_lower100 = zip(*archivos_y_pesos_lower100)

archivos_upper100_nombres = [os.path.basename(archivo) for archivo in archivos_upper100]
archivos_lower100_nombres = [os.path.basename(archivo) for archivo in archivos_lower100]

plt.figure(figsize=(12, 6))
plt.bar(archivos_upper100_nombres, pesos_upper100, color='green', label='Archivos >= 100')
plt.bar(archivos_lower100_nombres, pesos_lower100, color='blue', label='Archivos < 100')
plt.xlabel('Archivos')
plt.ylabel('Peso (KB)')
plt.title('Pesos de archivos (Fuente Normal)')
plt.xticks(rotation=45, ha='right', fontsize=10)
plt.legend()
plt.tight_layout()
plt.savefig('grafico_pesos_archivos_fuente_normal.png')
plt.show()

plt.figure(figsize=(12, 6))
plt.bar(archivos_upper100_nombres, pesos_upper100, color='green', label='Archivos >= 100')
plt.bar(archivos_lower100_nombres, pesos_lower100, color='blue', label='Archivos < 100')
plt.xlabel('Archivos')
plt.ylabel('Peso (KB)')
plt.title('Pesos de archivos (Fuente Pequeña)')
plt.xticks(rotation=45, ha='right', fontsize=8)
plt.legend()
plt.tight_layout()
plt.savefig('grafico_pesos_archivos_fuente_pequena.png')
plt.show()