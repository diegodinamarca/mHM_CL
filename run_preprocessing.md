# Documentación de `run_preprocessing.R`

Este script prepara todos los insumos necesarios para ejecutar el modelo mHM. A continuación se describen los pasos principales y la configuración requerida.

## Requisitos previos

- **R** con los paquetes indicados en el propio script (`reticulate`, `yaml`, `terra`, `sf`, `magrittr`, `tidyverse`, `whitebox`, `lubridate`).
- **Python** configurado dentro de un entorno virtual que incluya la instalación de mHM y los paquetes adicionales (`subprocess`, `os`, `pyyaml`, `re`).
- **mHM** instalado en un entorno **conda**. Se recomienda seguir la [guía oficial](https://mhm-ufz.org/guides/) para completar la instalación.
- Haber clonado el repositorio oficial de mHM para disponer del script `create_latlon.py` ubicado en `pre-proc`. En `preprocess_config.yaml` se debe definir la clave `"latlon_script_folder"` apuntando a esta carpeta.
- Un archivo `preprocess_config.yaml` dentro de la carpeta del dominio donde se guardarán las salidas.
- Para crear un nuevo dominio cree una carpeta (por ejemplo `new_domain`), copie allí `preprocess_config.yaml` con las rutas de las forzantes a procesar y especifique esta ruta en `run_preprocessing.R` mediante la variable `domain_path`.

## Uso básico

1. Defina la ruta del dominio modificando la variable `domain_path` al comienzo del script. En esta carpeta se escribirá todo el preprocesamiento.
2. Ajuste la ruta del intérprete de Python en la llamada `use_python()` para apuntar a su entorno virtual.
3. Ejecute el script desde R:

```R
source("R/run_preprocessing.R")
```

### Ejecución en modo batch

Si se desea lanzar el preprocesamiento directamente desde la línea de comandos,
se puede utilizar el script `run_preprocessing_batch.R` indicando el dominio
como argumento:

```bash
Rscript R/run_preprocessing_batch.R /ruta/al/dominio
```

Esta versión es equivalente al script original pero permite definir la ruta del
dominio sin editar el código.

El script creará automáticamente las carpetas definidas en el archivo de configuración y llamará a cada una de las funciones de preprocesamiento (clima, LAI, DEM, uso de suelo, suelos, geología, caudales, etc.). Además, ejecutará scripts de Python para generar los archivos `latlon` y actualizar los archivos de parámetros de mHM.

## Contenido del archivo de configuración

El archivo `preprocess_config.yaml` define rutas de entrada y de salida, variables climáticas a procesar y parámetros espaciales. Un ejemplo reducido es el siguiente:

```yaml
dem_file: ./DATA/RAST/Topography/DEM.tif
roi_file: ./DATA/SHP/roi.geojson
variables_clim:
  pr:
    input_dir: ./DATA/RAST/Clim/Pr
  tmax:
    input_dir: ./DATA/RAST/Clim/Tmax
out_folder: ./OUT
meteo_folder: ./meteo
lai_folder: ./lai
```

Cada función de preprocesamiento utilizará estas rutas para leer datos de entrada y escribir los resultados correspondientes.

## Salidas generadas

Al finalizar la ejecución se habrán creado, dentro de la carpeta del dominio, las siguientes subcarpetas principales (según lo definido en el YAML de configuración):

- `meteo/` – Archivos NetCDF con las variables meteorológicas recortadas.
- `lai/` – Series temporales de índice de área foliar.
- `morph/` – Productos derivados del DEM y de la red de drenaje.
- `landcover/` – Mapas de cobertura de suelo reproyectados.
- `gauges/` – Información de estaciones de aforo y archivos de parámetros actualizados.

Estas salidas son los insumos necesarios para lanzar corridas de mHM posteriormente.

## Consejos adicionales

- Revise los mensajes de consola para verificar el progreso de cada etapa.
- Si alguna función genera archivos temporales, puede eliminarlos automáticamente estableciendo `remove_temp = TRUE` en las llamadas del script.
- El script asume que las rutas de entrada definidas en `preprocess_config.yaml` existen y contienen los datos requeridos.
