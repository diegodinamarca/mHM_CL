# mHM_CL

Este repositorio contiene utilidades en R y Python para preparar los datos de entrada del modelo hidrológico distribuido mHM.

Para una explicación detallada de cómo ejecutar el script principal de preprocesamiento consulte [run_preprocessing.md](run_preprocessing.md).

Para visualizar promedios anuales de largo plazo de las salidas y forzantes del modelo se incluye el script [`R/visualize_annual_outputs.R`](R/visualize_annual_outputs.R).

La función `visualize_annual_outputs` permite opcionalmente entregar la ruta de un archivo ROI mediante el argumento `roi_file` para recortar y enmascarar los resultados a dicha zona en lugar de utilizar el ROI definido en la configuración.

En [`R/utils.R`](R/utils.R) se encuentran algunas funciones utiles para la ejecucion del pre y post procesamiento.

La función `write_output` permite exportar variables mensuales o anuales a TIFF o NetCDF desde el output del modelo `mHM_Fluxes_States.nc`. De forma opcional se puede especificar la carpeta de salida mediante el argumento `out.opt`, reemplazando la definida en el archivo `preprocess_config.yaml`. El resultado es un único archivo con una banda por cada periodo de tiempo solicitado.

Para visualizar promedios anuales de los archivos mosaiceados se puede utilizar el script [`R/visualize_mosaic_outputs.R`](R/visualize_mosaic_outputs.R). Este script permite seleccionar un archivo de ROI para enmascarar los resultados.

Para visualizar conjuntamente las salidas mosaiceadas y las forzantes se puede utilizar:

```r
domain_folders <- dir(pattern = "domain_zone", full.names = TRUE)
visualize_mosaic_full_outputs(out_dir = "domain_chile/OUT",
                              domains = domain_folders)
```

 La función `process_meteo_variable` permite procesar las forzantes meteorológicas diarias almacenadas en la carpeta `meteo`. El flujo calcula una sola vez los valores mensuales y luego utiliza ese resultado para obtener los datos anuales y el promedio anual de largo plazo. Durante la ejecución se muestran mensajes que indican cada etapa del proceso. Los archivos resultantes se escriben en `out/meteo/<var_name>`.
