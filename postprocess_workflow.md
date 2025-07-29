# Documentación de `postprocess_workflow.R`

Este documento explica el flujo básico de postprocesamiento de las salidas generadas por mHM. El script `postprocess_workflow.R` emplea funciones de [`R/utils.R`](R/utils.R) para generar archivos listos para su análisis.

## Requisitos previos

- Haber completado el preprocesamiento del dominio y contar con las salidas de mHM.
- Disponer de las librerías de R usadas en el script (`terra`, `tidyverse`, `stars`, `sf`, `yaml`).

## Uso básico

1. Ajuste `domain_path` y `config_name` al comienzo del script para apuntar al dominio y al archivo de configuración.
2. Ejecute desde R:

```R
source("R/postprocess_workflow.R")
```

Se creará la carpeta `FIGS/` en caso de no existir y se escribirán allí las figuras generadas.

## Tareas realizadas

- **write_output**: exporta variables mensuales desde `mHM_Fluxes_States.nc` recortando al ROI definido en la configuración.
- **write_clim**: procesa las forzantes meteorológicas a resolución mensual.
- **calculate_annual_mean** y **visualize_annual_mean**: calculan promedios anuales y guardan una imagen resumen.
- **get_qmm_table** y **get_qm3s_table**: extraen series de caudal en milímetros y metros cúbicos por segundo.
- Combina ambas tablas en `streamflow/streamflow_data.csv` dentro de la carpeta `OUT`.
- A partir de esos datos calcula métricas de desempeño (R2 y KGE) y genera mapas
  y gráficos comparativos del caudal simulado versus observado.
- Las funciones **read_and_merge_streamflow**, **setup_output_paths** y los
  distintos gráficos de evaluación de caudal se encuentran ahora en
  [`R/utils.R`](R/utils.R) y pueden reutilizarse en otros scripts.

## Salidas generadas

- Archivos NetCDF mensuales en `OUT`.
- Gráficos en `FIGS/annual_output_default.png`.
- Tabla de caudal observada y simulada en `OUT/streamflow/`.
- Mapas, boxplots y distribuciones acumuladas de las métricas de evaluación en `FIGS/`.

Estas utilidades permiten analizar rápidamente tanto las forzantes como las salidas del modelo.
