# mHM_CL

Este repositorio contiene utilidades en R y Python para preparar los datos de entrada del modelo hidrológico distribuido mHM.

Para una explicación detallada de cómo ejecutar el script principal de preprocesamiento consulte [run_preprocessing.md](run_preprocessing.md).

Para visualizar promedios anuales de las salidas y forzantes del modelo se incluye el script [`R/visualize_annual_outputs.R`](R/visualize_annual_outputs.R).  Las funciones auxiliares `annual_mean`, `daily_to_monthly` y `monthly_to_yearly` se encuentran en [`R/utils.R`](R/utils.R). `visualize_annual_outputs` puede opcionalmente enmascarar los resultados usando el `roi_file` definido en el archivo de configuración.
