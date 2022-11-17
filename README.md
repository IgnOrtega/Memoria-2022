# Memoria-2022
Memoria para optar al título de ingeniería civil matemática \\

Este repositorio consiste en un problema de grafos, el cuál consiste en estimar un subconjunto de nodos a través de un muestreo y un estimador. Los grafos pueden ser de dos tipos, aleatorio y escala libre. Además la formación del conjunto $H$, es decir el conjunto a determinar su tamaño, puede ser por diversos algoritmos, llamados algoritmos de infección. Por otro lado, existen varios tipos de muestreos programados, en total son 10. La lista de estimadores principales programados son 6. \\

Explicación del código:\\
El archivo "Creacion_grafos" esta configurado para crear un grafo de cada tipo, pero en los experimentos se usaron 100.\\

En la carpeta "Experimentacion" se tienen 3 archivos principales,
-El archivo "generar_muestras" sólo genera muestras, para que estas sean almacenadas en una ubicación.
-El archivo "Experimentacion_con_cargar_muestras_zip" realiza las experimentación cargando las muestras de una ubicación y los resultados de la experimentación de cada grafo la guarda en un dataframe, es decir crea un dataframe para cada grafo con que se realiza la experimentación.
-El archivo "experimentacion_sin_cargar_muestras" realiza las experimentación sin necesidad de cargar las muestras, pues este, también las realiza. Además, los resultados de la experimentación de cada grafo la guarda en un dataframe, es decir crea un dataframe para cada grafo con que se realiza la experimentación.

En la carpeta "graficos" se tienen programas que une los dataframe creados ("Descompresion_dataframe") y los demás muestra de forma gráfica la información obtenida de la fase de experimentación.

Herramientas:\\
La carpeta "librerias" contiene todas las librerias programadas por mi, a excepción de las funciónes para crear grafos, fueron obtenidas de la libreria "networkx" de python. Los archivos están ordenados en función de la funcionalidad, es decir, una librería para crear grafos, para crear muestras, estimar una propiedad del muestreo, conformación del conjunto $H$ y por último, los algoritmos estimadores.\\

