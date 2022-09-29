Métodos de Simulación
Taller 1 - Punto 5

Esta carpeta contiene el codigo necesario para obtener los resultados pedidos entre programas de C++ y gnuplot, un Makefile que facilita la compilación y los resultados ya obtenidos en la carpeta data.

IMPORTANTE: Para correr los puntos f en adelante es necesario tener instalada la librerai GNU Scientific Library (https://www.gnu.org/software/gsl/)

INSTRUCCIONES:
- el Makefile compila y corre el codigo, su uso es:
     make FILE_NAME (sin extension)
  y guarda la salida estandar (cout) en la carpeta data

-Para generar las gráficas:
  -Si en el código descomentó las opciones para animacion:
  gnuplot data/FILE_NAME.dat

  -Si la opcion para dibujar trayectoria está descomentada:
  gnuplot FILE_NAME.gp

  -el comando
    make clean
  busca y elimina todos los archivos con extensiones .x .dat y .log

  -Opcionalmente, se puede usar el comando
    make graph gp=FILE_NAME
  que compila FILE_NAME.cpp, corre el resultado y ejecuta
  gnuplot FILE_NAME.gp

Los archivos .gp guardan sus resultados en archivos .pdf en la carpeta data

Para hacer el histograma del punto f) hay un programa específico: punto_5f_histogram
