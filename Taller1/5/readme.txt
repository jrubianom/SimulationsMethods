Métodos de Simulación
Taller 1 - Punto 5
INSTRUCCIONES:
a)
- Se recomienda compilar así:
g++ -std=c++17 punto_5a.cpp -o punto_5a.x
./punto_5a.x > data/punto_5a.dat

-Para generar las gráficas:
  -Si en el código descomentó las opciones para animacion:
  gnuplot data/punto_5a.dat

  -Si la opcion para dibujar trayectoria está descomentada:
  gnuplot punto_5a.gp
