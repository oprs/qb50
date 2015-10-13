Table des champs magnétiques du NOAA (IGRF-mag.cpp)

    % make EIGEN=../eigen3
    g++ -Wall -Wextra -std=c++11 -I../eigen3 -o IGRF.o -c IGRF.cpp
    g++ -Wall -Wextra -std=c++11 -I../eigen3 -o IGRF-mag.o -c IGRF-mag.cpp
    g++ -Wall -Wextra -std=c++11 -I../eigen3 -o magtest.o -c magtest.cpp
    g++ -o magtest IGRF.o IGRF-mag.o magtest.o

    % ./magtest 48 2
    lat: 48 (norm: 137)
    lon: 2 (norm: 2)
    18792.4
     -144.4
      37128
