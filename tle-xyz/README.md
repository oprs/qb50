Conversion TLE vers XYZ

    % make
    g++ -Wall -Wextra -g -std=c++11 -o FormParser.o -c FormParser.cpp
    g++ -Wall -Wextra -g -std=c++11 -o main.o -c main.cpp
    g++ -Wall -Wextra -g -std=c++11 -o sgp4sdp4.o -c sgp4sdp4.cpp
    g++ -Wall -Wextra -g -std=c++11 -o sgp_math.o -c sgp_math.cpp
    g++ -Wall -Wextra -g -std=c++11 -o sgp_obs.o -c sgp_obs.cpp
    g++ -Wall -Wextra -g -std=c++11 -o sgp_time.o -c sgp_time.cpp
    g++ -o xyz FormParser.o main.o sgp4sdp4.o sgp_math.o sgp_obs.o sgp_time.o

    % ./xyz iss.txt  
                 Satellite number: 25544
                       Epoch year: 15
            Epoch day of the year: 264.555
          1st. derivative of mm/2: 7.247e-05
          2nd. derivative of mm/6: 0
                  BSTAR drag term: 0.00011673
                 Satellite number: 25544
            Inclination (degrees): 51.6466
                R.A.A.N (degrees): 322.924
                     Eccentricity: 0.0005168
        Arg. of perigee (degrees): 302.537
           Mean anomaly (degrees): 206.171
            Mean motion (rev/day): 15.541
        Revolution number @ epoch: 96303

    LAT: -21.0317
    LON: -106.683

    LAT: -21.0234
    LON: -106.678

    LAT: -21.0151
    LON: -106.674
