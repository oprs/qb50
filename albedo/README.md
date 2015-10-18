Adrian:

>En prenant l'hypothèse d'un albédo de 25% les estimations me semblent satisfaisantes (les pires erreurs que j'ai trouvé sont de 0.5 rad). Je vous envoie le code C++ si vous voulez essayer, il faut le compiler avec le dossier des headers Eigen, puis l'utiliser ainsi:
>
    ./testsalbedo <xSoleil> <ySoleil> <zSoleil> <xAlbedo> <yAlbedo> <zAlbedo>

>où les paramètres donnent la direction des rayons solaires et réfléchis dans le référentiel du satellite (renormalisés dans le programme).

    % make EIGEN=../eigen3/
    g++ -Wall -Wextra -std=c++11 -I../eigen3/ -o testsalbedo.o -c testsalbedo.cpp
    g++ -o testsalbedo testsalbedo.o

    % ./testsalbedo -1 1 0 1 1 0
    vmes
        0.25
           0
       0.125
       0.125
    0.883883
           0
           0
           0
    0.707107
    vest
    -0.316228
     0.948683
            0
    erreur : 0.463648 rad
