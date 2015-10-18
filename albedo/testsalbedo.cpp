#include <iostream>
#include <Eigen/Dense>
#include <cmath>

using namespace std;
using namespace Eigen;

double max(double a, double b){ return (a>b)? a : b;}

int main(int argc, char* argv[])
{
	MatrixXf h(3,9);
	h << 1,  1, 1,  1, 0,        0,       0,        0,      -sqrt(2),
		 1, -1, 0,  0, sqrt(2), -sqrt(2), 0,        0,       0,
		 0,  0, 1, -1, 0,        0,       sqrt(2), -sqrt(2), 0;
	h /= sqrt(2);
	
	if (argc < 7)
	{
		cout << "Utilisation : "<<argv[0]<<" <xSoleil> <ySoleil> <zSoleil> <xAlbedo> <yAlbedo> <zAlbedo>"<<endl;
		return 0;
	}
	
	//Partie simulation
	Vector3f vr(atof(argv[1]), atof(argv[2]),atof(argv[3]));  //Vecteur reel (pour simuler)
	vr = vr/vr.norm();
	Vector3f va(atof(argv[4]), atof(argv[5]),atof(argv[6]));  //Vecteur albedo
	va = va/(4.*va.norm());
	VectorXf vmes(9);  //Contient les mesures pour chaque capteur
	int i = 0;
	for (i=0; i<9; ++i)
		vmes(i) = max(0,vr.transpose()*h.col(i)) + max(0,va.transpose()*h.col(i));
	cout << "vmes" << endl << vmes << endl;
	//Fin de la partie simulation
	
	Vector3f vest(0,0,0);
	MatrixXf m(3,3);
	m << 1./3., 0, 0,
		 0, 1./3., 0,
		 0, 0, 1./3.;   //m = (h*h.transpose())^(-1)
	vest = h*vmes;
	vest = m*vest;
	vest = vest/vest.norm();
	cout << "vest" << endl << vest << endl;
	cout << "erreur : " << acos(vest.dot(vr)) << " rad" << endl;
    return 0;
}
