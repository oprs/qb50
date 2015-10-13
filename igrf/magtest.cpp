
#include <iostream>
#include "IGRF.h"

using namespace std;
using namespace qb50;


int main( int argc, char *argv[] )
{
	if( argc < 3 ) {
		cerr << "usage: " << argv[0] << " latitude longitude" << endl;
		exit( 1 );
	}

	int lat = atoi( argv[ 1 ]);
	int lon = atoi( argv[ 2 ]);

	cout << "lat: " << lat << " (norm: " << IGRF::nLat( lat ) << ")" << endl;
	cout << "lon: " << lon << " (norm: " << IGRF::nLon( lon ) << ")" << endl;

	cout << IGRF::mag( lat, lon ) << endl;

	return 0;
}


/*EoF*/
