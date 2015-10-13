
#include "IGRF.h"

using namespace qb50;


int IGRF::nLat( int lat )
{
	lat = ( lat + 89 ) % 179;
	return lat < 0 ? 179 + lat : lat;
}


int IGRF::nLon( int lon )
{
	lon %= 360;
	return lon < 0 ? 360 + lon : lon;
}


Eigen::Vector3d IGRF::mag( int lat, int lon )
{
	int off = 179 * nLon( lon ) + nLat( lat );

	return Eigen::Vector3d(
		magData[ off ][ 0 ],
		magData[ off ][ 1 ],
		magData[ off ][ 2 ]
	);
}


/*EoF*/
