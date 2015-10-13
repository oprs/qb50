
#ifndef _QB50_IGRF_H
#define _QB50_IGRF_H

#include <Eigen/Core>


namespace qb50 {

	class IGRF
	{
		public:

			static int nLat( int lat );
			static int nLon( int lon );

			static Eigen::Vector3d mag( int lat, int lon );

			static const float magData[ 360 * 179 ][ 3 ];
	};

} /* qb50 */


#endif /* _QB50_IGRF_H */

/*EoF*/
