
#include <iostream>
#include <fstream>

extern "C" {
 #include <unistd.h>
}

#include "FormParser.h"
#include "sgp4sdp4.h"

using namespace qb50;


static void handleT1( sat_t *sat, const uint8_t *x, size_t len )
{
   T1Form t1;

   FormParser::_parseT1Form( &t1, x, len );
   std::cout << "             Satellite number: " << t1.sat << "\r\n";
   std::cout << "                   Epoch year: " << t1.eyr << "\r\n";
   std::cout << "        Epoch day of the year: " << t1.edy << "\r\n";
   std::cout << "      1st. derivative of mm/2: " << t1.d1m << "\r\n";
   std::cout << "      2nd. derivative of mm/6: " << t1.d2m << "\r\n";
   std::cout << "              BSTAR drag term: " << t1.bdt << "\r\n";

   sat->tle.epoch_year = t1.eyr;
   sat->tle.epoch_day  = t1.edy;
   sat->tle.xndt2o     = t1.d1m;
   sat->tle.xndd6o     = t1.d2m;
   sat->tle.bstar      = t1.bdt;
}


static void handleT2( sat_t *sat, const uint8_t *x, size_t len )
{
   T2Form t2;

   FormParser::_parseT2Form( &t2, x, len );
   std::cout << "             Satellite number: " << t2.sat << "\r\n";
   std::cout << "        Inclination (degrees): " << t2.inc << "\r\n";
   std::cout << "            R.A.A.N (degrees): " << t2.ran << "\r\n";
   std::cout << "                 Eccentricity: " << t2.ecc << "\r\n";
   std::cout << "    Arg. of perigee (degrees): " << t2.aop << "\r\n";
   std::cout << "       Mean anomaly (degrees): " << t2.man << "\r\n";
   std::cout << "        Mean motion (rev/day): " << t2.mmo << "\r\n";
   std::cout << "    Revolution number @ epoch: " << t2.rev << "\r\n";

   sat->tle.xincl      = t2.inc;
   sat->tle.xnodeo     = t2.ran;
   sat->tle.eo         = t2.ecc;
   sat->tle.omegao     = t2.aop;
   sat->tle.xmo        = t2.man;
   sat->tle.xno        = t2.mmo;
   sat->tle.revnum     = t2.rev;
}


int main( int argc, char *argv[] )
{
   std::string tle1 = "TLE1=1";
   std::string tle2 = "TLE2=2";

   geodetic_t sat_geodetic;
   sat_t sat;

   std::string line;
   //int satnum = 0;

   if( argc < 2 ) {
      std::cerr << "usage: " << argv[0] << " tle.txt [satnum]" << std::endl;
      exit( 1 );
   }

   std::ifstream in( argv[1] );

   while( std::getline( in, line )) {
      std::string pfx = line.substr( 0, 6 );

      if( pfx == tle1 ) {
         handleT1( &sat, (const uint8_t*)line.data() + 5, line.length() - 5 );
      } else if( pfx == tle2 ) {
         handleT2( &sat, (const uint8_t*)line.data() + 5, line.length() - 5 );
         break;
      }
   }

   sat.flags = 0;

   double min = 0.0;

   for( ;; ) {
      std::cout << std::endl;

      SGP4( &sat, min );
      Convert_Sat_State( &sat.pos, &sat.vel );
      Magnitude( &sat.vel );
      sat.velo = sat.vel.w;
    //Calculate_Obs( sat.jul_utc, &sat.pos, &sat.vel, &obs_geodetic, &obs_set );
      Calculate_LatLonAlt( sat.jul_utc, &sat.pos, &sat_geodetic );

      while( sat_geodetic.lon < -pi )
         sat_geodetic.lon += twopi;

      while( sat_geodetic.lon > pi )
         sat_geodetic.lon -= twopi;

      std::cout << "LAT: " << Degrees( sat_geodetic.lat ) << std::endl;
      std::cout << "LON: " << Degrees( sat_geodetic.lon ) << std::endl;
/*
      std::cout << "posX: " << sat.pos.x << std::endl;
      std::cout << "posY: " << sat.pos.y << std::endl;
      std::cout << "posZ: " << sat.pos.z << std::endl;

      std::cout << "velX: " << sat.vel.x << std::endl;
      std::cout << "velY: " << sat.vel.y << std::endl;
      std::cout << "velZ: " << sat.vel.z << std::endl;
*/

      sleep( 1 );

      min += 1.0 / ( 60 * 3600 );
   }

   return 0;
}

/*EoF*/
