
#ifndef _QB50_FORM_H
#define _QB50_FORM_H

#include <cstdint>
#include <ctime>

extern "C" {
#include <sys/time.h>
}

namespace qb50 {

   struct Form
   {
      enum FormType {
         FORM_TYPE_0  = 0,
         FORM_TYPE_C  = 1,
         FORM_TYPE_H  = 2,
         FORM_TYPE_P  = 3,
         FORM_TYPE_T  = 4,
         FORM_TYPE_T1 = 5,
         FORM_TYPE_T2 = 6
      };

      FormType type;

      protected:

         Form( FormType ftype ) : type( ftype )
         { ; }
   };


   struct CForm : public Form
   {
      int  argc;
      long argv[4];

      CForm() : Form( FORM_TYPE_C )
      { argv[0] = argv[1] = argv[2] = argv[3] = 0; }
   };


   struct HForm : public Form
   {
      ::timeval tv;

      HForm() : Form( FORM_TYPE_H )
      { ; }
   };


   struct PForm : public Form
   {
      long pnum;
      long pval;

      PForm() : Form( FORM_TYPE_P )
      { ; }
   };


   struct T1Form : public Form
   {
      long   sat;  /* satellite number                    */
      long   eyr;  /* epoch year                          */
      double edy;  /* epoch day of the year               */
      double d1m;  /* 1st. time derivative of mean motion */
      double d2m;  /* 2nd. time derivative of mean motion */
      double bdt;  /* BSTAR drag term                     */

      T1Form() : Form( FORM_TYPE_T1 )
      { ; }
   };


   struct T2Form : public Form
   {
      long   sat;  /* satellite number                    */
      long   rev;  /* revolution number                   */
      double inc;  /* inclination (degrees)               */
      double ran;  /* R.A.A.N (degrees)                   */
      double ecc;  /* eccentricity                        */
      double aop;  /* argument of perigee (degrees)       */
      double man;  /* mean anomaly (degrees)              */
      double mmo;  /* mean motion (revolutions per day)   */

      T2Form() : Form( FORM_TYPE_T2 )
      { ; }
   };

}


#endif /*_QB50_FORM_H*/

/*EoF*/
