
#include <iostream>
#include <cmath>
#include "FormParser.h"

using namespace qb50;

namespace qb50 {
 extern void hexdump( const void *x, unsigned n );
}


static const int _mdays[ 12 ] = {
   0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334
};

#define ISDIGIT( c ) (((c) >= '0') && ((c) <= '9'))
#define ISWHITE( c ) ( (c) == ' '                 )
#define  ISSIGN( c ) (((c) == '+') || ((c) == '-'))


Form *FormParser::parse( const uint8_t *x, size_t n )
{
   Form *fp = (Form*)0;

   if( n < 1 )
      throw( "FormParser::parse()" );

   switch( *x ) {
      case 'C':
      case 'c':

         try {
            fp = new CForm();
            (void)_parseCForm( (CForm*)fp, x + 1, n - 1 );
         } catch( const char *e ) {
            delete fp;
            fp = (Form*)0;
         }

         break;

      case 'H':
      case 'h':

         fp = new HForm();
         (void)_parseHForm( (HForm*)fp, x + 1, n - 1 );

         break;

      case 'P':
      case 'p':

         fp = new PForm();
         (void)_parsePForm( (PForm*)fp, x + 1, n - 1 );

         break;

      case 'T':
      case 't':

         if( n < 2 )
            throw( "FormParser::parse()" );

         switch( x[1] ) {
            case '1':

               fp = new T1Form();
               (void)_parseT1Form( (T1Form*)fp, x + 1, n - 1 );

               break;

            case '2':

               fp = new T2Form();
               (void)_parseT2Form( (T2Form*)fp, x + 1, n - 1 );

               break;

            default:
               ;
         }

         break;

      default:
         ;
   }

   return fp;
}


size_t FormParser::_parseCForm( CForm *fp, const uint8_t *x, size_t n )
{
   size_t off = 0;
   int i;

   for( i = 0 ; i < 4 ; ++i ) {
      off += _parseInteger( fp->argv[i], x + off, n - off );
      if(( off >= n ) || ( x[off] != ',' )) break;
      ++off;
   }

   fp->argc = i + 1;

   return off;
}


size_t FormParser::_parseHForm( HForm *fp, const uint8_t *x, size_t n )
{
   size_t off = 0;
   size_t argc;
   long   argv[6];

   for( argc = 0 ; argc < 6 ; ++argc ) {
      if(( off >= n ) || ( x[off] != ',' )) break;
      ++off;
      off += _parseInteger( argv[argc], x + off, n - off );
   }

   if( argc < 3 )
      throw( "FormParser::_parseHForm()" );

   while( argc < 6 )
      argv[argc++] = 0;

   // see: https://github.com/git/git/blob/master/date.c
   //      https://en.wikipedia.org/wiki/Leap_year

   long dd = argv[0];
   long mm = argv[1] - 1;
   long yy = argv[2] - 1970;

   if(( yy < 0 ) || ( yy > 129 ))
      throw( "FormParser::_parseHForm()" );

   if(( mm < 0 ) || ( mm > 11 ))
      throw( "FormParser::_parseHForm()" );

   if(( mm < 2 ) || (( yy + 2 ) % 4 ))
      --dd;

   if(( argv[3] < 0 ) || ( argv[4] < 0 ) || ( argv[5] < 0 ))
      throw( "FormParser::_parseHForm()" );

   time_t ts = ( yy * 365 + ( yy + 1 ) / 4 + _mdays[mm] + dd ) * 24 * 60 * 60
             + argv[3] * 60 * 60
             + argv[4] * 60
             + argv[5];

   fp->tv.tv_sec  = ts;
   fp->tv.tv_usec = 0;

   return off;
}


size_t FormParser::_parsePForm( PForm *fp, const uint8_t *x, size_t n )
{
   size_t off = 0;

   off += _parseInteger( fp->pnum, x,       n       );
   off +=   _parseComma(           x + off, n - off ); /* XXX bof... */
   off += _parseInteger( fp->pval, x + off, n - off );

   return off;
}


size_t FormParser::_parseT1Form( T1Form *fp, const uint8_t *x, size_t n )
{
   _checkTLESum( x, n );

   (void)_parseInteger ( fp->sat, x + 2,   5 );
   (void)_parseInteger ( fp->eyr, x + 18,  2 );
   (void)_parseDouble  ( fp->edy, x + 20, 12 );
   (void)_parseDouble  ( fp->d1m, x + 33, 10 );
   (void)_parseNoradSci( fp->d2m, x + 44,  8 );
   (void)_parseNoradSci( fp->bdt, x + 53,  8 );

   return 69;
}


size_t FormParser::_parseT2Form( T2Form *fp, const uint8_t *x, size_t n )
{
   _checkTLESum( x, n );

   (void)_parseInteger( fp->sat, x +  2,  5 );
   (void)_parseDouble ( fp->inc, x +  8,  8 );
   (void)_parseDouble ( fp->ran, x + 17,  8 );
   (void)_parseDecimal( fp->ecc, x + 26,  7 );
   (void)_parseDouble ( fp->aop, x + 34,  8 );
   (void)_parseDouble ( fp->man, x + 43,  8 );
   (void)_parseDouble ( fp->mmo, x + 52, 11 );
   (void)_parseInteger( fp->rev, x + 63,  5 );

   return 69;
}


size_t FormParser::_parseInteger( long &p, const uint8_t *x, size_t n )
{
   size_t off = 0;
   int      m = 1; /* multiplier   */
   long     u = 0; /* integer part */

   uint8_t c;

   while(( off < n ) && ISWHITE( x[off] ))
      ++off;

   off += _parseSign( m, x + off, n - off );

   while( off < n ) {
      c = x[off];
      if( !ISDIGIT( c )) break;
      u = u * 10 + ( c - '0' );
      ++off;
   }

   p = u * m;

   return off;
}


size_t FormParser::_parseDecimal( double &p, const uint8_t *x, size_t n )
{
   size_t off = 0;
   long     d = 1; /* divisor      */
   long     u = 0; /* decimal part */

   uint8_t c;

   if(( n < 1 ) || !ISDIGIT( x[0] ))
      throw( "FormParser::_parseDecimal()" ); /* XXX */

   while( off < n ) {
      c = x[off];
      if( !ISDIGIT( c )) break;
      u = u * 10 + ( c - '0' );
      d *= 10;
      ++off;
   }

   p = (double)u / d;

   return off;
}


size_t FormParser::_parseDouble( double &p, const uint8_t *x, size_t n )
{
   size_t off = 0;
   int      m = 1; /* multiplier   */
   long    ip = 0; /* integer part */
   double  dp = 0; /* decimal part */

   while(( off < n ) && ISWHITE( x[off] ))
      ++off;

   off += _parseSign( m, x + off, n - off );

   if(( off < n ) && ISDIGIT( x[off] ))
      off += _parseInteger( ip, x + off, n - off );

   if(( off < n ) && ( x[off] == '.' )) {
      ++off;
      off += _parseDecimal( dp, x + off, n - off );
   }

   p = ( dp + ip ) * m;

   return off;
}


size_t FormParser::_parseSign( int &p, const uint8_t *x, size_t n )
{
   size_t off = 0;

   if( n < 1 )
      throw( "FormParser::_parseSign()" );

   switch( *x ) {
      case '-': p = -1; off = 1; break;
      case '+': p =  1; off = 1; break;
       default: p =  1; off = 0; break;
   }

   return off;
}


size_t FormParser::_parseExponent( double &p, const uint8_t *x, size_t n )
{
   size_t off = 0;
   long     u = 0; /* decimal part */

   if(( n < 1 ) || !ISSIGN( x[0] ))
      throw( "FormParser::_parseExponent()" );

   off += _parseInteger( u, x + off, n - off );

   //p = exp10( u );
   p = pow( 10, u );

   return off;
}


size_t FormParser::_parseNoradSci( double &p, const uint8_t *x, size_t n )
{
   size_t off = 0;
   int      m = 0; /* multiplier   */
   double  dp = 0; /* decimal part */
   double   e = 0; /* exponent     */

   while(( off < n ) && ISWHITE( x[off] ))
      ++off;

   off +=     _parseSign( m,  x + off, n - off );
   off +=  _parseDecimal( dp, x + off, n - off );
   off += _parseExponent( e,  x + off, n - off );

   p = dp * e * m;

   return off;
}


size_t FormParser::_parseComma( const uint8_t *x, size_t n )
{
   if(( n < 1 ) || ( *x != ',' ))
      throw( "FormParser::_parseComma()" );

   return 1;
}


void FormParser::_checkTLESum( const uint8_t *x, size_t n )
{
   size_t  off;
   uint8_t c;
   uint8_t s = 0;

   if(( n < 69 ) || ( !ISDIGIT( x[68] )))
      throw( "FormParser::_checkTLESum()" );

   for( off = 0 ; off < 68 ; ++off ) {
      c = x[off];
      if( ISDIGIT(c) ) {
         s = ( s + ( c - '0' )) % 10;
      } else if( c == '-' ) {
         s = ( s + 1 ) % 10;
      }
   }

   if( s != ( x[68] - '0' ))
      throw( "FormParser::_checkTLESum()" );
}


/*EoF*/
