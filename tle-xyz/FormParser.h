
#ifndef _QB50_FORM_PARSER_H
#define _QB50_FORM_PARSER_H

#include "Form.h"


namespace qb50 {

   class FormParser
   {
      public:

         static Form *parse ( const uint8_t *x, size_t n );

         static size_t _parseCForm    (  CForm *fp,  const uint8_t *x, size_t n );
         static size_t _parseHForm    (  HForm *fp,  const uint8_t *x, size_t n );
         static size_t _parsePForm    (  PForm *fp,  const uint8_t *x, size_t n );
         static size_t _parseT1Form   ( T1Form *fp,  const uint8_t *x, size_t n );
         static size_t _parseT2Form   ( T2Form *fp,  const uint8_t *x, size_t n );

         static size_t _parseUnsigned ( long     &p, const uint8_t *x, size_t n );
         static size_t _parseInteger  ( long     &p, const uint8_t *x, size_t n );
         static size_t _parseDecimal  ( double   &p, const uint8_t *x, size_t n );
         static size_t _parseDouble   ( double   &p, const uint8_t *x, size_t n );
         static size_t _parseSign     ( int      &p, const uint8_t *x, size_t n );
         static size_t _parseExponent ( double   &p, const uint8_t *x, size_t n );
         static size_t _parseNoradSci ( double   &p, const uint8_t *x, size_t n );
         static size_t _parseComma    (              const uint8_t *x, size_t n );

         static void   _checkTLESum   (              const uint8_t *x, size_t n );
   };

}


#endif /*_QB50_FORM_H*/

/*EoF*/
