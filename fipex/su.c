
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <time.h>

/* 01.01.2000 00:00:00 UTC in UNIX time */
#define VKI_EPOCH 946684800

#define SAFE_READ( x, n ) \
   if( read( 0, x, n ) != n ) die( "can't read from stdin" )


static int die( char *msg )
{
   (void)fprintf( stderr, "error: %s\n", msg );
   exit( 1 );
}


/* read unsigned int 32 little endian */
static uint32_t readu32le( void )
{
   uint8_t x[ 4 ];
   SAFE_READ( x, 4 );
   return x[0] | ( x[1] << 8 ) | ( x[2] << 16 ) | ( x[3] << 24 );
}


/* read unsigned int 16 little endian */
static uint16_t readu16le( void )
{
   uint8_t x[ 2 ];
   SAFE_READ( x, 2 );
   return x[0] | ( x[1] << 8 );
}


int main( int argc, char *argv[] )
{
   uint8_t x[ 16 ];
   uint32_t u32;
   uint16_t u16;
   uint8_t  u8;

   int i, j;

   (void)argc;
   (void)argv;

   SAFE_READ( &u8, 1 );
   (void)printf( "    script length: %d bytes\n", u8 );

   u32 = readu32le();
   time_t dt = VKI_EPOCH + u32;
   struct tm *utc = gmtime( &dt );
   (void)printf( " start time (UTC): %s", asctime( utc ));

   u16 = readu16le();
   (void)printf( "      repeat time: %ds\n", u16 );

   SAFE_READ( &u8, 1 );
   (void)printf( "    command count: %d\n", u8 );

   int cc = u8;
   for( i = 0 ; i < cc ; ++i ) {
      uint8_t xor = 0x00;

      SAFE_READ( &u8, 1 );
      if( u8 != 0x7e ) die( "garbage on stdin" );

      SAFE_READ( &u8, 1 );
      (void)printf( "CMD_ID: 0x%02x", u8 );
      if( u8 == 0xff ) {
         (void)printf( "\n" );
         break;
      }
      xor ^= u8;

      SAFE_READ( &u8, 1 );
      (void)printf( ", LEN: %d", u8 );
      xor ^= u8;

      if( u8 > 0 ) {
         SAFE_READ( x, u8 );
         xor ^= x[0];
         (void)printf( ", DATA: [ 0x%02x", x[0] );
         for( j = 1 ; j < u8 ; ++j ) {
            (void)printf( ", 0x%02x", x[j] );
            xor ^= x[j];
         }
         (void)printf( " ]" );
      }

      SAFE_READ( &u8, 1 );
      (void)printf( ", XOR: 0x%02x %s", u8, ( u8 == xor ) ? "(match)" : "MISMATCH" );

      u16 = readu16le();
      (void)printf( ", DELAY:" );
      if( u16 == 0xffff ) {
         (void)printf( " none\n" );
      } else {
         (void)printf( " %ds\n", u16 );
      }
   }

   return 0;
}

