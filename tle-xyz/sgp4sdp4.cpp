/* -*- Mode: C; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- */
/*
 *  Unit SGP4SDP4
 *           Author:  Dr TS Kelso
 * Original Version:  1991 Oct 30
 * Current Revision:  1992 Sep 03
 *          Version:  1.50
 *        Copyright:  1991-1992, All Rights Reserved
 *
 *   Ported to C by:  Neoklis Kyriazis  April 10  2001
 *   Reentrancy mods by Alexandru Csete OZ9AEC
 */

#include "sgp4sdp4.h"

/* SGP4 */
/* This function is used to calculate the position and velocity */
/* of near-earth (period < 225 minutes) satellites. tsince is   */
/* time since epoch in minutes, tle is a pointer to a tle_t     */
/* structure with Keplerian orbital elements and pos and vel    */
/* are vector_t structures returning ECI satellite position and */
/* velocity. Use Convert_Sat_State() to convert to km and km/s.*/
void
SGP4 (sat_t *sat, double tsince)
{
   double
      cosuk, sinuk, rfdotk, vx, vy, vz, ux, uy, uz, xmy, xmx,
      cosnok, sinnok, cosik, sinik, rdotk, xinck, xnodek, uk,
      rk, cos2u, sin2u, u, sinu, cosu, betal, rfdot, rdot, r, pl,
      elsq, esine, ecose, epw, cosepw, x1m5th, xhdot1, tfour,
      sinepw, capu, ayn, xlt, aynl, xll, axn, xn, beta, xl, e, a,
      tcube, delm, delomg, templ, tempe, tempa, xnode, tsq, xmp,
      omega, xnoddf, omgadf, xmdf, a1, a3ovk2, ao, betao, betao2,
      c1sq, c2, c3, coef, coef1, del1, delo, eeta, eosq, etasq,
      perige, pinvsq, psisq, qoms24, s4, temp, temp1, temp2,
      temp3, temp4, temp5, temp6, theta2, theta4, tsi;

   int i;

   /* Initialization */
   if (~sat->flags & SGP4_INITIALIZED_FLAG) {
   // (!(sat->flags & SGP4_INITIALIZED_FLAG))

      sat->flags |= SGP4_INITIALIZED_FLAG;

      //g_print ("SAT %d INITIALISED.\n", sat->tle.catnr);

      /* Recover original mean motion (xnodp) and   */
      /* semimajor axis (aodp) from input elements. */
      a1 = pow (xke/sat->tle.xno, tothrd);
      sat->sgps.cosio = cos (sat->tle.xincl);
      theta2 = sat->sgps.cosio * sat->sgps.cosio;
      sat->sgps.x3thm1 = 3 * theta2 - 1.0;
      eosq = sat->tle.eo * sat->tle.eo;
      betao2 = 1 - eosq;
      betao = sqrt (betao2);
      del1 = 1.5 * ck2 * sat->sgps.x3thm1 / (a1*a1*betao*betao2);
      ao = a1*(1-del1*(0.5*tothrd+del1*(1+134.0/81.0*del1)));
      delo = 1.5 * ck2 * sat->sgps.x3thm1 / (ao*ao*betao*betao2);
      sat->sgps.xnodp = sat->tle.xno / (1.0 + delo);
      sat->sgps.aodp = ao / (1.0 - delo);

      /* For perigee less than 220 kilometers, the "simple" flag is set */
      /* and the equations are truncated to linear variation in sqrt a  */
      /* and quadratic variation in mean anomaly.  Also, the c3 term,   */
      /* the delta omega term, and the delta m term are dropped.        */
      if ((sat->sgps.aodp * (1.0 - sat->tle.eo) / ae) < (220.0 / xkmper + ae))
         sat->flags |= SIMPLE_FLAG;
      else
         sat->flags &= ~SIMPLE_FLAG;

      /* For perigee below 156 km, the       */
      /* values of s and qoms2t are altered. */
      s4 = __s__;
      qoms24 = qoms2t;
      perige = (sat->sgps.aodp * (1 - sat->tle.eo) - ae) * xkmper;
      if (perige < 156.0) {
         if (perige <= 98.0)
            s4 = 20.0;
         else
            s4 = perige - 78.0;
         qoms24 = pow ((120.0 - s4) * ae / xkmper, 4);
         s4 = s4 / xkmper + ae;
      }; /* FIXME FIXME: End of if(perige <= 98) NO WAY!!!! */

      pinvsq = 1.0 / (sat->sgps.aodp * sat->sgps.aodp * betao2 * betao2);
      tsi = 1.0 / (sat->sgps.aodp - s4);
      sat->sgps.eta = sat->sgps.aodp * sat->tle.eo * tsi;
      etasq = sat->sgps.eta * sat->sgps.eta;
      eeta = sat->tle.eo * sat->sgps.eta;
      psisq = fabs (1.0 - etasq);
      coef = qoms24 * pow (tsi, 4);
      coef1 = coef / pow (psisq, 3.5);
      c2 = coef1 * sat->sgps.xnodp * (sat->sgps.aodp *
                  (1.0 + 1.5 * etasq + eeta * (4.0 + etasq)) +
                  0.75 * ck2 * tsi / psisq * sat->sgps.x3thm1 *
                  (8.0 + 3.0 * etasq * (8 + etasq)));
      sat->sgps.c1 = c2 * sat->tle.bstar;
      sat->sgps.sinio = sin (sat->tle.xincl);
      a3ovk2 = -xj3 / ck2 * pow (ae, 3);
      c3 = coef * tsi * a3ovk2 * sat->sgps.xnodp * ae * sat->sgps.sinio / sat->tle.eo;
      sat->sgps.x1mth2 = 1.0 - theta2;
      sat->sgps.c4 = 2.0 * sat->sgps.xnodp * coef1 * sat->sgps.aodp * betao2 *
         (sat->sgps.eta * (2.0 + 0.5 * etasq) +
          sat->tle.eo * (0.5 + 2.0 * etasq) -
          2.0 * ck2 * tsi / (sat->sgps.aodp * psisq) *
          (-3.0 * sat->sgps.x3thm1 * (1.0 - 2.0 * eeta + etasq * (1.5 - 0.5 * eeta)) +
           0.75 * sat->sgps.x1mth2 * (2.0 * etasq - eeta * (1.0 + etasq)) *
           cos (2.0 * sat->tle.omegao)));
      sat->sgps.c5 = 2.0 * coef1 * sat->sgps.aodp * betao2 *
         (1.0 + 2.75 * (etasq + eeta) + eeta * etasq);
      theta4 = theta2 * theta2;
      temp1 = 3.0 * ck2 * pinvsq * sat->sgps.xnodp;
      temp2 = temp1 * ck2 * pinvsq;
      temp3 = 1.25 * ck4 * pinvsq * pinvsq * sat->sgps.xnodp;
      sat->sgps.xmdot = sat->sgps.xnodp + 0.5 * temp1 * betao * sat->sgps.x3thm1 +
         0.0625 * temp2 * betao * (13.0 - 78.0 * theta2 + 137.0 * theta4);
      x1m5th = 1.0 - 5.0 * theta2;
      sat->sgps.omgdot = -0.5 * temp1 * x1m5th +
         0.0625 * temp2 * (7.0 - 114.0 * theta2 + 395.0 * theta4) +
         temp3 * (3.0 - 36.0 * theta2 + 49.0 * theta4);
      xhdot1 = -temp1 * sat->sgps.cosio;
      sat->sgps.xnodot = xhdot1 + (0.5 * temp2 * (4.0 - 19.0 * theta2) +
                    2.0 * temp3 * (3.0 - 7.0 * theta2)) * sat->sgps.cosio;
      sat->sgps.omgcof = sat->tle.bstar * c3 * cos (sat->tle.omegao);
      sat->sgps.xmcof = -tothrd * coef * sat->tle.bstar * ae / eeta;
      sat->sgps.xnodcf = 3.5 * betao2 * xhdot1 * sat->sgps.c1;
      sat->sgps.t2cof = 1.5 * sat->sgps.c1;
      sat->sgps.xlcof = 0.125 * a3ovk2 * sat->sgps.sinio *
         (3.0 + 5.0 * sat->sgps.cosio) / (1.0 + sat->sgps.cosio);
      sat->sgps.aycof = 0.25 * a3ovk2 * sat->sgps.sinio;
      sat->sgps.delmo = pow (1.0 + sat->sgps.eta * cos (sat->tle.xmo), 3);
      sat->sgps.sinmo = sin (sat->tle.xmo);
      sat->sgps.x7thm1 = 7.0 * theta2 - 1.0;
      if (~sat->flags & SIMPLE_FLAG) {
         c1sq = sat->sgps.c1 * sat->sgps.c1;
         sat->sgps.d2 = 4.0 * sat->sgps.aodp * tsi * c1sq;
         temp = sat->sgps.d2 * tsi * sat->sgps.c1 / 3.0;
         sat->sgps.d3 = (17.0 * sat->sgps.aodp + s4) * temp;
         sat->sgps.d4 = 0.5 * temp * sat->sgps.aodp * tsi *
            (221.0 * sat->sgps.aodp + 31.0 * s4) * sat->sgps.c1;
         sat->sgps.t3cof = sat->sgps.d2 + 2.0 * c1sq;
         sat->sgps.t4cof = 0.25 * (3.0 * sat->sgps.d3 + sat->sgps.c1 *
                    (12.0 * sat->sgps.d2 + 10.0 * c1sq));
         sat->sgps.t5cof = 0.2 * (3.0 * sat->sgps.d4 +
                   12.0 * sat->sgps.c1 * sat->sgps.d3 +
                   6.0 * sat->sgps.d2 * sat->sgps.d2 +
                   15.0 * c1sq * (2.0 * sat->sgps.d2 + c1sq));
      }; /* End of if (isFlagClear(SIMPLE_FLAG)) */
   }; /* End of SGP4() initialization */

   /* Update for secular gravity and atmospheric drag. */
   xmdf = sat->tle.xmo + sat->sgps.xmdot * tsince;
   omgadf = sat->tle.omegao + sat->sgps.omgdot * tsince;
   xnoddf = sat->tle.xnodeo + sat->sgps.xnodot * tsince;
   omega = omgadf;
   xmp = xmdf;
   tsq = tsince*tsince;
   xnode = xnoddf + sat->sgps.xnodcf * tsq;
   tempa = 1.0 - sat->sgps.c1 * tsince;
   tempe = sat->tle.bstar * sat->sgps.c4 * tsince;
   templ = sat->sgps.t2cof * tsq;
   if (~sat->flags & SIMPLE_FLAG) {
      delomg = sat->sgps.omgcof * tsince;
      delm = sat->sgps.xmcof * (pow (1 + sat->sgps.eta * cos (xmdf), 3) - sat->sgps.delmo);
      temp = delomg + delm;
      xmp = xmdf + temp;
      omega = omgadf - temp;
      tcube = tsq * tsince;
      tfour = tsince * tcube;
      tempa = tempa - sat->sgps.d2 * tsq - sat->sgps.d3 * tcube - sat->sgps.d4 * tfour;
      tempe = tempe + sat->tle.bstar * sat->sgps.c5 * (sin (xmp) - sat->sgps.sinmo);
      templ = templ + sat->sgps.t3cof * tcube + tfour *
         (sat->sgps.t4cof + tsince * sat->sgps.t5cof);
   }; /* End of if (isFlagClear(SIMPLE_FLAG)) */

   a = sat->sgps.aodp * pow (tempa, 2);
   e = sat->tle.eo - tempe;
   xl = xmp + omega + xnode + sat->sgps.xnodp * templ;
   beta = sqrt (1.0 - e*e);
   xn = xke / pow (a, 1.5);

   /* Long period periodics */
   axn = e * cos (omega);
   temp = 1.0 / (a * beta * beta);
   xll = temp * sat->sgps.xlcof * axn;
   aynl = temp * sat->sgps.aycof;
   xlt = xl + xll;
   ayn = e * sin (omega) + aynl;

   /* Solve Kepler's' Equation */
   capu = FMod2p (xlt - xnode);
   temp2 = capu;

   i = 0;
   do {
      sinepw = sin (temp2);
      cosepw = cos (temp2);
      temp3 = axn * sinepw;
      temp4 = ayn * cosepw;
      temp5 = axn * cosepw;
      temp6 = ayn * sinepw;
      epw = (capu - temp4 + temp3 - temp2) / (1.0 - temp5 - temp6) + temp2;
      if (fabs (epw - temp2) <= e6a)
         break;
      temp2 = epw;
   }
   while( i++ < 10 );

   /* Short period preliminary quantities */
   ecose = temp5 + temp6;
   esine = temp3 - temp4;
   elsq = axn*axn + ayn*ayn;
   temp = 1.0 - elsq;
   pl = a * temp;
   r = a * (1.0 - ecose);
   temp1 = 1.0 / r;
   rdot = xke * sqrt (a) * esine * temp1;
   rfdot = xke * sqrt (pl) * temp1;
   temp2 = a * temp1;
   betal = sqrt (temp);
   temp3 = 1.0 / (1.0 + betal);
   cosu = temp2 * (cosepw - axn + ayn * esine * temp3);
   sinu = temp2 * (sinepw - ayn - axn * esine * temp3);
   u = AcTan (sinu, cosu);
   sin2u = 2.0 * sinu * cosu;
   cos2u = 2.0 * cosu * cosu - 1.0;
   temp = 1.0 / pl;
   temp1 = ck2 * temp;
   temp2 = temp1 * temp;

   /* Update for short periodics */
   rk = r * (1.0 - 1.5 * temp2 * betal * sat->sgps.x3thm1) +
      0.5 * temp1 * sat->sgps.x1mth2 * cos2u;
   uk = u - 0.25 * temp2 * sat->sgps.x7thm1 * sin2u;
   xnodek = xnode + 1.5 * temp2 * sat->sgps.cosio * sin2u;
   xinck = sat->tle.xincl + 1.5 * temp2 * sat->sgps.cosio * sat->sgps.sinio * cos2u;
   rdotk = rdot - xn * temp1 * sat->sgps.x1mth2 * sin2u;
   rfdotk = rfdot + xn * temp1 * (sat->sgps.x1mth2 * cos2u + 1.5 * sat->sgps.x3thm1);


   /* Orientation vectors */
   sinuk = sin (uk);
   cosuk = cos (uk);
   sinik = sin (xinck);
   cosik = cos (xinck);
   sinnok = sin (xnodek);
   cosnok = cos (xnodek);
   xmx = -sinnok * cosik;
   xmy = cosnok * cosik;
   ux = xmx * sinuk + cosnok * cosuk;
   uy = xmy * sinuk + sinnok * cosuk;
   uz = sinik * sinuk;
   vx = xmx * cosuk - cosnok * sinuk;
   vy = xmy * cosuk - sinnok * sinuk;
   vz = sinik * cosuk;

   /* Position and velocity */
   sat->pos.x = rk*ux;
   sat->pos.y = rk*uy;
   sat->pos.z = rk*uz;
   sat->vel.x = rdotk*ux+rfdotk*vx;
   sat->vel.y = rdotk*uy+rfdotk*vy;
   sat->vel.z = rdotk*uz+rfdotk*vz;

   sat->phase = xlt - xnode - omgadf + twopi;
   if (sat->phase < 0)
      sat->phase += twopi;
   sat->phase = FMod2p (sat->phase);

   sat->tle.omegao1 = omega;
   sat->tle.xincl1  = xinck;
   sat->tle.xnodeo1 = xnodek;

} /*SGP4*/

static int Flags = 0;

void
SetFlag(int flag)
{
   Flags |= flag;
}

void
ClearFlag(int flag)
{
   Flags &= ~flag;
}

