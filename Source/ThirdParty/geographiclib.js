/*
 * Geodesic routines from GeographicLib translated to JavaScript.  See
 * http://geographiclib.sf.net/html/other.html#javascript
 *
 * The algorithms are derived in
 *
 *    Charles F. F. Karney,
 *    Algorithms for geodesics, J. Geodesy 87, 43-55 (2013),
 *    https://dx.doi.org/10.1007/s00190-012-0578-z
 *    Addenda: http://geographiclib.sf.net/geod-addenda.html
 *
 * This file is the concatenation and compression of the JavaScript files in
 * doc/scripts/GeographicLib in the source tree for GeographicLib.
 *
 * Copyright (c) Charles Karney (2011-2015) <charles@karney.com> and licensed
 * under the MIT/X11 License.  For more information, see
 * http://geographiclib.sf.net/
 *
 * Version: 1.45
 * File inventory:
 *   Math.js Geodesic.js GeodesicLine.js PolygonArea.js DMS.js
 */

(function(cb) {

/**************** Math.js ****************/
/*
 * Math.js
 * Transcription of Math.hpp, Constants.hpp, and Accumulator.hpp into
 * JavaScript.
 *
 * Copyright (c) Charles Karney (2011-2015) <charles@karney.com> and licensed
 * under the MIT/X11 License.  For more information, see
 * http://geographiclib.sourceforge.net/
 */

/**
 * @namespace GeographicLib
 * @description The parent namespace for the following modules:
 * - {@link module:GeographicLib/Geodesic GeographicLib/Geodesic} The main
 *   engine for solving geodesic problems via the
 *   {@link module:GeographicLib/Geodesic.Geodesic Geodesic} class.
 *
 * The following modules are used internally by the package:
 * - {@link module:GeographicLib/Math GeographicLib/Math} defines various
 *   mathematical functions.
 */
var GeographicLib = {};
GeographicLib.Math = {};

(function(
  /**
   * @exports GeographicLib/Math
   * @description Some useful mathematical constants and functions (mainly for
   *   internal use).
   */
  m) {
  "use strict";

  /**
   * @summary The number of digits of precision in floating-point numbers.
   * @constant {number}
   */
  m.digits = 53;
  /**
   * @summary The machine epsilon.
   * @constant {number}
   */
  m.epsilon = Math.pow(0.5, m.digits - 1);
  /**
   * @summary The factor to convert degrees to radians.
   * @constant {number}
   */
  m.degree = Math.PI/180;

  /**
   * @summary Square a number.
   * @param {number} x the number.
   * @returns {number} the square.
   */
  m.sq = function(x) { return x * x; };

  /**
   * @summary The hypotenuse function.
   * @param {number} x the first side.
   * @param {number} y the second side.
   * @returns {number} the hypotenuse.
   */
  m.hypot = function(x, y) {
    var a, b;
    x = Math.abs(x);
    y = Math.abs(y);
    a = Math.max(x, y); b = Math.min(x, y) / (a ? a : 1);
    return a * Math.sqrt(1 + b * b);
  };

  /**
   * @summary Cube root function.
   * @param {number} x the argument.
   * @returns {number} the real cube root.
   */
  m.cbrt = function(x) {
    var y = Math.pow(Math.abs(x), 1/3);
    return x < 0 ? -y : y;
  };

  /**
   * @summary The log1p function.
   * @param {number} x the argument.
   * @returns {number} log(1 + x).
   */
  m.log1p = function(x) {
    var y = 1 + x,
        z = y - 1;
    // Here's the explanation for this magic: y = 1 + z, exactly, and z
    // approx x, thus log(y)/z (which is nearly constant near z = 0) returns
    // a good approximation to the true log(1 + x)/x.  The multiplication x *
    // (log(y)/z) introduces little additional error.
    return z === 0 ? x : x * Math.log(y) / z;
  };

  /**
   * @summary Inverse hyperbolic tangent.
   * @param {number} x the argument.
   * @returns {number} tanh<sup>&minus;1</sup> x.
   */
  m.atanh = function(x) {
    var y = Math.abs(x);          // Enforce odd parity
    y = m.log1p(2 * y/(1 - y))/2;
    return x < 0 ? -y : y;
  };

  /**
   * @summary An error-free sum.
   * @param {number} u
   * @param {number} v
   * @returns {object} sum with sum.s = round(u + v) and sum.t is u + v &minus;
   *   round(u + v)
   */
  m.sum = function(u, v) {
    var s = u + v,
        up = s - v,
        vpp = s - up,
        t;
    up -= u;
    vpp -= v;
    t = -(up + vpp);
    // u + v =       s      + t
    //       = round(u + v) + t
    return {s: s, t: t};
  };

  /**
   * @summary Evaluate a polynomial.
   * @param {integer} N the order of the polynomial.
   * @param {array} p the coefficient array (of size N + 1) (leading
   *   order coefficient first)
   * @param {number} x the variable.
   * @returns {number} the value of the polynomial.
   */
  m.polyval = function(N, p, s, x) {
    var y = N < 0 ? 0 : p[s++];
    while (--N >= 0) y = y * x + p[s++];
    return y;
  };

  /**
   * @summary Coarsen a value close to zero.
   * @param {number} x
   * @returns {number} the coarsened value.
   */
  m.AngRound = function(x) {
    // The makes the smallest gap in x = 1/16 - nextafter(1/16, 0) = 1/2^57 for
    // reals = 0.7 pm on the earth if x is an angle in degrees.  (This is about
    // 1000 times more resolution than we get with angles around 90 degrees.)
    // We use this to avoid having to deal with near singular cases when x is
    // non-zero but tiny (e.g., 1.0e-200).  This also converts -0 to +0.
    var z = 1/16,
        y = Math.abs(x);
    // The compiler mustn't "simplify" z - (z - y) to y
    y = y < z ? z - (z - y) : y;
    return x < 0 ? 0 - y : y;
  };

  /**
   * @summary Normalize an angle.
   * @param {number} x the angle in degrees.
   * @returns {number} the angle reduced to the range [&minus;180&deg;,
   *   180&deg;).
   */
  m.AngNormalize = function(x) {
    // Place angle in [-180, 180).
    x = x % 360;
    return x < -180 ? x + 360 : (x < 180 ? x : x - 360);
  };

  /**
   * @summary Normalize a latitude.
   * @param {number} x the angle in degrees.
   * @returns {number} x if it is in the range [&minus;90&deg;, 90&deg;],
   *   otherwise return NaN.
   */
  m.LatFix = function(x) {
    // Replace angle with NaN if outside [-90, 90].
    return Math.abs(x) > 90 ? Number.NaN : x;
  };

  /**
   * @summary Difference of two angles reduced to [&minus;180&deg;,
   *   180&deg;]
   * @param {number} x the first angle in degrees.
   * @param {number} y the second angle in degrees.
   * @return {number} y &minus; x, reduced to the range [&minus;180&deg;,
   *   180&deg;].
   */
  m.AngDiff = function(x, y) {
    // Compute y - x and reduce to [-180,180] accurately.
    var r = m.sum(m.AngNormalize(x), m.AngNormalize(-y)),
        d = - m.AngNormalize(r.s),
        t = r.t;
    return (d == 180 && t < 0 ? -180 : d) - t;
  };

  /**
   * @summary Evaluate the sine and cosine function with the argument in
   *   degrees
   * @param {number} x in degrees.
   * @returns {object} r with r.s = sin(x) and r.c = cos(x).
   */
  m.sincosd = function(x) {
    // In order to minimize round-off errors, this function exactly reduces
    // the argument to the range [-45, 45] before converting it to radians.
    var r, q, s, c, sinx, cosx;
    r = x % 360;
    q = Math.floor(r / 90 + 0.5);
    r -= 90 * q;
    // now abs(r) <= 45
    r *= this.degree;
    // Possibly could call the gnu extension sincos
    s = Math.sin(r); c = Math.cos(r);
    switch (q & 3) {
    case  0: sinx =     s; cosx =     c; break;
    case  1: sinx =     c; cosx = 0 - s; break;
    case  2: sinx = 0 - s; cosx = 0 - c; break;
    default: sinx = 0 - c; cosx =     s; break; // case 3
    }
    return {s: sinx, c: cosx};
  };

  /**
   * @summary Evaluate the atan2 function with the result in degrees
   * @param {number} y
   * @param {number} x
   * @returns atan2(y, x) in degrees, in the range [&minus;180&deg;
   *   180&deg;).
   */
  m.atan2d = function(y, x) {
    // In order to minimize round-off errors, this function rearranges the
    // arguments so that result of atan2 is in the range [-pi/4, pi/4] before
    // converting it to degrees and mapping the result to the correct
    // quadrant.
    var q = 0, t, ang;
    if (Math.abs(y) > Math.abs(x)) { t = x; x = y; y = t; q = 2; }
    if (x < 0) { x = -x; ++q; }
    // here x >= 0 and x >= abs(y), so angle is in [-pi/4, pi/4]
    ang = Math.atan2(y, x) / this.degree;
    switch (q) {
      // Note that atan2d(-0.0, 1.0) will return -0.  However, we expect that
      // atan2d will not be called with y = -0.  If need be, include
      //
      //   case 0: ang = 0 + ang; break;
      //
      // and handle mpfr as in AngRound.
    case 1: ang = (y > 0 ? 180 : -180) - ang; break;
    case 2: ang =  90 - ang; break;
    case 3: ang = -90 + ang; break;
    }
    return ang;
  };
})(GeographicLib.Math);

/**************** Geodesic.js ****************/
/*
 * Geodesic.js
 * Transcription of Geodesic.[ch]pp into JavaScript.
 *
 * See the documentation for the C++ class.  The conversion is a literal
 * conversion from C++.
 *
 * The algorithms are derived in
 *
 *    Charles F. F. Karney,
 *    Algorithms for geodesics, J. Geodesy 87, 43-55 (2013);
 *    https://dx.doi.org/10.1007/s00190-012-0578-z
 *    Addenda: http://geographiclib.sf.net/geod-addenda.html
 *
 * Copyright (c) Charles Karney (2011-2015) <charles@karney.com> and licensed
 * under the MIT/X11 License.  For more information, see
 * http://geographiclib.sourceforge.net/
 */

// Load AFTER Math.js

GeographicLib.Geodesic = {};

(function(
  /**
   * @exports GeographicLib/Geodesic
   * @description Solve geodesic problems via the
   *   {@link module:GeographicLib/Geodesic.Geodesic Geodesic} class.
   */
  g, l, p, m, c) {
  "use strict";

  var GEOGRAPHICLIB_GEODESIC_ORDER = 6,
      nA1_ = GEOGRAPHICLIB_GEODESIC_ORDER,
      nA2_ = GEOGRAPHICLIB_GEODESIC_ORDER,
      nA3_ = GEOGRAPHICLIB_GEODESIC_ORDER,
      nA3x_ = nA3_,
      nC3x_, nC4x_,
      maxit1_ = 20,
      maxit2_ = maxit1_ + m.digits + 10,
      tol0_ = m.epsilon,
      tol1_ = 200 * tol0_,
      tol2_ = Math.sqrt(tol0_),
      tolb_ = tol0_ * tol1_,
      xthresh_ = 1000 * tol2_,
      CAP_NONE = 0,
      CAP_ALL  = 0x1F,
      CAP_MASK = CAP_ALL,
      OUT_ALL  = 0x7F80,
      Astroid,
      A1m1f_coeff, C1f_coeff, C1pf_coeff,
      A2m1f_coeff, C2f_coeff,
      A3_coeff, C3_coeff, C4_coeff;

  g.tiny_ = Math.sqrt(Number.MIN_VALUE);
  g.nC1_ = GEOGRAPHICLIB_GEODESIC_ORDER;
  g.nC1p_ = GEOGRAPHICLIB_GEODESIC_ORDER;
  g.nC2_ = GEOGRAPHICLIB_GEODESIC_ORDER;
  g.nC3_ = GEOGRAPHICLIB_GEODESIC_ORDER;
  g.nC4_ = GEOGRAPHICLIB_GEODESIC_ORDER;
  nC3x_ = (g.nC3_ * (g.nC3_ - 1)) / 2;
  nC4x_ = (g.nC4_ * (g.nC4_ + 1)) / 2,
  g.CAP_C1   = 1<<0;
  g.CAP_C1p  = 1<<1;
  g.CAP_C2   = 1<<2;
  g.CAP_C3   = 1<<3;
  g.CAP_C4   = 1<<4;

  g.SinCosSeries = function(sinp, sinx, cosx, c) {
    // Evaluate
    // y = sinp ? sum(c[i] * sin( 2*i    * x), i, 1, n) :
    //            sum(c[i] * cos((2*i+1) * x), i, 0, n-1)
    // using Clenshaw summation.  N.B. c[0] is unused for sin series
    // Approx operation count = (n + 5) mult and (2 * n + 2) add
    var k = c.length,           // Point to one beyond last element
        n = k - (sinp ? 1 : 0),
        ar = 2 * (cosx - sinx) * (cosx + sinx), // 2 * cos(2 * x)
        y0 = n & 1 ? c[--k] : 0, y1 = 0;        // accumulators for sum
    // Now n is even
    n = Math.floor(n/2);
    while (n--) {
      // Unroll loop x 2, so accumulators return to their original role
      y1 = ar * y0 - y1 + c[--k];
      y0 = ar * y1 - y0 + c[--k];
    }
    return (sinp ? 2 * sinx * cosx * y0 : // sin(2 * x) * y0
            cosx * (y0 - y1));            // cos(x) * (y0 - y1)
  };

  Astroid = function(x, y) {
    // Solve k^4+2*k^3-(x^2+y^2-1)*k^2-2*y^2*k-y^2 = 0 for positive
    // root k.  This solution is adapted from Geocentric::Reverse.
    var k,
        p = m.sq(x),
        q = m.sq(y),
        r = (p + q - 1) / 6,
        S, r2, r3, disc, u, T3, T, ang, v, uv, w;
    if ( !(q === 0 && r <= 0) ) {
      // Avoid possible division by zero when r = 0 by multiplying
      // equations for s and t by r^3 and r, resp.
      S = p * q / 4;            // S = r^3 * s
      r2 = m.sq(r);
      r3 = r * r2
      // The discriminant of the quadratic equation for T3.  This is
      // zero on the evolute curve p^(1/3)+q^(1/3) = 1
      disc = S * (S + 2 * r3);
      u = r;
      if (disc >= 0) {
        T3 = S + r3;
        // Pick the sign on the sqrt to maximize abs(T3).  This
        // minimizes loss of precision due to cancellation.  The
        // result is unchanged because of the way the T is used
        // in definition of u.
        T3 += T3 < 0 ? -Math.sqrt(disc) : Math.sqrt(disc);    // T3 = (r * t)^3
        // N.B. cbrt always returns the real root.  cbrt(-8) = -2.
        T = m.cbrt(T3);     // T = r * t
        // T can be zero; but then r2 / T -> 0.
        u += T + (T !== 0 ? r2 / T : 0);
      } else {
        // T is complex, but the way u is defined the result is real.
        ang = Math.atan2(Math.sqrt(-disc), -(S + r3));
        // There are three possible cube roots.  We choose the
        // root which avoids cancellation.  Note that disc < 0
        // implies that r < 0.
        u += 2 * r * Math.cos(ang / 3);
      }
      v = Math.sqrt(m.sq(u) + q);       // guaranteed positive
      // Avoid loss of accuracy when u < 0.
      uv = u < 0 ? q / (v - u) : u + v; // u+v, guaranteed positive
      w = (uv - q) / (2 * v);           // positive?
      // Rearrange expression for k to avoid loss of accuracy due to
      // subtraction.  Division by 0 not possible because uv > 0, w >= 0.
      k = uv / (Math.sqrt(uv + m.sq(w)) + w); // guaranteed positive
    } else {                                  // q == 0 && r <= 0
      // y = 0 with |x| <= 1.  Handle this case directly.
      // for y small, positive root is k = abs(y)/sqrt(1-x^2)
      k = 0;
    }
    return k;
  };

  A1m1f_coeff = [
    // (1-eps)*A1-1, polynomial in eps2 of order 3
      +1, 4, 64, 0, 256,
  ];

  // The scale factor A1-1 = mean value of (d/dsigma)I1 - 1
  g.A1m1f = function(eps) {
    var p = Math.floor(nA1_/2),
        t = m.polyval(p, A1m1f_coeff, 0, m.sq(eps)) / A1m1f_coeff[p + 1];
    return (t + eps) / (1 - eps);
  };

  C1f_coeff = [
    // C1[1]/eps^1, polynomial in eps2 of order 2
      -1, 6, -16, 32,
    // C1[2]/eps^2, polynomial in eps2 of order 2
      -9, 64, -128, 2048,
    // C1[3]/eps^3, polynomial in eps2 of order 1
      +9, -16, 768,
    // C1[4]/eps^4, polynomial in eps2 of order 1
      +3, -5, 512,
    // C1[5]/eps^5, polynomial in eps2 of order 0
      -7, 1280,
    // C1[6]/eps^6, polynomial in eps2 of order 0
      -7, 2048,
  ];

  // The coefficients C1[l] in the Fourier expansion of B1
  g.C1f = function(eps, c) {
    var eps2 = m.sq(eps),
        d = eps,
        o = 0,
        l, p;
    for (l = 1; l <= g.nC1_; ++l) {     // l is index of C1p[l]
      p = Math.floor((g.nC1_ - l) / 2); // order of polynomial in eps^2
      c[l] = d * m.polyval(p, C1f_coeff, o, eps2) / C1f_coeff[o + p + 1];
      o += p + 2;
      d *= eps;
    }
  };

  C1pf_coeff = [
    // C1p[1]/eps^1, polynomial in eps2 of order 2
      +205, -432, 768, 1536,
    // C1p[2]/eps^2, polynomial in eps2 of order 2
      +4005, -4736, 3840, 12288,
    // C1p[3]/eps^3, polynomial in eps2 of order 1
      -225, 116, 384,
    // C1p[4]/eps^4, polynomial in eps2 of order 1
      -7173, 2695, 7680,
    // C1p[5]/eps^5, polynomial in eps2 of order 0
      +3467, 7680,
    // C1p[6]/eps^6, polynomial in eps2 of order 0
      +38081, 61440,
  ];

  // The coefficients C1p[l] in the Fourier expansion of B1p
  g.C1pf = function(eps, c) {
    var eps2 = m.sq(eps),
        d = eps,
        o = 0,
        l, p;
    for (l = 1; l <= g.nC1p_; ++l) {     // l is index of C1p[l]
      p = Math.floor((g.nC1p_ - l) / 2); // order of polynomial in eps^2
      c[l] = d * m.polyval(p, C1pf_coeff, o, eps2) / C1pf_coeff[o + p + 1];
      o += p + 2;
      d *= eps;
    }
  };

  A2m1f_coeff = [
    // (eps+1)*A2-1, polynomial in eps2 of order 3
      -11, -28, -192, 0, 256,
  ];

  // The scale factor A2-1 = mean value of (d/dsigma)I2 - 1
  g.A2m1f = function(eps) {
    var p = Math.floor(nA2_/2),
        t = m.polyval(p, A2m1f_coeff, 0, m.sq(eps)) / A2m1f_coeff[p + 1];
    return (t - eps) / (1 + eps);
  };

  C2f_coeff = [
    // C2[1]/eps^1, polynomial in eps2 of order 2
      +1, 2, 16, 32,
    // C2[2]/eps^2, polynomial in eps2 of order 2
      +35, 64, 384, 2048,
    // C2[3]/eps^3, polynomial in eps2 of order 1
      +15, 80, 768,
    // C2[4]/eps^4, polynomial in eps2 of order 1
      +7, 35, 512,
    // C2[5]/eps^5, polynomial in eps2 of order 0
      +63, 1280,
    // C2[6]/eps^6, polynomial in eps2 of order 0
      +77, 2048,
  ];

  // The coefficients C2[l] in the Fourier expansion of B2
  g.C2f = function(eps, c) {
    var eps2 = m.sq(eps),
        d = eps,
        o = 0,
        l, p;
    for (l = 1; l <= g.nC2_; ++l) {     // l is index of C2[l]
      p = Math.floor((g.nC2_ - l) / 2); // order of polynomial in eps^2
      c[l] = d * m.polyval(p, C2f_coeff, o, eps2) / C2f_coeff[o + p + 1];
      o += p + 2;
      d *= eps;
    }
  };

  /**
   * @class
   * @property {number} a the equatorial radius (meters).
   * @property {number} f the flattening.
   * @summary Initialize a Geodesic object for a specific ellipsoid.
   * @classdesc Performs geodesic calculations on an ellipsoid of revolution.
   *   The routines for solving the direct and inverse problems return an
   *   object with some of the following fields set: lat1, lon1, azi1, lat2,
   *   lon2, azi2, s12, a12, m12, M12, M21, S12.  See {@tutorial 2-interface},
   *   "The results".
   * @example
   * var GeographicLib = require("geographiclib"),
   *     geod = GeographicLib.Geodesic.WGS84;
   * var inv = geod.Inverse(1,2,3,4);
   * console.log("lat1 = " + inv.lat1 + ", lon1 = " + inv.lon1 +
   *             ", lat2 = " + inv.lat2 + ", lon2 = " + inv.lon2 +
   *             ",\nazi1 = " + inv.azi1 + ", azi2 = " + inv.azi2 +
   *             ", s12 = " + inv.s12);
   * @param {number} a the equatorial radius of the ellipsoid (meters).
   * @param {number} f the flattening of the ellipsoid.  Setting f = 0 gives
   *   a sphere (on which geodesics are great circles).  Negative f gives a
   *   prolate ellipsoid.
   * @throws an error if the parameters are illegal.
   */
  g.Geodesic = function(a, f) {
    this.a = a;
    this.f = f;
    this._f1 = 1 - this.f;
    this._e2 = this.f * (2 - this.f);
    this._ep2 = this._e2 / m.sq(this._f1); // e2 / (1 - e2)
    this._n = this.f / ( 2 - this.f);
    this._b = this.a * this._f1;
    // authalic radius squared
    this._c2 = (m.sq(this.a) + m.sq(this._b) *
                (this._e2 === 0 ? 1 :
                 (this._e2 > 0 ? m.atanh(Math.sqrt(this._e2)) :
                  Math.atan(Math.sqrt(-this._e2))) /
                 Math.sqrt(Math.abs(this._e2))))/2;
    // The sig12 threshold for "really short".  Using the auxiliary sphere
    // solution with dnm computed at (bet1 + bet2) / 2, the relative error in
    // the azimuth consistency check is sig12^2 * abs(f) * min(1, 1-f/2) / 2.
    // (Error measured for 1/100 < b/a < 100 and abs(f) >= 1/1000.  For a given
    // f and sig12, the max error occurs for lines near the pole.  If the old
    // rule for computing dnm = (dn1 + dn2)/2 is used, then the error increases
    // by a factor of 2.)  Setting this equal to epsilon gives sig12 = etol2.
    // Here 0.1 is a safety factor (error decreased by 100) and max(0.001,
    // abs(f)) stops etol2 getting too large in the nearly spherical case.
    this._etol2 = 0.1 * tol2_ /
      Math.sqrt( Math.max(0.001, Math.abs(this.f)) *
                 Math.min(1.0, 1 - this.f/2) / 2 );
    if (!(isFinite(this.a) && this.a > 0))
      throw new Error("Major radius is not positive");
    if (!(isFinite(this._b) && this._b > 0))
      throw new Error("Minor radius is not positive");
    this._A3x = new Array(nA3x_);
    this._C3x = new Array(nC3x_);
    this._C4x = new Array(nC4x_);
    this.A3coeff();
    this.C3coeff();
    this.C4coeff();
  };

  A3_coeff = [
    // A3, coeff of eps^5, polynomial in n of order 0
      -3, 128,
    // A3, coeff of eps^4, polynomial in n of order 1
      -2, -3, 64,
    // A3, coeff of eps^3, polynomial in n of order 2
      -1, -3, -1, 16,
    // A3, coeff of eps^2, polynomial in n of order 2
      +3, -1, -2, 8,
    // A3, coeff of eps^1, polynomial in n of order 1
      +1, -1, 2,
    // A3, coeff of eps^0, polynomial in n of order 0
      +1, 1,
  ];

  // The scale factor A3 = mean value of (d/dsigma)I3
  g.Geodesic.prototype.A3coeff = function() {
    var o = 0, k = 0,
        j, p;
    for (j = nA3_ - 1; j >= 0; --j) { // coeff of eps^j
      p = Math.min(nA3_ - j - 1, j);  // order of polynomial in n
      this._A3x[k++] = m.polyval(p, A3_coeff, o, this._n)
        / A3_coeff[o + p + 1];
      o += p + 2;
    }
  };

  C3_coeff = [
    // C3[1], coeff of eps^5, polynomial in n of order 0
      +3, 128,
    // C3[1], coeff of eps^4, polynomial in n of order 1
      +2, 5, 128,
    // C3[1], coeff of eps^3, polynomial in n of order 2
      -1, 3, 3, 64,
    // C3[1], coeff of eps^2, polynomial in n of order 2
      -1, 0, 1, 8,
    // C3[1], coeff of eps^1, polynomial in n of order 1
      -1, 1, 4,
    // C3[2], coeff of eps^5, polynomial in n of order 0
      +5, 256,
    // C3[2], coeff of eps^4, polynomial in n of order 1
      +1, 3, 128,
    // C3[2], coeff of eps^3, polynomial in n of order 2
      -3, -2, 3, 64,
    // C3[2], coeff of eps^2, polynomial in n of order 2
      +1, -3, 2, 32,
    // C3[3], coeff of eps^5, polynomial in n of order 0
      +7, 512,
    // C3[3], coeff of eps^4, polynomial in n of order 1
      -10, 9, 384,
    // C3[3], coeff of eps^3, polynomial in n of order 2
      +5, -9, 5, 192,
    // C3[4], coeff of eps^5, polynomial in n of order 0
      +7, 512,
    // C3[4], coeff of eps^4, polynomial in n of order 1
      -14, 7, 512,
    // C3[5], coeff of eps^5, polynomial in n of order 0
      +21, 2560,
  ];

  // The coefficients C3[l] in the Fourier expansion of B3
  g.Geodesic.prototype.C3coeff = function() {
    var o = 0, k = 0,
        l, j, p;
    for (l = 1; l < g.nC3_; ++l) {        // l is index of C3[l]
      for (j = g.nC3_ - 1; j >= l; --j) { // coeff of eps^j
        p = Math.min(g.nC3_ - j - 1, j);  // order of polynomial in n
        this._C3x[k++] = m.polyval(p, C3_coeff, o, this._n) /
          C3_coeff[o + p + 1];
        o += p + 2;
      }
    }
  };

  C4_coeff = [
    // C4[0], coeff of eps^5, polynomial in n of order 0
      +97, 15015,
    // C4[0], coeff of eps^4, polynomial in n of order 1
      +1088, 156, 45045,
    // C4[0], coeff of eps^3, polynomial in n of order 2
      -224, -4784, 1573, 45045,
    // C4[0], coeff of eps^2, polynomial in n of order 3
      -10656, 14144, -4576, -858, 45045,
    // C4[0], coeff of eps^1, polynomial in n of order 4
      +64, 624, -4576, 6864, -3003, 15015,
    // C4[0], coeff of eps^0, polynomial in n of order 5
      +100, 208, 572, 3432, -12012, 30030, 45045,
    // C4[1], coeff of eps^5, polynomial in n of order 0
      +1, 9009,
    // C4[1], coeff of eps^4, polynomial in n of order 1
      -2944, 468, 135135,
    // C4[1], coeff of eps^3, polynomial in n of order 2
      +5792, 1040, -1287, 135135,
    // C4[1], coeff of eps^2, polynomial in n of order 3
      +5952, -11648, 9152, -2574, 135135,
    // C4[1], coeff of eps^1, polynomial in n of order 4
      -64, -624, 4576, -6864, 3003, 135135,
    // C4[2], coeff of eps^5, polynomial in n of order 0
      +8, 10725,
    // C4[2], coeff of eps^4, polynomial in n of order 1
      +1856, -936, 225225,
    // C4[2], coeff of eps^3, polynomial in n of order 2
      -8448, 4992, -1144, 225225,
    // C4[2], coeff of eps^2, polynomial in n of order 3
      -1440, 4160, -4576, 1716, 225225,
    // C4[3], coeff of eps^5, polynomial in n of order 0
      -136, 63063,
    // C4[3], coeff of eps^4, polynomial in n of order 1
      +1024, -208, 105105,
    // C4[3], coeff of eps^3, polynomial in n of order 2
      +3584, -3328, 1144, 315315,
    // C4[4], coeff of eps^5, polynomial in n of order 0
      -128, 135135,
    // C4[4], coeff of eps^4, polynomial in n of order 1
      -2560, 832, 405405,
    // C4[5], coeff of eps^5, polynomial in n of order 0
      +128, 99099,
  ];

  g.Geodesic.prototype.C4coeff = function() {
    var o = 0, k = 0,
        l, j, p;
    for (l = 0; l < g.nC4_; ++l) {        // l is index of C4[l]
      for (j = g.nC4_ - 1; j >= l; --j) { // coeff of eps^j
        p = g.nC4_ - j - 1;               // order of polynomial in n
        this._C4x[k++] = m.polyval(p, C4_coeff, o, this._n)
          / C4_coeff[o + p + 1];
        o += p + 2;
      }
    }
  };

  g.Geodesic.prototype.A3f = function(eps) {
    // Evaluate A3
    return m.polyval(nA3x_ - 1, this._A3x, 0, eps);
  };

  g.Geodesic.prototype.C3f = function(eps, c) {
    // Evaluate C3 coeffs
    // Elements c[1] thru c[nC3_ - 1] are set
    var mult = 1,
        o = 0,
        l, p;
    for (l = 1; l < g.nC3_; ++l) { // l is index of C3[l]
      p = g.nC3_ - l - 1;          // order of polynomial in eps
      mult *= eps;
      c[l] = mult * m.polyval(p, this._C3x, o, eps);
      o += p + 1;
    }
  };

  g.Geodesic.prototype.C4f = function(eps, c) {
    // Evaluate C4 coeffs
    // Elements c[0] thru c[g.nC4_ - 1] are set
    var mult = 1,
        o = 0,
        l, p;
    for (l = 0; l < g.nC4_; ++l) { // l is index of C4[l]
      p = g.nC4_ - l - 1;          // order of polynomial in eps
      c[l] = mult * m.polyval(p, this._C4x, o, eps);
      o += p + 1;
      mult *= eps;
    }
  };

  // return s12b, m12b, m0, M12, M21
  g.Geodesic.prototype.Lengths = function(eps, sig12,
                                          ssig1, csig1, dn1, ssig2, csig2, dn2,
                                          cbet1, cbet2,
                                          C1a, C2a) {
    // Return m12b = (reduced length)/_b; also calculate s12b =
    // distance/_b, and m0 = coefficient of secular term in
    // expression for reduced length.
    var vals = {},
        m0x = 0, J12 = 0, A1 = 0, A2 = 0,
        B1, B2, l, csig12, t;
    A1 = g.A1m1f(eps);
    g.C1f(eps, C1a);
    A1 = 1 + A1;
    B1 = g.SinCosSeries(true, ssig2, csig2, C1a) -
        g.SinCosSeries(true, ssig1, csig1, C1a);
    // Missing a factor of _b
    vals.s12b = A1 * (sig12 + B1);
    return vals;
  };

  // return sig12, salp1, calp1, salp2, calp2, dnm
  g.Geodesic.prototype.InverseStart = function(sbet1, cbet1, dn1,
                                               sbet2, cbet2, dn2, lam12,
                                               C1a, C2a) {
    // Return a starting point for Newton's method in salp1 and calp1
    // (function value is -1).  If Newton's method doesn't need to be
    // used, return also salp2 and calp2 and function value is sig12.
    // salp2, calp2 only updated if return val >= 0.
    var vals = {},
        // bet12 = bet2 - bet1 in [0, pi); bet12a = bet2 + bet1 in (-pi, 0]
        sbet12 = sbet2 * cbet1 - cbet2 * sbet1,
        cbet12 = cbet2 * cbet1 + sbet2 * sbet1,
        sbet12a, shortline, omg12, sbetm2, somg12, comg12, t, ssig12, csig12,
        x, y, lamscale, betscale, k2, eps, cbet12a, bet12a, m12b, m0, nvals,
        k, omg12a;
    vals.sig12 = -1;        // Return value
    // Volatile declaration needed to fix inverse cases
    // 88.202499451857 0 -88.202499451857 179.981022032992859592
    // 89.262080389218 0 -89.262080389218 179.992207982775375662
    // 89.333123580033 0 -89.333123580032997687 179.99295812360148422
    // which otherwise fail with g++ 4.4.4 x86 -O3
    sbet12a = sbet2 * cbet1;
    sbet12a += cbet2 * sbet1;

    shortline = cbet12 >= 0 && sbet12 < 0.5 && cbet2 * lam12 < 0.5;
    omg12 = lam12;
    if (shortline) {
      sbetm2 = m.sq(sbet1 + sbet2);
      // sin((bet1+bet2)/2)^2
      // =  (sbet1 + sbet2)^2 / ((sbet1 + sbet2)^2 + (cbet1 + cbet2)^2)
      sbetm2 /= sbetm2 + m.sq(cbet1 + cbet2);
      vals.dnm = Math.sqrt(1 + this._ep2 * sbetm2);
      omg12 /= this._f1 * vals.dnm;
    }
    somg12 = Math.sin(omg12); comg12 = Math.cos(omg12);

    vals.salp1 = cbet2 * somg12;
    vals.calp1 = comg12 >= 0 ?
      sbet12 + cbet2 * sbet1 * m.sq(somg12) / (1 + comg12) :
      sbet12a - cbet2 * sbet1 * m.sq(somg12) / (1 - comg12);

    ssig12 = m.hypot(vals.salp1, vals.calp1);
    csig12 = sbet1 * sbet2 + cbet1 * cbet2 * comg12;

    if (shortline && ssig12 < this._etol2) {
      // really short lines
      vals.salp2 = cbet1 * somg12;
      vals.calp2 = sbet12 - cbet1 * sbet2 *
        (comg12 >= 0 ? m.sq(somg12) / (1 + comg12) : 1 - comg12);
      // norm(vals.salp2, vals.calp2);
      t = m.hypot(vals.salp2, vals.calp2); vals.salp2 /= t; vals.calp2 /= t;
      // Set return value
      vals.sig12 = Math.atan2(ssig12, csig12);
    } else if (Math.abs(this._n) > 0.1 || // Skip astroid calc if too eccentric
               csig12 >= 0 ||
               ssig12 >= 6 * Math.abs(this._n) * Math.PI * m.sq(cbet1)) {
      // Nothing to do, zeroth order spherical approximation is OK
    } else {
      // Scale lam12 and bet2 to x, y coordinate system where antipodal
      // point is at origin and singular point is at y = 0, x = -1.
      if (this.f >= 0) {       // In fact f == 0 does not get here
        // x = dlong, y = dlat
        k2 = m.sq(sbet1) * this._ep2;
        eps = k2 / (2 * (1 + Math.sqrt(1 + k2)) + k2);
        lamscale = this.f * cbet1 * this.A3f(eps) * Math.PI;
        betscale = lamscale * cbet1;

        x = (lam12 - Math.PI) / lamscale;
        y = sbet12a / betscale;
      } else {                  // f < 0
        // x = dlat, y = dlong
        cbet12a = cbet2 * cbet1 - sbet2 * sbet1;
        bet12a = Math.atan2(sbet12a, cbet12a);
        // In the case of lon12 = 180, this repeats a calculation made
        // in Inverse.
        nvals = this.Lengths(this._n, Math.PI + bet12a,
                             sbet1, -cbet1, dn1, sbet2, cbet2, dn2,
                             cbet1, cbet2, g.REDUCEDLENGTH, C1a, C2a);
        m12b = nvals.m12b; m0 = nvals.m0;
        x = -1 + m12b / (cbet1 * cbet2 * m0 * Math.PI);
        betscale = x < -0.01 ? sbet12a / x :
          -this.f * m.sq(cbet1) * Math.PI;
        lamscale = betscale / cbet1;
        y = (lam12 - Math.PI) / lamscale;
      }

      if (y > -tol1_ && x > -1 - xthresh_) {
        // strip near cut
        if (this.f >= 0) {
          vals.salp1 = Math.min(1, -x);
          vals.calp1 = - Math.sqrt(1 - m.sq(vals.salp1));
        } else {
          vals.calp1 = Math.max(x > -tol1_ ? 0 : -1, x);
          vals.salp1 = Math.sqrt(1 - m.sq(vals.calp1));
        }
      } else {
        // Estimate alp1, by solving the astroid problem.
        //
        // Could estimate alpha1 = theta + pi/2, directly, i.e.,
        //   calp1 = y/k; salp1 = -x/(1+k);  for f >= 0
        //   calp1 = x/(1+k); salp1 = -y/k;  for f < 0 (need to check)
        //
        // However, it's better to estimate omg12 from astroid and use
        // spherical formula to compute alp1.  This reduces the mean number of
        // Newton iterations for astroid cases from 2.24 (min 0, max 6) to 2.12
        // (min 0 max 5).  The changes in the number of iterations are as
        // follows:
        //
        // change percent
        //    1       5
        //    0      78
        //   -1      16
        //   -2       0.6
        //   -3       0.04
        //   -4       0.002
        //
        // The histogram of iterations is (m = number of iterations estimating
        // alp1 directly, n = number of iterations estimating via omg12, total
        // number of trials = 148605):
        //
        //  iter    m      n
        //    0   148    186
        //    1 13046  13845
        //    2 93315 102225
        //    3 36189  32341
        //    4  5396      7
        //    5   455      1
        //    6    56      0
        //
        // Because omg12 is near pi, estimate work with omg12a = pi - omg12
        k = Astroid(x, y);
        omg12a = lamscale * ( this.f >= 0 ? -x * k/(1 + k) : -y * (1 + k)/k );
        somg12 = Math.sin(omg12a); comg12 = -Math.cos(omg12a);
        // Update spherical estimate of alp1 using omg12 instead of
        // lam12
        vals.salp1 = cbet2 * somg12;
        vals.calp1 = sbet12a -
          cbet2 * sbet1 * m.sq(somg12) / (1 - comg12);
      }
    }
    // Sanity check on starting guess.  Backwards check allows NaN through.
    if (!(vals.salp1 <= 0)) {
      // norm(vals.salp1, vals.calp1);
      t = m.hypot(vals.salp1, vals.calp1); vals.salp1 /= t; vals.calp1 /= t;
    } else {
      vals.salp1 = 1; vals.calp1 = 0;
    }
    return vals;
  };

  // return lam12, salp2, calp2, sig12, ssig1, csig1, ssig2, csig2, eps,
  // domg12, dlam12,
  g.Geodesic.prototype.Lambda12 = function(sbet1, cbet1, dn1, sbet2, cbet2, dn2,
                                           salp1, calp1, diffp,
                                           C1a, C2a, C3a) {
    var vals = {},
        t, salp0, calp0,
        somg1, comg1, somg2, comg2, omg12, B312, h0, k2, nvals;
    if (sbet1 === 0 && calp1 === 0)
      // Break degeneracy of equatorial line.  This case has already been
      // handled.
      calp1 = -g.tiny_;

    // sin(alp1) * cos(bet1) = sin(alp0)
    salp0 = salp1 * cbet1;
    calp0 = m.hypot(calp1, salp1 * sbet1); // calp0 > 0

    // tan(bet1) = tan(sig1) * cos(alp1)
    // tan(omg1) = sin(alp0) * tan(sig1) = tan(omg1)=tan(alp1)*sin(bet1)
    vals.ssig1 = sbet1; somg1 = salp0 * sbet1;
    vals.csig1 = comg1 = calp1 * cbet1;
    // norm(vals.ssig1, vals.csig1);
    t = m.hypot(vals.ssig1, vals.csig1); vals.ssig1 /= t; vals.csig1 /= t;
    // norm(somg1, comg1); -- don't need to normalize!

    // Enforce symmetries in the case abs(bet2) = -bet1.  Need to be careful
    // about this case, since this can yield singularities in the Newton
    // iteration.
    // sin(alp2) * cos(bet2) = sin(alp0)
    vals.salp2 = cbet2 !== cbet1 ? salp0 / cbet2 : salp1;
    // calp2 = sqrt(1 - sq(salp2))
    //       = sqrt(sq(calp0) - sq(sbet2)) / cbet2
    // and subst for calp0 and rearrange to give (choose positive sqrt
    // to give alp2 in [0, pi/2]).
    vals.calp2 = cbet2 !== cbet1 || Math.abs(sbet2) !== -sbet1 ?
      Math.sqrt(m.sq(calp1 * cbet1) + (cbet1 < -sbet1 ?
                                       (cbet2 - cbet1) * (cbet1 + cbet2) :
                                       (sbet1 - sbet2) * (sbet1 + sbet2))) /
      cbet2 : Math.abs(calp1);
    // tan(bet2) = tan(sig2) * cos(alp2)
    // tan(omg2) = sin(alp0) * tan(sig2).
    vals.ssig2 = sbet2; somg2 = salp0 * sbet2;
    vals.csig2 = comg2 = vals.calp2 * cbet2;
    // norm(vals.ssig2, vals.csig2);
    t = m.hypot(vals.ssig2, vals.csig2); vals.ssig2 /= t; vals.csig2 /= t;
    // norm(somg2, comg2); -- don't need to normalize!

    // sig12 = sig2 - sig1, limit to [0, pi]
    vals.sig12 = Math.atan2(Math.max(0, vals.csig1 * vals.ssig2 -
                                     vals.ssig1 * vals.csig2),
                            vals.csig1 * vals.csig2 + vals.ssig1 * vals.ssig2);

    // omg12 = omg2 - omg1, limit to [0, pi]
    omg12 = Math.atan2(Math.max(0, comg1 * somg2 - somg1 * comg2),
                       comg1 * comg2 + somg1 * somg2);
    k2 = m.sq(calp0) * this._ep2;
    vals.eps = k2 / (2 * (1 + Math.sqrt(1 + k2)) + k2);
    this.C3f(vals.eps, C3a);
    B312 = (g.SinCosSeries(true, vals.ssig2, vals.csig2, C3a) -
            g.SinCosSeries(true, vals.ssig1, vals.csig1, C3a));
    h0 = -this.f * this.A3f(vals.eps);
    vals.domg12 = salp0 * h0 * (vals.sig12 + B312);
    vals.lam12 = omg12 + vals.domg12;

    if (diffp) {
      if (vals.calp2 === 0)
        vals.dlam12 = - 2 * this._f1 * dn1 / sbet1;
      else {
        nvals = this.Lengths(vals.eps, vals.sig12,
                             vals.ssig1, vals.csig1, dn1,
                             vals.ssig2, vals.csig2, dn2,
                             cbet1, cbet2, g.REDUCEDLENGTH, C1a, C2a);
        vals.dlam12 = nvals.m12b;
        vals.dlam12 *= this._f1 / (vals.calp2 * cbet2);
      }
    }
    return vals;
  };

  /**
   * @summary Solve the inverse geodesic problem.
   * @param {number} lat1 the latitude of the first point in degrees.
   * @param {number} lon1 the longitude of the first point in degrees.
   * @param {number} lat2 the latitude of the second point in degrees.
   * @param {number} lon2 the longitude of the second point in degrees.
   * @returns {object} the requested results
   * @description The lat1, lon1, lat2, lon2, and a12 fields of the result are
   *   always set.
   */
  g.Geodesic.prototype.Inverse = function(lat1, lon1, lat2, lon2) {
    var vals = {},
        lon12, lonsign, t, swapp, latsign,
        sbet1, cbet1, sbet2, cbet2, s12x, m12x,
        dn1, dn2, lam12, slam12, clam12,
        sig12, calp1, salp1, calp2, salp2, C1a, C2a, C3a, meridian, nvals,
        ssig1, csig1, ssig2, csig2, eps, omg12, dnm,
        numit, salp1a, calp1a, salp1b, calp1b,
        tripn, tripb, v, dv, dalp1, sdalp1, cdalp1, nsalp1,
        salp0, calp0, alp12, k2, A4, C4a, B41, B42,
        somg12, domg12, dbet1, dbet2, salp12, calp12;
    // Compute longitude difference (AngDiff does this carefully).  Result is
    // in [-180, 180] but -180 is only for west-going geodesics.  180 is for
    // east-going and meridional geodesics.
    vals.lat1 = lat1 = m.LatFix(lat1); vals.lat2 = lat2 = m.LatFix(lat2);
    lon12 = m.AngDiff(lon1, lon2);
    vals.lon1 = m.AngNormalize(lon1); vals.lon2 = m.AngNormalize(lon2);
    // If very close to being on the same half-meridian, then make it so.
    lon12 = m.AngRound(lon12);
    // Make longitude difference positive.
    lonsign = lon12 >= 0 ? 1 : -1;
    lon12 *= lonsign;
    // If really close to the equator, treat as on equator.
    lat1 = m.AngRound(lat1);
    lat2 = m.AngRound(lat2);
    // Swap points so that point with higher (abs) latitude is point 1
    // If one latitude is a nan, then it becomes lat1.
    swapp = Math.abs(lat1) < Math.abs(lat2) ? -1 : 1;
    if (swapp < 0) {
      lonsign *= -1;
      t = lat1;
      lat1 = lat2;
      lat2 = t;
      // swap(lat1, lat2);
    }
    // Make lat1 <= 0
    latsign = lat1 < 0 ? 1 : -1;
    lat1 *= latsign;
    lat2 *= latsign;
    // Now we have
    //
    //     0 <= lon12 <= 180
    //     -90 <= lat1 <= 0
    //     lat1 <= lat2 <= -lat1
    //
    // longsign, swapp, latsign register the transformation to bring the
    // coordinates to this canonical form.  In all cases, 1 means no change was
    // made.  We make these transformations so that there are few cases to
    // check, e.g., on verifying quadrants in atan2.  In addition, this
    // enforces some symmetries in the results returned.

    t = m.sincosd(lat1); sbet1 = this._f1 * t.s; cbet1 = t.c;
    // norm(sbet1, cbet1);
    t = m.hypot(sbet1, cbet1); sbet1 /= t; cbet1 /= t;
    // Ensure cbet1 = +epsilon at poles
    cbet1 = Math.max(g.tiny_, cbet1);

    t = m.sincosd(lat2); sbet2 = this._f1 * t.s; cbet2 = t.c;
    // norm(sbet2, cbet2);
    t = m.hypot(sbet2, cbet2); sbet2 /= t; cbet2 /= t;
    // Ensure cbet2 = +epsilon at poles
    cbet2 = Math.max(g.tiny_, cbet2);

    // If cbet1 < -sbet1, then cbet2 - cbet1 is a sensitive measure of the
    // |bet1| - |bet2|.  Alternatively (cbet1 >= -sbet1), abs(sbet2) + sbet1 is
    // a better measure.  This logic is used in assigning calp2 in Lambda12.
    // Sometimes these quantities vanish and in that case we force bet2 = +/-
    // bet1 exactly.  An example where is is necessary is the inverse problem
    // 48.522876735459 0 -48.52287673545898293 179.599720456223079643
    // which failed with Visual Studio 10 (Release and Debug)

    if (cbet1 < -sbet1) {
      if (cbet2 === cbet1)
        sbet2 = sbet2 < 0 ? sbet1 : -sbet1;
    } else {
      if (Math.abs(sbet2) === -sbet1)
        cbet2 = cbet1;
    }

    dn1 = Math.sqrt(1 + this._ep2 * m.sq(sbet1));
    dn2 = Math.sqrt(1 + this._ep2 * m.sq(sbet2));

    lam12 = lon12 * m.degree;
    t = m.sincosd(lon12); slam12 = t.s; clam12 = t.c;

    // index zero elements of these arrays are unused
    C1a = new Array(g.nC1_ + 1);
    C2a = new Array(g.nC2_ + 1)
    C3a = new Array(g.nC3_);

    meridian = lat1 === -90 || slam12 === 0;
    if (meridian) {

      // Endpoints are on a single full meridian, so the geodesic might
      // lie on a meridian.

      calp1 = clam12; salp1 = slam12; // Head to the target longitude
      calp2 = 1; salp2 = 0;           // At the target we're heading north

      // tan(bet) = tan(sig) * cos(alp)
      ssig1 = sbet1; csig1 = calp1 * cbet1;
      ssig2 = sbet2; csig2 = calp2 * cbet2;

      // sig12 = sig2 - sig1
      sig12 = Math.atan2(Math.max(0, csig1 * ssig2 - ssig1 * csig2),
                         csig1 * csig2 + ssig1 * ssig2);
      nvals = this.Lengths(this._n, sig12,
                           ssig1, csig1, dn1, ssig2, csig2, dn2, cbet1, cbet2,
                           C1a, C2a);
      s12x = nvals.s12b;
      m12x = nvals.m12b;
      // Ignore m0
      // Add the check for sig12 since zero length geodesics might yield
      // m12 < 0.  Test case was
      //
      //    echo 20.001 0 20.001 0 | GeodSolve -i
      //
      // In fact, we will have sig12 > pi/2 for meridional geodesic
      // which is not a shortest path.
      if (sig12 < 1 || m12x >= 0) {
        // Need at least 2, to handle 90 0 90 180
        if (sig12 < 3 * g.tiny_)
          sig12 = m12x = s12x = 0;
        m12x *= this._b;
        s12x *= this._b;
        vals.a12 = sig12 / m.degree;
      } else
        // m12 < 0, i.e., prolate and too close to anti-podal
        meridian = false;
    }

    if (!meridian &&
        sbet1 === 0 &&           // and sbet2 == 0
        // Mimic the way Lambda12 works with calp1 = 0
        (this.f <= 0 || lam12 <= Math.PI - this.f * Math.PI)) {

      // Geodesic runs along equator
      calp1 = calp2 = 0; salp1 = salp2 = 1;
      s12x = this.a * lam12;
      sig12 = omg12 = lam12 / this._f1;
      m12x = this._b * Math.sin(sig12);
      vals.a12 = lon12 / this._f1;

    } else if (!meridian) {

      // Now point1 and point2 belong within a hemisphere bounded by a
      // meridian and geodesic is neither meridional or equatorial.

      // Figure a starting point for Newton's method
      nvals = this.InverseStart(sbet1, cbet1, dn1, sbet2, cbet2, dn2, lam12,
                                C1a, C2a);
      sig12 = nvals.sig12;
      salp1 = nvals.salp1;
      calp1 = nvals.calp1;

      if (sig12 >= 0) {
        salp2 = nvals.salp2;
        calp2 = nvals.calp2;
        // Short lines (InverseStart sets salp2, calp2, dnm)

        dnm = nvals.dnm;
        s12x = sig12 * this._b * dnm;
        m12x = m.sq(dnm) * this._b * Math.sin(sig12 / dnm);
        vals.a12 = sig12 / m.degree;
        omg12 = lam12 / (this._f1 * dnm);
      } else {

        // Newton's method.  This is a straightforward solution of f(alp1) =
        // lambda12(alp1) - lam12 = 0 with one wrinkle.  f(alp) has exactly one
        // root in the interval (0, pi) and its derivative is positive at the
        // root.  Thus f(alp) is positive for alp > alp1 and negative for alp <
        // alp1.  During the course of the iteration, a range (alp1a, alp1b) is
        // maintained which brackets the root and with each evaluation of
        // f(alp) the range is shrunk if possible.  Newton's method is
        // restarted whenever the derivative of f is negative (because the new
        // value of alp1 is then further from the solution) or if the new
        // estimate of alp1 lies outside (0,pi); in this case, the new starting
        // guess is taken to be (alp1a + alp1b) / 2.
        numit = 0;
        // Bracketing range
        salp1a = g.tiny_; calp1a = 1; salp1b = g.tiny_; calp1b = -1;
        for (tripn = false, tripb = false; numit < maxit2_; ++numit) {
          // the WGS84 test set: mean = 1.47, sd = 1.25, max = 16
          // WGS84 and random input: mean = 2.85, sd = 0.60
          nvals = this.Lambda12(sbet1, cbet1, dn1, sbet2, cbet2, dn2,
                                salp1, calp1, numit < maxit1_,
                                C1a, C2a, C3a);
          v = nvals.lam12 - lam12;
          salp2 = nvals.salp2;
          calp2 = nvals.calp2;
          sig12 = nvals.sig12;
          ssig1 = nvals.ssig1;
          csig1 = nvals.csig1;
          ssig2 = nvals.ssig2;
          csig2 = nvals.csig2;
          eps = nvals.eps;
          omg12 = nvals.domg12;
          dv = nvals.dlam12;

          // 2 * tol0 is approximately 1 ulp for a number in [0, pi].
          // Reversed test to allow escape with NaNs
          if (tripb || !(Math.abs(v) >= (tripn ? 8 : 2) * tol0_))
            break;
          // Update bracketing values
          if (v > 0 && (numit < maxit1_ || calp1/salp1 > calp1b/salp1b)) {
              salp1b = salp1; calp1b = calp1;
          } else if (v < 0 &&
                     (numit < maxit1_ || calp1/salp1 < calp1a/salp1a)) {
            salp1a = salp1; calp1a = calp1;
          }
          if (numit < maxit1_ && dv > 0) {
            dalp1 = -v/dv;
            sdalp1 = Math.sin(dalp1); cdalp1 = Math.cos(dalp1);
            nsalp1 = salp1 * cdalp1 + calp1 * sdalp1;
            if (nsalp1 > 0 && Math.abs(dalp1) < Math.PI) {
              calp1 = calp1 * cdalp1 - salp1 * sdalp1;
              salp1 = nsalp1;
              // norm(salp1, calp1);
              t = m.hypot(salp1, calp1); salp1 /= t; calp1 /= t;
              // In some regimes we don't get quadratic convergence because
              // slope -> 0.  So use convergence conditions based on epsilon
              // instead of sqrt(epsilon).
              tripn = Math.abs(v) <= 16 * tol0_;
              continue;
            }
          }
          // Either dv was not postive or updated value was outside legal
          // range.  Use the midpoint of the bracket as the next estimate.
          // This mechanism is not needed for the WGS84 ellipsoid, but it does
          // catch problems with more eccentric ellipsoids.  Its efficacy is
          // such for the WGS84 test set with the starting guess set to alp1 =
          // 90deg:
          // the WGS84 test set: mean = 5.21, sd = 3.93, max = 24
          // WGS84 and random input: mean = 4.74, sd = 0.99
          salp1 = (salp1a + salp1b)/2;
          calp1 = (calp1a + calp1b)/2;
          // norm(salp1, calp1);
          t = m.hypot(salp1, calp1); salp1 /= t; calp1 /= t;
          tripn = false;
          tripb = (Math.abs(salp1a - salp1) + (calp1a - calp1) < tolb_ ||
                   Math.abs(salp1 - salp1b) + (calp1 - calp1b) < tolb_);
        }
        nvals = this.Lengths(eps, sig12,
                             ssig1, csig1, dn1, ssig2, csig2, dn2, cbet1, cbet2,
                             C1a, C2a);
        s12x = nvals.s12b;
        m12x = nvals.m12b;
        // Ignore m0
        m12x *= this._b;
        s12x *= this._b;
        vals.a12 = sig12 / m.degree;
        omg12 = lam12 - omg12;
      }
    }

    vals.s12 = 0 + s12x;      // Convert -0 to 0

    // Returned value in [0, 180]
    return vals;
  };
})(GeographicLib.Geodesic, GeographicLib.Math, GeographicLib.Constants);

cb(GeographicLib);

})(function(geo) {
  if (typeof module === 'object' && module.exports) {
    /******** support loading with node's require ********/
    module.exports = geo;
  } else if (typeof define === 'function' && define.amd) {
    /******** support loading with AMD ********/
    define(/*'geographiclib',*/ [], function() { return geo; });
  } else {
    /******** otherwise just pollute our global namespace ********/
    window.GeographicLib = geo;
  }
});
