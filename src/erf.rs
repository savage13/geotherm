
/// Error Function Trait
///
/// Implementation based on [Netlib Specfun 2.5](http://www.netlib.org/specfun/)
///
/// Translated from Fortran March 2019
///
/// ## References
///  - Cody, W. (1990), Performance evaluation of programs for the error
///      and complementary error functions, ACM Transactions on Mathematical
///      Software (TOMS), 16(1), 29–37.
///
pub trait ErrorFunction {
    /// Computes the Error Function
    fn erf(self) -> Self ;
    /// Computes Complementary Error Function
    fn erfc(self) -> Self ;
    /// Computes exp(x^2) * erfc(x).
    fn erfcx(self) -> Self ;
}

impl ErrorFunction for f64 {
    fn erf(self)   -> Self { crate::fad::erf(self)   }
    fn erfc(self)  -> Self { crate::fad::erfc(self)  }
    fn erfcx(self) -> Self { crate::fad::erfcx(self) }
}

#[derive(PartialEq)]
enum ErfKind {
    Erf,
    Erfc,
    Erfcx,
}

fn calerf(arg: f64, kind: ErfKind) -> f64 { 
    //------------------------------------------------------------------
    //
    // This packet evaluates  erf(x),  erfc(x),  and  exp(x*x)*erfc(x)
    //   for a real argument  x.  It contains three FUNCTION type
    //   subprograms: ERF, ERFC, and ERFCX (or DERF, DERFC, and DERFCX),
    //   and one SUBROUTINE type subprogram, CALERF.  The calling
    //   statements for the primary entries are:
    //
    //                   Y=ERF(X)     (or   Y=DERF(X)),
    //
    //                   Y=ERFC(X)    (or   Y=DERFC(X)),
    //   and
    //                   Y=ERFCX(X)   (or   Y=DERFCX(X)).
    //
    //   The routine  CALERF  is intended for internal packet use only,
    //   all computations within the packet being concentrated in this
    //   routine.  The function subprograms invoke  CALERF  with the
    //   statement
    //
    //          CALL CALERF(ARG,RESULT,JINT)
    //
    //   where the parameter usage is as follows
    //
    //      Function                     Parameters for CALERF
    //       call              ARG                  Result          JINT
    //
    //     ERF(ARG)      ANY REAL ARGUMENT         ERF(ARG)          0
    //     ERFC(ARG)     ABS(ARG) .LT. XBIG        ERFC(ARG)         1
    //     ERFCX(ARG)    XNEG .LT. ARG .LT. XMAX   ERFCX(ARG)        2
    //
    //   The main computation evaluates near-minimax approximations
    //   from "Rational Chebyshev approximations for the error function"
    //   by W. J. Cody, Math. Comp., 1969, PP. 631-638.  This
    //   transportable program uses rational functions that theoretically
    //   approximate  erf(x)  and  erfc(x)  to at least 18 significant
    //   decimal digits.  The accuracy achieved depends on the arithmetic
    //   system, the compiler, the intrinsic functions, and proper
    //   selection of the machine-dependent constants.
    //
    //*******************************************************************
    //*******************************************************************
    //
    // Explanation of machine-dependent constants
    //
    //   XMIN   = the smallest positive floating-point number.
    //   XINF   = the largest positive finite floating-point number.
    //   XNEG   = the largest negative argument acceptable to ERFCX;
    //            the negative of the solution to the equation
    //            2*exp(x*x) = XINF.
    //   XSMALL = argument below which erf(x) may be represented by
    //            2*x/sqrt(pi)  and above which  x*x  will not underflow.
    //            A conservative value is the largest machine number X
    //            such that   1.0 + X = 1.0   to machine precision.
    //   XBIG   = largest argument acceptable to ERFC;  solution to
    //            the equation:  W(x) * (1-0.5/x**2) = XMIN,  where
    //            W(x) = exp(-x*x)/[x*sqrt(pi)].
    //   XHUGE  = argument above which  1.0 - 1/(2*x*x) = 1.0  to
    //            machine precision.  A conservative value is
    //            1/[2*sqrt(XSMALL)]
    //   XMAX   = largest acceptable argument to ERFCX; the minimum
    //            of XINF and 1/[sqrt(pi)*XMIN].
    //
    //   Approximate values for some important machines are:
    //
    //                          XMIN       XINF        XNEG     XSMALL
    //
    //  CDC 7600      (S.P.)  3.13E-294   1.26E+322   -27.220  7.11E-15
    //  CRAY-1        (S.P.)  4.58E-2467  5.45E+2465  -75.345  7.11E-15
    //  IEEE (IBM/XT,
    //    SUN, etc.)  (S.P.)  1.18E-38    3.40E+38     -9.382  5.96E-8
    //  IEEE (IBM/XT,
    //    SUN, etc.)  (D.P.)  2.23D-308   1.79D+308   -26.628  1.11D-16
    //  IBM 195       (D.P.)  5.40D-79    7.23E+75    -13.190  1.39D-17
    //  UNIVAC 1108   (D.P.)  2.78D-309   8.98D+307   -26.615  1.73D-18
    //  VAX D-Format  (D.P.)  2.94D-39    1.70D+38     -9.345  1.39D-17
    //  VAX G-Format  (D.P.)  5.56D-309   8.98D+307   -26.615  1.11D-16
    //
    //
    //                          XBIG       XHUGE       XMAX
    //
    //  CDC 7600      (S.P.)  25.922      8.39E+6     1.80X+293
    //  CRAY-1        (S.P.)  75.326      8.39E+6     5.45E+2465
    //  IEEE (IBM/XT,
    //    SUN, etc.)  (S.P.)   9.194      2.90E+3     4.79E+37
    //  IEEE (IBM/XT,
    //    SUN, etc.)  (D.P.)  26.543      6.71D+7     2.53D+307
    //  IBM 195       (D.P.)  13.306      1.90D+8     7.23E+75
    //  UNIVAC 1108   (D.P.)  26.582      5.37D+8     8.98D+307
    //  VAX D-Format  (D.P.)   9.269      1.90D+8     1.70D+38
    //  VAX G-Format  (D.P.)  26.569      6.71D+7     8.98D+307
    //
    //*******************************************************************
    //*******************************************************************
    //
    // Error returns
    //
    //  The program returns  ERFC = 0      for  ARG .GE. XBIG;
    //
    //                       ERFCX = XINF  for  ARG .LT. XNEG;
    //      and
    //                       ERFCX = 0     for  ARG .GE. XMAX.
    //
    // Intrinsic functions required are:
    //
    //     ABS, AINT, EXP
    //
    //  Author: W. J. Cody
    //          Mathematics and Computer Science Division
    //          Argonne National Laboratory
    //          Argonne, IL 60439
    //
    //  Latest modification: March 19, 1990
    //
    //      INTEGER I,JINT
    //      DOUBLE PRECISION
    //     1     A,ARG,B,C,D,DEL,FOUR,HALF,P,ONE,Q,RESULT,SIXTEN,SQRPI,
    //     2     TWO,THRESH,X,XBIG,XDEN,XHUGE,XINF,XMAX,XNEG,XNUM,XSMALL,
    //     3     Y,YSQ,ZERO
    //      DIMENSION A(5),B(4),C(9),D(8),P(6),Q(5)
    //------------------------------------------------------------------
    //  Mathematical constants
    //------------------------------------------------------------------
    //S    DATA FOUR,ONE,HALF,TWO,ZERO/4.0E0,1.0E0,0.5E0,2.0E0,0.0E0/,
    //S   1     SQRPI/5.6418958354775628695E-1/,THRESH/0.46875E0/,
    //S   2     SIXTEN/16.0E0/
    //      DATA FOUR,ONE,HALF,TWO,ZERO/4.0D0,1.0D0,0.5D0,2.0D0,0.0D0/,
    //     1     SQRPI/5.6418958354775628695D-1/,THRESH/0.46875D0/,
    //     2     SIXTEN/16.0D0/
    let one    = 1.0;
    let two    = 2.0;
    let zero   = 0.0;
    let half   = 0.5;
    let thresh = 0.46875;
    let four   = 4.0;
    let sixten = 16.0;
    let sqrpi  = 5.6418958354775628695e-1; // 1/sqrt(pi)
    //------------------------------------------------------------------
    //  Machine-dependent constants
    //------------------------------------------------------------------
    //S    DATA XINF,XNEG,XSMALL/3.40E+38,-9.382E0,5.96E-8/,
    //S   1     XBIG,XHUGE,XMAX/9.194E0,2.90E3,4.79E37/
    //      DATA XINF,XNEG,XSMALL/1.79D308,-26.628D0,1.11D-16/,
    //     1     XBIG,XHUGE,XMAX/26.543D0,6.71D7,2.53D307/
    let xinf = 1.79e308;
    let xneg = -26.628e0;
    let xsmall = 1.11e-16;
    let xbig = 26.543e0;
    let xhuge = 6.71e7;
    let xmax = 2.53e307;
    
    //------------------------------------------------------------------
    //  Coefficients for approximation to  erf  in first interval
    //------------------------------------------------------------------
    //S    DATA A/3.16112374387056560E00,1.13864154151050156E02,
    //S   1       3.77485237685302021E02,3.20937758913846947E03,
    //S   2       1.85777706184603153E-1/
    //S    DATA B/2.36012909523441209E01,2.44024637934444173E02,
    //S   1       1.28261652607737228E03,2.84423683343917062E03/
    let a = [3.16112374387056560e00,1.13864154151050156e02,
             3.77485237685302021e02,3.20937758913846947e03,
             1.85777706184603153e-1];
    let b = [2.36012909523441209e01,2.44024637934444173e02,
             1.28261652607737228e03,2.84423683343917062e03];
    //------------------------------------------------------------------
    //  Coefficients for approximation to  erfc  in second interval
    //------------------------------------------------------------------
    //S    DATA C/5.64188496988670089E-1,8.88314979438837594E0,
    //S   1       6.61191906371416295E01,2.98635138197400131E02,
    //S   2       8.81952221241769090E02,1.71204761263407058E03,
    //S   3       2.05107837782607147E03,1.23033935479799725E03,
    //S   4       2.15311535474403846E-8/
    //S    DATA D/1.57449261107098347E01,1.17693950891312499E02,
    //S   1       5.37181101862009858E02,1.62138957456669019E03,
    //S   2       3.29079923573345963E03,4.36261909014324716E03,
    //S   3       3.43936767414372164E03,1.23033935480374942E03/
    let c = [5.64188496988670089e-1,8.88314979438837594e0,
             6.61191906371416295e01,2.98635138197400131e02,
             8.81952221241769090e02,1.71204761263407058e03,
             2.05107837782607147e03,1.23033935479799725e03,
             2.15311535474403846e-8];
    let d = [1.57449261107098347e01,1.17693950891312499e02,
             5.37181101862009858e02,1.62138957456669019e03,
             3.29079923573345963e03,4.36261909014324716e03,
             3.43936767414372164e03,1.23033935480374942e03];
    //------------------------------------------------------------------
    //  Coefficients for approximation to  erfc  in third interval
    //------------------------------------------------------------------
    //S    DATA P/3.05326634961232344E-1,3.60344899949804439E-1,
    //S   1       1.25781726111229246E-1,1.60837851487422766E-2,
    //S   2       6.58749161529837803E-4,1.63153871373020978E-2/
    //S    DATA Q/2.56852019228982242E00,1.87295284992346047E00,
    //S   1       5.27905102951428412E-1,6.05183413124413191E-2,
    //S   2       2.33520497626869185E-3/
    let p = [3.05326634961232344e-1,3.60344899949804439e-1,
             1.25781726111229246e-1,1.60837851487422766e-2,
             6.58749161529837803e-4,1.63153871373020978e-2];
    let q = [2.56852019228982242e00,1.87295284992346047e00,
             5.27905102951428412e-1,6.05183413124413191e-2,
             2.33520497626869185e-3];
    //------------------------------------------------------------------
    let x = arg;
    let y = x.abs();
    let mut result = 0.0;
    if y <= thresh {
        //
        //  Evaluate  erf  for  |X| <= 0.46875
        //
        let mut ysq = 0.0;
        if y > xsmall {
            ysq = y * y
        }
        let mut xnum = a[5-1]*ysq;
        let mut xden = ysq;
        for i in 0 .. 3 {
            xnum = (xnum + a[i]) * ysq;
            xden = (xden + b[i]) * ysq;
        }
        result = x * (xnum + a[4-1]) / (xden + b[4-1]);

        if kind != ErfKind::Erf {
            result = one - result;
        }
        if kind == ErfKind::Erfcx  {
            result = ysq.exp() * result;
        }
        return result;
    } else if y <= four {
        //
        //  Evaluate  erfc  for 0.46875 <= |X| <= 4.0
        //
        let mut xnum = c[9-1]*y;
        let mut xden = y;
        for i in 0 .. 7 {
            xnum = (xnum + c[i]) * y;
            xden = (xden + d[i]) * y;
        }
        result = (xnum + c[8-1]) / (xden + d[8-1]);
        if kind != ErfKind::Erfcx {
            let ysq = (y*sixten).trunc()/sixten;
            let del = (y-ysq)*(y+ysq);
            result = (-ysq*ysq).exp() * (-del).exp() * result;
        }
    } else {
        //
        //  Evaluate  erfc  for |X| > 4.0
        //
        let mut goto_300 = false;
        if y >= xbig {
            if kind != ErfKind::Erfcx || y >= xmax {
                goto_300 = true;
            }
            if y >= xhuge {
                result = sqrpi / y;
                //go to 300
                goto_300 = true;
            }
        }
        if ! goto_300 {
            let ysq = one / (y * y);
            let mut xnum = p[6-1]*ysq;
            let mut xden = ysq;
            for i in 0 .. 4 {
                xnum = (xnum + p[i]) * ysq;
                xden = (xden + q[i]) * ysq;
            }
            result = ysq *(xnum + p[5-1]) / (xden + q[5-1]);
            result = (sqrpi -  result) / y;
            if kind != ErfKind::Erfcx {
                let ysq = (y*sixten).trunc()/sixten;
                let del = (y-ysq)*(y+ysq);
                result = (-ysq*ysq).exp() * (-del).exp() * result;
            }
        }
    }
    //
    //  Fix up for negative argument, erf, etc.
    //
    if kind == ErfKind::Erf {
        // line 300
        result = (half - result) + half;
        if x < zero {
            result = -result;
        }
    } else if kind == ErfKind::Erfc {
        if x < zero {
            result = two - result
        }
    } else { // ErfKind::Erfcx
        if x < zero {
            if x < xneg {
                result = xinf;
            } else {
                let ysq = (x*sixten).trunc()/sixten;
                let del = (x-ysq)*(x+ysq);
                let y = (ysq*ysq).exp() * del.exp();
                result = (y+y) - result;
            }
        }
    }
    result
}

// This subprogram computes approximate values for erf(x).
//   (see comments heading CALERF).
//
//   Author/date: W. J. Cody, January 8, 1985
pub fn erf(x: f64) -> f64 {
    calerf(x,ErfKind::Erf)
}

// This subprogram computes approximate values for erfc(x).
//   (see comments heading CALERF).
//
//   Author/date: W. J. Cody, January 8, 1985
pub fn erfc(x: f64) -> f64 {
    calerf(x, ErfKind::Erfc)
}

// This subprogram computes approximate values for exp(x*x) * erfc(x).
//   (see comments heading CALERF).
//
//   Author/date: W. J. Cody, March 30, 1987
pub fn erfcx(x: f64) -> f64 {
    calerf(x, ErfKind::Erfcx)
}

#[cfg(test)]
mod tests {
    //use crate::erf::ErrorFunction;
    // use special_fun::FloatSpecial;
    // #[test]
    // fn erf_cephes_test() {
    //     let eps = 1.5 * std::f64::EPSILON;
    //     // Check vs special_fun::erf() (Cephes C library Wrapper)
    //     let scale = 1000000;

    //     for i in 0 .. 27 * scale {
    //         let v = i as f64 / scale as f64;
    //         if v > 26.641 {
    //             break;
    //         }
    //         let a = v.erf();
    //         let b = crate::erf::erf( v );
    //         assert!((a-b).abs() <= eps, "{} != {} [{:e}] {:e}", a,b,(a-b).abs(), eps);
    //     }

    // }
    // #[test]
    // fn erfc_cephes_test() {
    //     let eps = 2.0 * std::f64::EPSILON;
    //     // Check vs special_fun::erfc() (Cephes C library Wrapper)
    //     let scale = 1000000;

    //     for i in 0 .. 27 * scale {
    //         let v = i as f64 / scale as f64;
    //         if v > 26.641 {
    //             break;
    //         }
    //         let a = v.erfc();
    //         let b = crate::erf::erfc( v );
    //         assert!((a-b).abs() <= eps, "{} != {} [{:e}] {:e}", a,b,(a-b).abs(), eps);
    //     }
    // }
    #[inline]
    fn conv<T>(v: T) -> f64 where f64: std::convert::From<T> {
        f64::from(v)
    }
    #[inline]
    fn max<T: PartialOrd>(a: T, b: T) -> T { if a >= b { a } else { b } }
    #[inline]
    fn min<T: PartialOrd>(a: T, b: T) -> T { if a <= b { a } else { b } }

    fn ren(k: i64) -> f64 {
        //---------------------------------------------------------------------
        //  random number generator - based on algorithm 266 by pike and
        //   hill (modified by hansson), communications of the acm,
        //   vol. 8, no. 10, october 1965.
        //
        //  this subprogram is intended for use on computers with
        //   fixed point wordlength of at least 29 bits.  it is
        //   best if the floating-point significand has at most
        //   29 bits.
        //
        //  latest modification: may 30, 1989
        //
        //  author: w. j. cody
        //          mathematics and computer science division
        //          argonne national laboratory
        //          argonne, il 60439

        let one = 1.0;
        let c1 = 2796203.0e0;
        let c2 = 1.0e-6;
        let c3 = 1.0e-12;
        let iy = 100001;
        let _j = k;
        let iy = iy * 125;
        let iy = iy - (iy/2796203) * 2796203;
        let ren = conv(iy) / c1 * (one + c2 + c3);
        ren
    }

    #[test]
    fn erftst() {
        let ibeta = std::f64::RADIX;
        let it = std::f64::MANTISSA_DIGITS;
        let eps = std::f64::EPSILON;
        let xmin = std::f64::MIN_POSITIVE;
        let c1 = 5.6418958354775628695e-1;
        let one = 1.0_f64;
        let half = 0.5;
        let ten = 10.0;
        let zero = 0.0;
        let sixten = 16.0;
        let x99 = -999.0;
        let thresh = 0.46875e0;
        let two = 2.0;
        let beta = conv(ibeta);
        let albeta = (beta).ln();
        let ait = conv(it);
        let c = (ait*albeta).abs() + ten;
        let mut a = zero;
        let mut b = thresh;
        let n = 2000;
        let xn = conv(n);
        let jt = 0;
        let mut n0 = ( (ibeta/2)*(it+5)/6+4 ) as usize;
        let mut r1 = vec![0.0; 501];
        let xmax = std::f64::MAX;
        let func1 = crate::fad::erf;
        let func2 = crate::fad::erfc;
        let func3 = crate::fad::erfcx;

        assert_eq!(xmin, 2.2250738585072014E-308);  // Fragile
        //-----------------------------------------------------------------
        //  Determine largest argument for ERFC test by Newton iteration
        //-----------------------------------------------------------------
        let c2 = xmin.ln() + (one/c1).ln();
        let mut xbig = (-c2).sqrt();
        loop {
            let x = xbig; // 50
            let f0 = x*x;
            let ff = f0 + half/f0 + x.ln() + c2;
            let f0 = x+x + one/x - one/(x*f0);
            xbig = x - ff/f0;
            if (x-xbig).abs()/x <= ten*eps {//go to 50
                break;
            }
        }
        let mut w;
        let mut v;
        let mut u;
        let mut z;
        let mut r;
        let mut f0;
        let mut ff;
        assert_eq!(xbig, 26.543258430632299);
        for j in 1 ..= 5 { // 300 /// Tests 1-5
            let mut k1 = 0;
            let mut k3 = 0;
            let mut xc = zero;
            let mut r6 = zero;
            let mut r7 = zero;
            let del = (b - a) / xn;
            let mut xl = a;
            for _i in 0 .. n { // 200 // Number of Tests to conduct
                let mut x = del * ren(jt) + xl;
                if j == 1 {
                    //-----------------------------------------------------------------
                    //  test erf against double series expansion
                    //-----------------------------------------------------------------
                    f0 = conv(n0 as i32);
                    ff = f0+f0+one; 
                    z = x*x;
                    w = z+z;
                    u = zero;
                    v = zero;
                    for _k in 1 ..= n0  { // 60
                        u = -z/f0*(one+u);
                        v = w/ff*(one+v);
                        f0 = f0 - one;
                        ff = ff - two;
                    }
                    v = c1*(x+x)*(((u*v+(u+v))+half)+half);
                    u = func1(x);
                } else {
                    //-----------------------------------------------------------------
                    //  test erfc or scaled erfc against expansion in repeated
                    //   integrals of the coerror function.
                    //-----------------------------------------------------------------
                    z = x + half;
                    x = z - half;
                    r = zero;
                    if x <= one {
                        n0 = 499;
                    } else {
                        //n0 = min(499,int(c/(abs(log(z)))));
                        n0 = min(499.0, (c / z.ln().abs()).floor() ) as usize;
                    }
                    let mut n1 = n0 ;
                    let mut xn1 = conv((n1+1) as i32);
                    for _k in 1 ..= n0 {
                        r = half/(z+xn1*r);
                        r1[n1] = r;
                        n1 = n1 - 1;
                        xn1 = xn1 - one;
                    }//  100             continue
                    let mut k = n0;
                    ff = c1/(z+r1[1]);
                    if (j/2)*2 == j {
                        f0 = func2(z) * (x+half*half).exp();
                        u = func2(x);
                    } else {
                        f0 = func3(z);
                        u = func3(x);
                    }
                    let sc = f0/ff;
                    //-----------------------------------------------------------------
                    //  scale things to avoid premature underflow
                    //-----------------------------------------------------------------
                    let epscon = f0 ;
                    ff = sixten*ff/eps;
                    for n1 in 1 ..= n0 {
                        ff = r1[n1]*ff;
                        r1[n1] = ff * sc;
                        if r1[n1] < epscon {
                            k = n1;
                            break; //go to 111
                        }
                    } //110             continue
                    v = r1[k]; // 111
                    for n1 in 1 ..= k-1 {
                        v = v + r1[k-n1];
                    }
                    //120             continue
                    //-----------------------------------------------------------------
                    //  remove scaling here
                    //-----------------------------------------------------------------
                    v = v*eps/sixten + f0;
                }
                //--------------------------------------------------------------------
                //  accumulate results
                //--------------------------------------------------------------------
                w = (u - v) / u; // Relative Error
                let fac = 450.0;
                assert!(w.abs() < fac*std::f64::EPSILON, "{:e} != {:e} diff {:e} {:e} {}",
                        u,v, w.abs(), fac*std::f64::EPSILON, w.abs()/std::f64::EPSILON);
                if w > zero {
                    k1 = k1 + 1;
                } else if w < zero {
                    k3 = k3 + 1;
                }
                w = w.abs();
                if w > r6 {
                    r6 = w;
                    xc = x;
                }
                r7 = r7 + w * w;
                xl = xl + del;
            }//  200    continue
            //------------------------------------------------------------------
            //  gather and print statistics for test
            //------------------------------------------------------------------
            let k2 = n - k3 - k1;
            r7 = (r7/xn).sqrt();
            println!("TEST #{}", j);
            if j == 1 {
                println!(" Test of erf(x) vs double series expansion\n\n");
                println!("   {} Random arguments were tested from the interval ({}, {})", n,a,b);
                println!("   ERF(X) was larger {} times", k1);
            } else if (j/2)*2 == j {
                println!(" Test of erfc(x) vs exp(x+1/4) SUM i^n erfc(x+1/2)");
                println!("   {} Random arguments were tested from the interval ({}, {})", n,a,b);
                println!("ERFC(X) was larger {} times", k1);
            } else {
                println!(" Test of exp(x*x) erfc(x) vs SUM i^n erfc(x+1/2)");
                println!("{} Random arguments were tested from the interval ({}, {})", n,a,b);
                println!("ERFCX(X) was larger {} times", k1);
            }
            println!("              agreed {} times, and",k2);
            println!("         was smaller {} times.\n", k3);
            println!(" There are {} base {} significant digits in a floating-point number", it,ibeta);
            if r6 != zero {
                w = r6.abs().ln() / albeta;
            } else {
                w = x99;
            }
            println!(" The maximum relative error of {:.4e} = {} ** {}\n    occurred for X = {:.8e}",
                     r6, ibeta, w, xc);
            w = max(ait+w,zero);
            println!(" The estimated loss of base {} significant digits is {}", ibeta, w);
            if r7 != zero {
                w = r7.abs().ln() / albeta;
            } else {
                w = x99;
            }
            println!(" The root mean square relative error was {:.4e} = {} ** {}", r7, ibeta, w);
            w = max(ait+w,zero);
            println!(" The estimated loss of base {} significant digits is {}", ibeta, w);
            println!("");
            // -----------------------------------------------------------------
            //  initialize for next test
            //------------------------------------------------------------------
            if j == 1 {
                a = b;
                b = two;
            } else if j == 3 {
                a = b;
                b = (xbig*sixten).trunc()/sixten-half;
            } else if j == 4 {
                b = ten + ten;
            }
        } //300 continue

        println!("Special Tests");
        println!(" Estimated loss of base {} significant digits in", ibeta);
        println!("'Erf(x)+Erf(-x)   Erf(x)+Erfc(x)-1   Erfcx(x)-exp(x*x)*erfc(x)");

        let ans = [0.2763263901682369,
                   -0.2763263901682369,
                   0.7236736098317631,
                    1.276326390168237];
        assert!((func1(0.25) - ans[0]).abs() < 1e-15);
        assert!((func1(-0.25) - ans[1]).abs() < 1e-15);
        assert!((func2(0.25) - ans[2]).abs() < 1e-15);
        assert!((func2(-0.25) - ans[3]).abs() < 1e-15, "{} != {}", func2(-0.25),ans[3]);

        let arg = 0.5;
        let ans = [0.5204998778130465,
                   -0.5204998778130465,
                   0.4795001221869535,
                   1.5204998778130465];
        assert!((func1(arg) - ans[0]).abs() < 1e-15);
        assert!((func1(-arg) - ans[1]).abs() < 1e-15);
        assert!((func2(arg) - ans[2]).abs() < 1e-15);
        assert!((func2(-arg) - ans[3]).abs() < 1e-15, "{} != {}", func2(-arg),ans[3]);

        let mut x = zero;
        let del = -half;
        for _i in 1 ..= 10 {
            let u = func1(x);
            let mut a = u + func1(-x);
            assert!(a.abs() < 2.0 * std::f64::EPSILON);
            if a*u != zero {
                a = ait + (a/u).abs().ln() / albeta;
            }
            let v = func2(x);
            let mut b = u + v - one;
            assert!(b.abs() < 2.0 * std::f64::EPSILON);
            if b != zero {
                b = ait + b.abs().ln() / albeta;
            }
            let w = func3(x);
            let mut c = (x*sixten).trunc() / sixten;
            let r = (x-c)*(x+c);
            c = ((c*c).exp() * r.exp() * v-w)/w;
            assert!(c.abs() < 2.0 * std::f64::EPSILON);
            if c != zero {
                c = max(zero,ait + c.abs().ln() / albeta);
            }
            println!("{:7.3} {:16.2} {:16.2} {:16.2}", x,a,b,c);
            //write (iout,1031) x,a,b,c
            x = x + del;
        }

        println!("Test of special arguments");
        let z = xmax;
        let zz = func1(z);
        assert!((zz-1.0).abs() < std::f64::EPSILON, "{} != {}", zz, 1.0);
        println!("   ERF ({:e}) = {:e}", z,zz);
        let z = zero;
        let zz = func1(z);
        assert!((zz-0.0).abs() < std::f64::EPSILON);
        println!("   ERF ({:e}) = {:e}", z,zz);
        let zz = func2(z);
        assert!((zz-1.0).abs() < std::f64::EPSILON);
        println!("  ERFC ({:e}) = {:e}", z,zz);
        let z = -xmax;
        let zz = func2(z);
        assert!((zz-2.0).abs() < std::f64::EPSILON);
        println!("  ERFC ({:e}) = {:e}", z,zz);


        println!("Test of Error Returns");
        let w = xbig;
        let z = w * (one - half * half);
        println!("ERFC will be called with the argument {:e}", z);
        println!("   This should **not** underflow");
        let zz = func2(z);
        assert!((zz-0.217879e-173).abs() < std::f64::EPSILON);
        println!("  ERFC ({:e}) = {:e}", z,zz);
        println!("");

        let z = w * (one + ten * eps);
        println!("ERFC will be called with the argument {:e}", z);
        println!("   This **may** underflow");
        let zz = func2(z);
        assert!((zz-0.222508e-307).abs() < std::f64::EPSILON);
        println!("  ERFC ({:e}) = {:e}", z,zz);
        println!("");

        let mut w = xmax;
        if c1 < xmax*xmin {
            w = c1/xmin;
        }
        let z = w * (one - one/sixten);
        println!("ERFCX will be called with the argument {:e}", z);
        println!("   This should **not** underflow");
        let zz = func3(z);
        println!(" ERFCX ({:e}) = {:e}", z,zz);
        assert!((zz-0.237341e-307).abs() < std::f64::EPSILON);
        println!("");

        let w = -(xmax/two).ln().sqrt();
        let z = w * (one-one/ten);
        println!("ERFCX will be called with the argument {:e}", z);
        println!("   This should **not** overflow");
        let zz = func3(z);
        println!(" ERFCX ({:e}) = {:e}", z,zz);
        let ans = 0.5540070644707187e+250;
        let ans = 0.5540070644707037e+250;
        assert!((zz-ans).abs() < 1e-15,
                "{:e} != {:e} diff {:e}",
                zz,ans,(zz-ans).abs() );
        println!("");

        let z = w * (one + ten*eps);
        println!("ERFCX will be called with the argument {:e}", z);
        println!("   This **may** overflow");
        let zz = func3(z);
        println!(" ERFCX ({:e}) = {:e}", z,zz);
        let ans = 0.179000e+309;
        assert!(zz.is_infinite());
        println!("");

    }
}
