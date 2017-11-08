cimport besselc

def I(int n, double x):
    """I(n,x): Compute modified Bessel function I_n(x)"""
    return besselc.besselI(n, x)

def I0e(double x):
    """I0e(x): Compute exponentially scaled modified Bessel function I_0(x)"""
    return besselc.besselI0e(x)

def I1e(double x):
    """I1e(x): Compute exponentially scaled modified Bessel function I_1(x)"""
    return besselc.besselI1e(x)

def lnInu(int nu, double x):
    """lnInu(nu,x): Compute logarithm of modified Bessel function I_(ν+1/2)(x)"""
    return besselc.bessel_lnInu(nu, x)

def lnKnu(int nu, double x):
    """lnKnu(nu,x): Compute logarithm of modified Bessel function K_(ν+1/2)(x)"""
    return besselc.bessel_lnKnu(nu, x)

def lnInuKnu(int nu, double x):
    """lnInuKnu(nu,x): Compute logarithm of modified Bessel functions I and K (see lnInu and lnKnu)"""
    cdef double Inu, Knu
    besselc.bessel_lnInuKnu(nu, x, &Inu, &Knu)
    return Inu, Knu

def continued_fraction(int nu, double x):
    """continued_fraction(nu,x): Compute continued fraction

    Compute the ratio of the modified Bessel functions of the first kind
    I_{ν+1/2}(x)/I_{ν+3/2}(x) using a continued fraction.
    """
    return besselc.bessel_continued_fraction(nu, x)
