cdef extern from "bessel.h":
    double besselI0e(double x)
    double besselI1e(double x)
    double besselI(int n, double x)
    double bessel_continued_fraction(int nu, double x)
    double bessel_lnInu(int nu, double x)
    double bessel_lnKnu(int nu, double x)
    void bessel_lnInuKnu(int nu, const double x, double *lnInu_p, double *lnKnu_p)
