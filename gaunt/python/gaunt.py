# library to evaluate Gaunt coefficients
# Copyright (C) 2016 Michael Hartmann <michael.hartmann@physik.uni-augsburg.de>
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.


from __future__ import division
from math import lgamma,exp

def gaunt_a0(n,nu,m,mu):
    # eq. (20)
    return exp( lgamma(2*n+1)-lgamma(n+1)+lgamma(2*nu+1)-lgamma(1+nu)+lgamma(n+nu+1)-lgamma(2*n+2*nu+1) + lgamma(1+n+nu-m-mu)-lgamma(1+n-m)-lgamma(1+nu-mu) )


def gaunt(n,nu,m,mu):
    """
    Determine Gaunt coefficients a(m, n, mu, nu, p) for m, n, mu and nu fixed.
    These coefficients can be used to express the product of two associated
    Legendre polynomials:

    P_n^m(x)*P_{nu}^{mu}(x) = a0 sum_{q=0}^{qmax} aq_tilde P_{n+nu-2q}^(m+mu)(x)

    Returns: qmax, a0, aq_tilde
    qmax is the upper bound of summation, a0 is the prefactor and aq_tilde is a
    list of normalized Gaunt coefficients.

    See [1] for more information, especially chapter 3. There is a brief
    outline how to calculate Gaunt coefficients at the end of the chapter.

    Ref.: [1] Y.-L. Xu, J. Comp. Appl. Math. 85, 53 (1997)
    """

    # eq. (28)
    Ap = lambda p:p*(p-1)*(m-mu)-(m+mu)*(n-nu)*(n+nu+1)

    # eq. (3)
    alpha = lambda p: ((p**2-(n+nu+1)**2)*(p**2-(n-nu)**2))/(4*p**2-1)

    # eq. (24)
    qmax = min(n,nu,(n+nu-abs(m+mu))//2)

    a0 = gaunt_a0(n,nu,m,mu)

    n4 = n+nu-m-mu
    a_tilde = [0]*(qmax+1)

    a_tilde[0] = 1
    if qmax == 0:
        return qmax,a0,a_tilde

    # eq. (29)
    a_tilde[1] = (n+nu-1.5)*(1-(2*n+2*nu-1)/(n4*(n4-1))*((m-n)*(m-n+1)/(2*n-1)+(mu-nu)*(mu-nu+1)/(2*nu-1)))
    if qmax == 1:
        return qmax,a0,a_tilde

    # eq. (35)
    a_tilde[2] = (2*n+2*nu-1)*(2*n+2*nu-7)/4*( (2*n+2*nu-3)/(n4*(n4-1)) * ( (2*n+2*nu-5)/(2*(n4-2)*(n4-3)) \
                * ( (m-n)*(m-n+1)*(m-n+2)*(m-n+3)/((2*n-1)*(2*n-3)) \
                + 2*(m-n)*(m-n+1)*(mu-nu)*(mu-nu+1)/((2*n-1)*(2*nu-1)) \
                + (mu-nu)*(mu-nu+1)*(mu-nu+2)*(mu-nu+3)/((2*nu-1)*(2*nu-3)) ) - (m-n)*(m-n+1)/(2*n-1) \
                - (mu-nu)*(mu-nu+1)/(2*nu-1) ) +0.5)


    for q in range(3,qmax+1):
        p = n+nu-2*q
        p1 = p-m-mu
        p2 = p+m+mu

        if Ap(p+4) != 0:
            # eqs. (26), (27)
            c0 = (p+2)*(p+3)*(p1+1)*(p1+2)*Ap(p+4)*alpha(p+1)
            c1 = Ap(p+2)*Ap(p+3)*Ap(p+4) \
               + (p+1)*(p+3)*(p1+2)*(p2+2)*Ap(p+4)*alpha(p+2) \
               + (p+2)*(p+4)*(p1+3)*(p2+3)*Ap(p+2)*alpha(p+3)
            c2 = -(p+2)*(p+3)*(p2+3)*(p2+4)*Ap(p+2)*alpha(p+4)
            a_tilde[q] = (c1*a_tilde[q-1] + c2*a_tilde[q-2])/c0
        else:
            if Ap(p+6) == 0:
                # eq. (30)
                a_tilde[q] = (p+1)*(p2+2)*alpha(p+2)*a_tilde[q-1] / ((p+2)*(p1+1)*alpha(p+1))
            else:
                # eq. (32), (33)
                c0 = (p+2)*(p+3)*(p+5)*(p1+1)*(p1+2)*(p1+4)*Ap(p+6)*alpha(p+1)
                c1 = (p+5)*(p1+4)*Ap(p+6)*(Ap(p+2)*Ap(p+3)+(p+1)*(p+3)*(p1+2)*(p2+2)*alpha(p+2))
                c2 = (p+2)*(p2+3)*Ap(p+2)*(Ap(p+5)*Ap(p+6)+(p+4)*(p+6)*(p1+5)*(p2+5)*alpha(p+5))
                c3 = -(p+2)*(p+4)*(p+5)*(p2+3)*(p2+5)*(p2+6)*Ap(p+2)*alpha(p+6)
                a_tilde[q] = (c1*a_tilde[q-1] + c2*a_tilde[q-2] + c3*a_tilde[q-3])/c0

    return qmax,a0,a_tilde


if __name__ == "__main__":
    # XXX add tests
    l1,l2,M = 500,500,400

    qmax,a0,g = gaunt(l1,l2,M,M)

    print(g)
