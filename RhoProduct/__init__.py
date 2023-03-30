#!/usr/bin/env python3

import sage.all
from sage.arith.functions import lcm
from sage.arith.misc import gcd
from sage.arith.misc import kronecker_symbol
from sage.arith.misc import xgcd
from sage.functions.other import ceil
from sage.matrix.constructor import Matrix
from sage.misc.functional import sqrt
from sage.misc.misc_c import prod
from sage.rings.big_oh import O
from sage.rings.infinity import Infinity
from sage.rings.integer_ring import ZZ
from sage.rings.laurent_series_ring import LaurentSeriesRing
from sage.rings.universal_cyclotomic_field import UniversalCyclotomicField


def simplify(a, b):
    g = gcd(a, b)
    return (a/g, b/g)


def coefficients(series):
    coeffs = dict()
    # Extract coefficients and exponents from the expansion
    # Whilst there is a simple list() interface for PowerSeriesRing, there
    # is not (a priori) one for LaurentSeriesRing
    # list() on a LaurentSeriesRing element will return the list of
    # coefficients, starting at the smallest non-zero exponent
    # Compare
    # R.<q> = LaurentSeriesRing(ZZ); list(q**5 + q**3 + q)
    # R.<q> = PolynomialRing(ZZ); list(q**5 + q**3 + q)
    # We use .exponents()[0] to extract the smallest non-zero exponent
    # So we can shift the coefficient
    lowest_exponent = series.exponents()[0]
    for exponent, coefficient in enumerate(list(series)):
        if coefficient:
            coeffs[exponent+lowest_exponent] = coefficient
    # Sort dictionary to ensure that we get the powers in ascending order
    return sorted(coeffs.items(), key=lambda i: [0])


def repr_fractional_power(series, k=24):
    """Print series in fractional exponents from whole exponents

    This is to overcome a Sage limitation of not allowing fractional powers in
    LaurentSeriesRings.
    """

    # Set precision of the input series
    # If a series is not truncated by O(z^k), then series.prec() returns
    # infinity. This is not what we want in thise case
    if series.prec() == Infinity:
        # Coefficients returns a sorted list of pairs (exponent, coefficient)
        prec = max(coefficients(series), key=lambda i: i[0])[0]
    else:
        prec = series.prec()

    # If integral exponents, we can return this as a LaurentSeriesRing object
    if all([i % k == 0 for i in series.exponents()]):
        reduced_series = 0
        R = LaurentSeriesRing(ZZ, names=('q',))
        (q,) = R._first_ngens(1)
        for coef, exp in zip(series.coefficients(),
                             [exp//k for exp in series.exponents()]):
            reduced_series += coef*q**exp

        reduced_series += O(q**ceil(prec/k))
        return reduced_series

    # Else, the exponents are not integral, and we cannot return a
    # LaurentSeriesRing, so we return a string
    fmt_str = []
    for exponent, coefficient in coefficients(series):
        # Ensuring that the fraciton n**2/k is simplified
        new_exp_numerator, new_exp_denominator = simplify(exponent, k)
        new_exponent = new_exp_numerator / new_exp_denominator

        # This case distinction is ugly, but required to avoid
        # 4x^2 + 3x^1 + 4x^0
        if exponent == 0:
            fmt_str += [f'{coefficient}']
        else:
            if coefficient == 1:
                c_str = ''
            elif coefficient == -1:
                c_str = '-'
            else:
                c_str = f'{coefficient}*'

            if new_exponent == 1:
                q_str = 'q'
            elif new_exp_denominator == 1:
                q_str = f'q^{new_exp_numerator}'
            else:
                # t = round(float(new_exp_numerator/new_exp_denominator), 2)
                t = new_exp_numerator/new_exp_denominator
                q_str = f'q^({t})'

            fmt_str += [c_str + q_str]

    # Sage cannot do calculations inside fstrings
    prec = ceil(prec/k)
    if series.prec() != Infinity:
        fmt_str += [f'O(q^{prec})']
    fmt_str = " + ".join(fmt_str)
    # Clean up negative coefficients and multiplying by 1
    fmt_str = fmt_str.replace('+ -', '- ')
    return fmt_str


class NonSlashableEtaProduct(Exception):
    pass


class NotEtaProduct(Exception):
    pass


class RhoProduct:
    """RhoProduct: Ryan's Eta Product

    Calling this RhoProduct allows us to compare to Sage's implementation of
    EtaProduct"""

    def __init__(self, level, exponents):
        self.level = level
        # Alias
        self.N = level
        self.exponents = exponents

        ndiv = [delta for delta in self.exponents.keys() if level % delta]
        if ndiv:
            if len(ndiv) == 1:
                err = f"{ndiv[0]} does not divide level ({level})"
            else:
                ndiv_str = ", ".join(map(str, ndiv))
                err = f"{ndiv_str} do not divide level ({level})"

            rhoproduct = f"RhoProudct({level}, {repr(exponents)})"
            raise NotEtaProduct(f"Wrong level/exponents in {rhoproduct}: {err}")

        self.k = sum(exponents.values())//2
        # Requires int() call to prevent sage from casting the elements as
        # floats and runing into overflow errors (prod is a sage function)
        self.char = prod([int(d**exp) for (d, exp) in self.exponents.items()])
        self.char *= (-1)**self.k

    def __eta_function(self, prec, alpha=1, beta=0):
        # Alias for debugging performance
        level = self.level
        if beta:
            U = UniversalCyclotomicField()
            E = U.gen
            R = LaurentSeriesRing(U, names=('q',))
            (q,) = R._first_ngens(1)
            self.expansion = 0
            for n in range(ceil(sqrt(prec*level*24))):
                # Here, use q, since we want a power of q^(1/24)
                print(n**2*beta)
                self.expansion += kronecker_symbol(12, n) \
                                  * E(24)**(n**2*beta)    \
                                  * q**(alpha*level*n**2)
        else:
            R = LaurentSeriesRing(ZZ, names=('q',))
            (q,) = R._first_ngens(1)
            self.expansion = 0
            for n in range(ceil(sqrt(prec*level*24))):
                self.expansion += kronecker_symbol(12, n) * q**(alpha*level*n**2)

        # Here we do need q**24 since we want prec to represent in integral
        # powers of q
        self.expansion += O(q**(24*level*prec))
        return self.expansion

    def __coefficients(self):
        return coefficients(self.expansion)

    def __repr__(self):
        # Repr in the style of Sage's EtaProduct with addl true repr
        product = " ".join([f'(eta_{d})^{e}' for d, e in self.exponents.items()])
        interpretation = f'Rho Product of level {self.level} : {product}'
        true_repr = f'RhoProduct({self.N}, {repr(self.exponents)})'
        return interpretation + f' ({true_repr})'

        # # Informal __repr__, in the style of Sage
        # product_terms = [f'n(z)^{e}' if d == 1 else f'n({d}z)^{e}'
        #                  for d, e in self.exponents.items()]
        # product = "*".join(product_terms)
        # interpretation = f' is the RhoProduct {product} of level {self.N}'

        # # True __repr__, i.e. one can type this get the object
        # true_repr = f'RhoProduct({self.N}, {repr(self.exponents)})'
        # return true_repr + interpretation

    def __pow__(self, power):
        exponents = {d: e*power for d, e in self.exponents.items()}
        return RhoProduct(self.level, exponents)

    def __mul__(self, other):
        level = lcm(self.level, other.level)
        # Create new exponents dict. We do not want to change self or other
        exponents = dict()
        for delta, exponent in self.exponents.items():
            exponents[delta] = exponent
        for delta, exponent in other.exponents.items():
            if delta in exponents:
                exponents[delta] += exponent
            else:
                exponents[delta] = exponent
        return RhoProduct(level, exponents)

    def isetaproduct(self):
        # For integral weight for \Gamma_0(N)
        if sum(self.exponents.values()) % 2:
            # print("The sum of exponents is not divisible by 2")
            return False
        elif sum([d*exp for d, exp in self.exponents.items()]) % 24:
            # print("The sum of d*a_d is not divisible by 24")
            return False
        elif sum([self.N/d*exp for d, exp in self.exponents.items()]) % 24:
            # print("The sum of N/d*a_d is not divisible by 24")
            return False

        return True

    def isslashable(self):
        # For integral weight for \Gamma_0(N)
        if sum(self.exponents.values()) % 2:
            print("The sum of exponents is not divisible by 2")
            return False
        return True

    def q_expansion(self, prec=10):
        self.expansion = prod([self.__eta_function(prec, alpha=delta)**exponent
                               for delta, exponent in self.exponents.items()])
        return repr_fractional_power(self.expansion, 24*self.level)

    def expansionatcusp(self, cusp):
        pass

    def order(self, cusp):
        d = cusp.denominator()
        s = sum([gcd(d, delta)*exp/delta for delta, exp in self.exponents.items()])
        order = self.N/24 * 1/d * 1/gcd(d, self.N/d) * s
        return order

    def slashby(self, matrix, prec=10):
        a = matrix[0][0]
        b = matrix[0][1]
        c = matrix[1][0]
        d = matrix[1][1]

        if not a*d - b*c == 1:
            raise NonSlashableEtaProduct("Matrix is not in SL_2(Z)")

        if not self.isslashable():
            raise NonSlashableEtaProduct("Cannot slash by matrix")

        U = UniversalCyclotomicField()
        R = LaurentSeriesRing(U, names=('q',))
        (q,) = R._first_ngens(1)
        p = 1

        for delta, exponent in self.exponents.items():
            g = gcd(c, delta)
            _, y, x = xgcd(-c/delta, delta*a/g)
            nu = x*delta*b + y*d
            p *= self.__eta_function(prec,
                                     alpha=1/delta,
                                     beta=nu*g/delta)**int(exponent)
                                     # alpha=1,
                                     # alpha=g**2/delta,

        # Need to cast this properly
        return repr_fractional_power(p, 24*self.level)

    # Since the fricke involution is only available as an expansion for now, we
    # must pass a precision paramater to it
    def al_involution(self, Q, matrix=None, prec=10):
        if gcd(self.level, Q) != Q:
            print("Q does not divide level!")
            exit(1)
        elif gcd(self.level/Q, Q) != 1:
            print("Q is not an _exact_ divisor of level!")
            exit(1)

        if matrix is None:
            a = Q
            b = 1
            _, d, c = xgcd(Q, -self.level/Q)
            c *= self.level
            d *= Q
        else:
            a = matrix[0][0]
            b = matrix[0][1]
            c = matrix[1][0]
            d = matrix[1][1]

        # First "slash" with (Q, 0; 0 ,1)
        # Do this by replacing z with Qz
        exponents = {Q * delta: exponent for delta, exponent in
                     self.exponents.items()}
        # Now try to slash the new "eta product" by the reduced matrix
        return RhoProduct(self.level*Q, exponents
                          ).slashby(Matrix([[a/Q, b], [c/Q, d]]))
