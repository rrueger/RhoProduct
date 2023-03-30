#!/usr/bin/env python3

from RhoProduct import RhoProduct, repr_fractional_power, NotEtaProduct
import datetime
import sage.all
from sage.rings.infinity import Infinity
from sage.misc.misc_c import prod
from sage.rings.laurent_series_ring import LaurentSeriesRing
from sage.rings.integer_ring import ZZ
from sage.rings.big_oh import O
from sage.modular.etaproducts import EtaProduct
from sage.matrix.constructor import Matrix

# print("----------------------------------------------------------------------")
# print("Testing __repr__ of RhoProduct")
# print("----------------------------------------------------------------------")
# print(RhoProduct(8, {1: 24, 2: -24}))
# print(RhoProduct(8, {1: 24, 2: -24})*RhoProduct(6, {2: 24, 3: -24}))

# print("----------------------------------------------------------------------")
# print("Testing NotEtaProduct level/exponents check")
# print("----------------------------------------------------------------------")
# try:
#     RhoProduct(8, {3: 24, 7: -24})
# except NotEtaProduct as err:
#     print("RhoProduct(8, {3: 24, 7: -24})")
#     print(err)

# print("----------------------------------------------------------------------")
# print("Check that the Eta Function is correctly implemented")
# print("----------------------------------------------------------------------")
# # This is surprisingly difficult to do, for two reasons
# # 1. Sage's EtaProducts have no way of just printing the expansion of the
# #    EtaFunction
# # 2. Naively calculating (1-q^n) converges at a different rate that the fourier
# #    expansion. This means that one has to calculate both terms to a high
# #    precision to get comparable results.
# # Must be LaurentSeriesRing for coefficients() interface to work
# R = LaurentSeriesRing(ZZ, names=('q',)); (q,) = R._first_ngens(1)
# # Implicitly use q^24, so that we can properly divide later.
# p = prod([1 - q**(24*n) for n in range(1, 10)])
# # "Multiply" by q^(1/24)
# print("Naive product expansion:", repr_fractional_power(q*p + O(q**(10*24))))
# print("RhoProduct             :", RhoProduct(1, {1: 1}).q_expansion())

# print("----------------------------------------------------------------------")
# print("Directly comparing expansions of EtaProduct with RhoProduct")
# print("----------------------------------------------------------------------")
# print("EtaProduct:", EtaProduct(8, {1: 24, 2: -24}).q_expansion(10))
# print("RhoProduct:", RhoProduct(8, {1: 24, 2: -24}).q_expansion(10))

# print("----------------------------------------------------------------------")
# print("Calculating non-EtaProduct examples from Web of Modularity")
# print("----------------------------------------------------------------------")
# print("Ex. 1.66 Web of Modularity:", RhoProduct(5, {5: 5, 1: -1}).q_expansion())
# print("Ex. 1.66 Web of Modularity:", RhoProduct(8, {4: 2, 8: 2}).q_expansion())

# print("----------------------------------------------------------------------")
# print("Slashing by matrices")
# print("----------------------------------------------------------------------")
# print(RhoProduct(5, {5: 5, 1: -1}).slashby(Matrix([[4, 7], [1, 2]])))
# print()
# print(RhoProduct(1, {1: 24}).slashby(Matrix([[0, 1], [-1, 0]])))
# print(RhoProduct(1, {1: 24}).slashby(Matrix([[4, 7], [1, 2]])))
# print(RhoProduct(1, {1: 24}).q_expansion())
# print(RhoProduct(11, {1: 2, 11: 2}).q_expansion())
# print(RhoProduct(11, {1: 2, 11: 2}).isetaproduct())
# print(RhoProduct(11, {1: 2, 11: 2}).slashby(Matrix([[0, -1], [1, 0]])))
# print("AL Involution:", RhoProduct(8, {1: 24, 2: -24}).al_involution(8))


# print("----------------------------------------------------------------------")
# print("Test Accuracy and benchmarking speed of RhoProduct vs EtaProduct")
# print("----------------------------------------------------------------------")
# max_prec = 250
# level = 12
# exponents = {1: -336, 2: 576, 3: 696, 4: -216, 6: -576, 12: -144}
# print(":: Comparing EtaProduct and RhoProducts")
# print(EtaProduct(level, exponents))
# print(RhoProduct(level, exponents))
# print(f":: Calculating EtaProduct expansion up to precision {max_prec}...")
# etaproducts = []
# start = datetime.datetime.now()
# for prec in range(1, max_prec+1):
#     print(f'\r{prec}/{max_prec}', end='')
#     etaproducts.append(EtaProduct(level, exponents).q_expansion(prec))
# print()
# eta_diff = (datetime.datetime.now() - start)
# eta_time = eta_diff.seconds + float(eta_diff.microseconds / 10**9)
# eta_timeper = round(float(eta_time/max_prec), 5)
# print(f':: Time taken: {eta_time}s ({eta_timeper}s per calculation)')

# print(f":: Calculating RhoProduct expansion up to precision {max_prec}...")
# rhoproducts = []
# start = datetime.datetime.now()
# for prec in range(1, max_prec+1):
#     if datetime.datetime.now() - start > eta_diff:
#         factor = round((datetime.datetime.now() - start)/eta_diff * max_prec/prec, 2)
#         print(f'\r{prec}/{max_prec} (Slower than EtaProduct by factor {factor})', end='')
#     else:
#         print(f'\r{prec}/{max_prec}', end='')
#     rhoproducts.append(RhoProduct(level, exponents).q_expansion(prec))
# print()
# rho_diff = (datetime.datetime.now() - start)
# rho_time = rho_diff.seconds + float(rho_diff.microseconds / 10**9)
# rho_timeper = round(float(rho_time/max_prec), 5)
# print(f':: Time taken: {rho_time}s ({rho_timeper}s per calculation)')

# if not etaproducts == rhoproducts:
#     print(":: Not Equal")
#     print(":: First 5 product expansions")
#     for etaproduct, rhoproduct in zip(etaproducts[:5], rhoproducts[:5]):
#         print("EtaProduct:",  etaproduct)
#         print("RhoProduct:",  rhoproduct)
# else:
#     print(":: Both results are equal at all precision levels")

# print(f'EtaProduct was {round(float(rho_time/eta_time), 5)}x faster')

R = LaurentSeriesRing(ZZ, names=('q',))
(q,) = R._first_ngens(1)
print(repr_fractional_power(q, 24))
