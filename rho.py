# Dickman's rho function (to compute probability of success of ecm).
#
# Copyright 2004, 2005, 2006, 2008, 2009, 2010, 2011, 2012, 2013
# Alexander Kruppa, Paul Zimmermann.
#
# This file translated directly from a part of the ECM Library.
#
# The ECM Library is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation; either version 3 of the License, or (at your
# option) any later version.
#
# The ECM Library is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
# License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with the ECM Library; see the file COPYING.LIB.  If not, see
# http://www.gnu.org/licenses/ or write to the Free Software Foundation, Inc.,
# 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.


import math
import os
import sys
from functools import cache

import gmpy2
import logging

EXTRA_SMOOTHNESS_SQUARE = 0.416384512396064
EXTRA_SMOOTHNESS_32BITS_D = 0.330484606500389
ECM_EXTRA_SMOOTHNESS = 3.134
DEFAULT_B2_EXPONENT = 1.43
ECM_COST = 11.0 / 6.0

ECM_PARAM_DEFAULT = -1
ECM_PARAM_SUYAMA = 0
ECM_PARAM_BATCH_SQUARE = 1
ECM_PARAM_BATCH_2 = 2
ECM_PARAM_BATCH_32BITS_D = 3
ECM_PARAM_WEIERSTRASS = 5
ECM_PARAM_HESSIAN = 6
ECM_PARAM_TWISTED_HESSIAN = 7
ECM_PARAM_TORSION = 8

NTT_SIZE_THRESHOLD = 30

ECM_ERROR = -1
ECM_DEFAULT_K = 0
MPZSPV_NORMALISE_STRIDE = 512
SP_NUMB_BITS = 62

M_PI_SQR_6 = 1.644934066848226436
M_PI_SQR = 9.869604401089358619

M_EULER = 0.577215664901532861
M_EULER_1 = 0.422784335098467139

primemap = [
  254, 223, 239, 126, 182, 219, 61, 249, 213, 79, 30, 243, 234, 166, 237, 158,
  230, 12, 211, 211, 59, 221, 89, 165, 106, 103, 146, 189, 120, 30, 166, 86,
  86, 227, 173, 45, 222, 42, 76, 85, 217, 163, 240, 159, 3, 84, 161, 248, 46,
  253, 68, 233, 102, 246, 19, 58, 184, 76, 43, 58, 69, 17, 191, 84, 140, 193,
  122, 179, 200, 188, 140, 79, 33, 88, 113, 113, 155, 193, 23, 239, 84, 150,
  26, 8, 229, 131, 140, 70, 114, 251, 174, 101, 146, 143, 88, 135, 210, 146,
  216, 129, 101, 38, 227, 160, 17, 56, 199, 38, 60, 129, 235, 153, 141, 81,
  136, 62, 36, 243, 51, 77, 90, 139, 28, 167, 42, 180, 88, 76, 78, 38, 246,
  25, 130, 220, 131, 195, 44, 241, 56, 2, 181, 205, 205, 2, 178, 74, 148, 12,
  87, 76, 122, 48, 67, 11, 241, 203, 68, 108, 36, 248, 25, 1, 149, 168, 92,
  115, 234, 141, 36, 150, 43, 80, 166, 34, 30, 196, 209, 72, 6, 212, 58, 47,
  116, 156, 7, 106, 5, 136, 191, 104, 21, 46, 96, 85, 227, 183, 81, 152, 8,
  20, 134, 90, 170, 69, 77, 73, 112, 39, 210, 147, 213, 202, 171, 2, 131, 97,
  5, 36, 206, 135, 34, 194, 169, 173, 24, 140, 77, 120, 209, 137, 22, 176, 87,
  199, 98, 162, 192, 52, 36, 82, 174, 90, 64, 50, 141, 33, 8, 67, 52, 182,
  210, 182, 217, 25, 225, 96, 103, 26, 57, 96, 208, 68, 122, 148, 154, 9, 136,
  131, 168, 116, 85, 16, 39, 161, 93, 104, 30, 35, 200, 50, 224, 25, 3, 68,
  115, 72, 177, 56, 195, 230, 42, 87, 97, 152, 181, 28, 10, 104, 197, 129,
  143, 172, 2, 41, 26, 71, 227, 148, 17, 78, 100, 46, 20, 203, 61, 220, 20,
  197, 6, 16, 233, 41, 177, 130, 233, 48, 71, 227, 52, 25, 195, 37, 10, 48,
  48, 180, 108, 193, 229, 70, 68, 216, 142, 76, 93, 34, 36, 112, 120, 146,
  137, 129, 130, 86, 38, 27, 134, 233, 8, 165, 0, 211, 195, 41, 176, 194, 74,
  16, 178, 89, 56, 161, 29, 66, 96, 199, 34, 39, 140, 200, 68, 26, 198, 139,
  130, 129, 26, 70, 16, 166, 49, 9, 240, 84, 47, 24, 210, 216, 169, 21, 6, 46,
  12, 246, 192, 14, 80, 145, 205, 38, 193, 24, 56, 101, 25, 195, 86, 147, 139,
  42, 45, 214, 132, 74, 97, 10, 165, 44, 9, 224, 118, 196, 106, 60, 216, 8,
  232, 20, 102, 27, 176, 164, 2, 99, 54, 16, 49, 7, 213, 146, 72, 66, 18, 195,
  138, 160, 159, 45, 116, 164, 130, 133, 120, 92, 13, 24, 176, 97, 20, 29, 2,
  232, 24, 18, 193, 1, 73, 28, 131, 48, 103, 51, 161, 136, 216, 15, 12, 244,
  152, 136, 88, 215, 102, 66, 71, 177, 22, 168, 150, 8, 24, 65, 89, 21, 181,
  68, 42, 82, 225, 179, 170, 161, 89, 69, 98, 85, 24, 17, 165, 12, 163, 60,
  103, 0, 190, 84, 214, 10, 32, 54, 107, 130, 12, 21, 8, 126, 86, 145, 1, 120,
  208, 97, 10, 132, 168, 44, 1, 87, 14, 86, 160, 80, 11, 152, 140, 71, 108,
  32, 99, 16, 196, 9, 228, 12, 87, 136, 11, 117, 11, 194, 82, 130, 194, 57,
  36, 2, 44, 86, 37, 122, 49, 41, 214, 163, 32, 225, 177, 24, 176, 12, 138,
  50, 193, 17, 50, 9, 197, 173, 48, 55, 8, 188, 145, 130, 207, 32, 37, 107,
  156, 48, 143, 68, 38, 70, 106, 7, 73, 142, 9, 88, 16, 2, 37, 197, 196, 66,
  90, 128, 160, 128, 60, 144, 40, 100, 20, 225, 3, 132, 81, 12, 46, 163, 138,
  164, 8, 192, 71, 126, 211, 43, 3, 205, 84, 42, 0, 4, 179, 146, 108, 66, 41,
  76, 131, 193, 146, 204, 28
]


invh = 0
tablemax = 0
h = 0.
rhotable = None

# void ret
# int parm_invh
# int parm_tablemax
def rhoinit(parm_invh, parm_tablemax):
    global invh, tablemax, rhotable, h

    if parm_invh == invh and parm_tablemax == tablemax:
        return

    if rhotable is not None:
        rhotable = None
        invh = 0
        h = 0.
        tablemax = 0

    # The integration below expects 3 * invh > 4
    if parm_tablemax == 0 or parm_invh < 2:
        return

    invh = parm_invh
    h = 1. / float(invh)
    tablemax = parm_tablemax

    rhotable = [None] * (parm_invh * parm_tablemax)  # (double *) malloc (parm_invh * parm_tablemax * sizeof (double))

    i = 0
    while i < (3 if 3 < parm_tablemax else parm_tablemax) * invh:
        rhotable[i] = rhoexact (i * h)
        i += 1

    i = 3 * invh
    while i < parm_tablemax * invh:
        # rho(i*h) = 1 - \int_{1}^{i*h} rho(x-1)/x dx
        #            = rho((i-4)*h) - \int_{(i-4)*h}^{i*h} rho(x-1)/x dx

        rhotable[i] = rhotable[i - 4] - 2. / 45. * (
            7. * rhotable[i - invh - 4] / float(i - 4)
          + 32. * rhotable[i - invh - 3] / float(i - 3)
          + 12. * rhotable[i - invh - 2] / float(i - 2)
          + 32. * rhotable[i - invh - 1] / float(i - 1)
          + 7. * rhotable[i - invh]  / float(i) )
        if rhotable[i] < 0.:
            rhotable[i] = 0.
        i += 1


param = 1
n_default = pow(2, 200)


# public api
@cache
def ecmprob(digits, B1, B2=None, param=1, n=n_default):
    if B2 is None:
        # Default B2
        B2 = pow(B1 * ECM_COST, DEFAULT_B2_EXPONENT)

    if param == ECM_PARAM_SUYAMA or param == ECM_PARAM_BATCH_2:
        smoothness_correction = 1.0
    elif param == ECM_PARAM_BATCH_SQUARE:
        smoothness_correction = EXTRA_SMOOTHNESS_SQUARE
    elif param == ECM_PARAM_BATCH_32BITS_D:
        smoothness_correction = EXTRA_SMOOTHNESS_32BITS_D
    else:  # This case should never happen
        smoothness_correction = 0.0

    B2len = B2 - B1
    if B2len < 1e7:
        S = 1    # x^1
    elif B2len < 1e8:
        S = 2    # x^2
    elif B2len < 1e9:
        S = -3   # Dickson(3)
    elif B2len < 1e10:
        S = -6   # Dickson(6)
    elif B2len < 3e11:
        S = -12  # Dickson(12)
    else:
        S = -30  # Dickson(30)

    # estimating dF from https://www.mersenneforum.org/showthread.php?p=641459#post641459
    # until I convert bestd.c to python
    # dF = 1 << math.floor(math.log2(math.sqrt(B2)))
    # k is optimally 4 I believe
    # k = 4

    N = pow(10, digits - 0.5) / smoothness_correction
    B2min = B1
    use_ntt = 1 if gmpy2.bit_length(n)//gmpy2.mp_limbsize() <= NTT_SIZE_THRESHOLD else 0
    po2 = 1 if use_ntt else 0
    maxmem = 0
    treefile = None
    modulus = None
    dF, k, B2_ = bestD(B2min, B2, po2, use_ntt, maxmem, treefile, modulus)
    # logging.debug(f"dF = {dF}, k = {k}, B2 = {B2}, B2_ = {B2_}")
    nr = dF * dF * k

    if rhotable is None:
        rhoinit(256, 10)

    return _ecmprob(B1, int(B2_), N, nr, S)



# ret double
# double B1
# double B2
# double N
# double nr
# int S
@cache
def _ecmprob (B1, B2, N, nr, S):
    return prob(B1, B2, N, nr, S, ECM_EXTRA_SMOOTHNESS)


# Assume N is as likely smooth as a number around N/exp(delta)
# ret double
# double B1
# double B2
# double N
# double nr
# int S
# double delta
@cache
def prob(B1, B2, N, nr, S, delta):
    sumthresh = 20000
    effN = N / math.exp(delta)

    assert(rhotable is not None)

    # What to do if rhotable is not initialised and asserting is not enabled?
    # For now, bail out with 0. result. Not really pretty, either
    if rhotable is None:
        return 0.

    if B1 < 2. or N <= 1.:
        return 0.

    if (effN <= B1):
        return 1.

    # logging.debug(f"B1 = {B1:f}, B2 = {B2:f}, N = {N:.0f}, nr = {nr:f}, S = {S}")

    alpha = math.log(effN) / math.log(B1)
    stage1 = dickmanlocal(alpha, effN)
    stage2 = 0.
    if B2 > B1:
        if B1 < sumthresh:
            stage2 += dickmanmu_sum(B1, min(B2, sumthresh), effN)
            beta = math.log(B2) / math.log(min(B2, sumthresh))
        else:
            beta = math.log(B2) / math.log(B1)

        if beta > 1.:
            stage2 += dickmanmu(alpha, beta, effN)
    brsu = 0.
    if S < -1:
        brsu = brsudickson(B1, B2, effN, nr, -S * 2)
    if S > 1:
        brsu = brsupower(B1, B2, effN, nr, S * 2)


    # logging.getLogger().debug(f"stage 1 : {stage1:f}, stage 2 : {stage2:f}, Brent-Suyama : {brsu:f}")

    return (stage1 + stage2 + brsu) if (stage1 + stage2 + brsu) > 0. else 0.


# return the value of the "local" Dickman rho function, for numbers near x
# (as opposed to numbers <= x for the original Dickman rho function).
# Reference: PhD thesis of Alexander Kruppa,
# http://docnum.univ-lorraine.fr/public/SCD_T_2010_0054_KRUPPA.pdf,
# equation (5.6) page 100

# ret double
# double alpha
# double x
@cache
def dickmanlocal(alpha, x):
    if alpha <= 1.:
        return rhoexact(alpha)
    if alpha < tablemax:
        return dickmanrho(alpha) - M_EULER * dickmanrho(alpha - 1.) / math.log(x)
    return 0.


# assumes alpha < tablemax
# ret double
# double alpha
@cache
def dickmanrho(alpha):
    assert(alpha < tablemax)

    if alpha <= 3.:
        return rhoexact (alpha)
    a = math.floor(alpha * invh)
    rho1 = rhotable[a]
    rho2 = rhotable[a + 1] if (a + 1) < tablemax * invh else 0
    return rho1 + (rho2 - rho1) * (alpha * invh - float(a))


# return the sum in Equation (5.10) page 102 of Alexander Kruppa's PhD thesis
# ret double
# const unsigned long B1
# const unsigned long B2
# const double x
@cache
def dickmanmu_sum(B1, B2, x):
    s = 0.
    inv_logB1 = 1. / math.log(B1)
    logx = math.log(x)
    p = B1 + 1

    while p <= B2:
        if isprime(p):
            s += dickmanlocal((logx - math.log(p)) * inv_logB1, x // p) // p
        p += 1

    return s



# return the probability that a number < x has its 2nd largest prime factor
# less than x^(1/alpha) and its largest prime factor less than x^(beta/alpha)
# ret double
# double alpha
# double beta
# double x
@cache
def dickmanmu(alpha, beta, x):
    #double a, b, sum
    #int ai, bi, i
    ai = math.ceil((alpha - beta) * invh)
    if ai > tablemax * invh:
      ai = tablemax * invh
    a = float(ai) * h
    bi = math.floor((alpha - 1.) * invh)
    if bi > tablemax * invh:
        bi = tablemax * invh
    b = float(bi) * h
    sum_ = 0.
    i = ai + 1
    while i < bi:
      sum_ += dickmanlocal_i(i, x) / (alpha - i * h)
      i += 1
    sum_ += 0.5 * dickmanlocal_i(ai, x) / (alpha - a)
    sum_ += 0.5 * dickmanlocal_i(bi, x) / (alpha - b)
    sum_ *= h
    sum_ += (a - alpha + beta) * 0.5 * (dickmanlocal_i(ai, x) / (alpha - a) + dickmanlocal (alpha - beta, x) / beta)
    sum_ += (alpha - 1. - b) * 0.5 * (dickmanlocal(alpha - 1., x) + dickmanlocal_i(bi, x) / (alpha - b))

    return sum_


# ret double
# int ai
# double x
@cache
def dickmanlocal_i(ai, x):
    if ai <= 0:
        return 0.
    if ai <= invh:
        return 1.
    if ai <= 2 * invh and ai < tablemax * invh:
        return rhotable[ai] - M_EULER / gmpy2.log(x)
    if ai < tablemax * invh:
        logx = math.log(x)
        return rhotable[ai] - (M_EULER * rhotable[ai - invh] + M_EULER_1 * rhotable[ai - 2 * invh] / logx) / logx

    return 0.


# ret double
# double x
@cache
def rhoexact(x):
    assert(x <= 3.)
    if x <= 0.:
        return 0.
    elif x <= 1.:
        return 1.
    elif x <= 2.:
        return 1. - math.log(x)
    else:  # 2 < x <= 3 thus -2 <= 1-x < -1
        return 1. - math.log(x) * (1. - math.log(x - 1.)) + dilog (1. - x) + 0.5 * M_PI_SQR_6


# ret double
# const double z
def dilog_series (z):
    r = 0.0
    k = 1
    k2 = 1
    zk = z
    # Doubles have 53 bits in significand, with |z| <= 0.5 the k+1-st term
    # is <= 1/(2^k k^2) of the result, so 44 terms should do
    while k <= 44:
        r += zk / float(k2)
        k2 += 2 * k + 1
        k += 1
        zk *= z
    return r


# ret double
# double x
def dilog (x):
    assert(x <= -1.0) # dilog(1-x) is called from rhoexact for 2 < x <= 3

    if (x <= -2.0):
        return -dilog_series (1./x) - M_PI_SQR_6 - 0.5 * math.log(-1./x) * math.log(-1./x)
    else: # x <= -1.0
        # L2(z) = -L2(1 - z) + 1/6 * Pi^2 - ln(1 - z)*ln(z)
        # L2(z) = -L2(1/z) - 1/6 * Pi^2 - 0.5*ln^2(-1/z)
        # ->
        # L2(z) = -(-L2(1/(1-z)) - 1/6 * Pi^2 - 0.5*ln^2(-1/(1-z))) + 1/6 * Pi^2 - ln(1 - z)*ln(z)
        #       = L2(1/(1-z)) - 1/6 * Pi^2 + 0.5*ln(1 - z)^2 - ln(1 - z)*ln(-z)
        # z in [-1, -2) -> 1/(1-z) in [1/2, 1/3)
        log1x = math.log(1. - x)
        return dilog_series(1. / (1. - x)) - M_PI_SQR_6 + log1x * (0.5 * log1x - math.log (-x))


# return the sum in Equation (5.10) page 102 of Alexander Kruppa's PhD thesis
# ret double
# const unsigned long B1
# const unsigned long B2
# const double x
@cache
def dickmanmu_sum(B1, B2, x):
    B1 = int(B1)
    B2 = int(B2)
    s = 0.
    inv_logB1 = 1. / math.log(B1)
    logx = math.log(x)
    p = B1 + 1

    while p <= B2:
        if isprime(p):
            s += dickmanlocal((logx - math.log(p)) * inv_logB1, x // p) // p
        p += 1
    return s


# ret int
# unsigned long n
@cache
def isprime(n):
    if n % 2 == 0:
        return n == 2
    if n % 3 == 0:
        return n == 3
    if n % 5 == 0:
        return n == 5

    if n // 30 >= len(primemap):
        # is this enough? do we need to exit?
        os.abort()

    r = n % 30 # 8 possible values: 1,7,11,13,17,19,23,29
    r = (r * 16 + r) // 64 # maps the 8 values onto 0, ..., 7

    return (primemap[n // 30] & (1 << r)) != 0


# ret double
# double B1
# double B2
# double N
# double nr
@cache
def brentsuyama(B1, B2, N, nr):
    #double a, alpha, beta, sum_;
    #int ai, i;
    alpha = math.log(N) / math.log(B1)
    beta = math.log(B2) / math.log(B1)
    ai = math.floor((alpha - beta) * invh)
    if ai > tablemax * invh:
        ai = tablemax * invh
    a = float(ai) * h
    sum_ = 0.
    i = 1
    while i < ai:
        sum_ += dickmanlocal_i(i, N) / (alpha - i * h) * (1 - math.exp(-nr * pow(B1, (-alpha + i * h))))
        i += 1
    sum_ += 0.5 * (1 - math.exp(-nr / pow(B1, alpha)))
    sum_ += 0.5 * dickmanlocal_i(ai, N) / (alpha - a) * (1 - math.exp(-nr * pow(B1, (-alpha + a))))
    sum_ *= h
    sum_ += 0.5 * (alpha - beta - a) * (dickmanlocal_i (ai, N) / (alpha - a) + dickmanlocal(alpha - beta, N) / beta)
    return sum_


# ret double
# double B1
# double B2
# double N
# double nr
# int S
@cache
def brsudickson(B1, B2, N, nr, S):
    sum_ = 0
    f = eulerphi(S) // 2
    i = 1
    while i <= S // 2:
        if gcd(i, S) == 1:
            sum_ += brentsuyama(B1, B2, N, nr * (gcd(i - 1, S) + gcd(i + 1, S) - 4) // 2)
        i += 1
    return sum_ / float(f)


# ret double
# double B1
# double B2
# double N
# double nr
# int S
@cache
def brsupower(B1, B2, N, nr, S):
    # int i, f;
    # double sum;
    sum_ = 0
    f = eulerphi (S);
    i = 1
    while i < S:
        if gcd (i, S) == 1:
            sum_ += brentsuyama (B1, B2, N, nr * (gcd (i - 1, S) - 2))
        i += 1

    return sum_ / float(f);


# ret long
# unsigned long n
@cache
def eulerphi(n):
    phi = 1
    p = 2

    while p * p <= n:
        if n % p == 0:
            phi *= p - 1
            n //= p
            while n % p == 0:
                phi *= p
                n //= p
        if p == 2:
            p -= 1
        p += 2
    # now n is prime
    return phi if n == 1 else (phi * (n - 1))


# ret unsigned long
# unsigned long a
# unsigned long b
def gcd(a, b):
    while b != 0:
        t = a % b
        a = b
        b = t
    return a


def exp2(n):
    pow(2, n)


# Compute probability for primes p == r (mod m)
# ret double
# double B1
# double B2
# double N
# double nr
# int S
# unsigned long r
# unsigned long m
def pm1prob_rm(B1, B2, N, nr, S, r, m):
    # unsigned long cof
    smoothness = 1.2269688
    # unsigned long p
    cof = m
    p = 2
    while p < 100:
        if cof % p == 0:  # For each prime in m
            # unsigned long cof_r, k, i;
            # Divisibility by i is determined by r and m. We need to
            # adjust the smoothness parameter. In P-1, we had estimated the
            # expected value for the exponent of p as p/(p-1)^2. Undo that.
            smoothness -= float(p) / ((p-1)*(p-1)) * math.log(float(p))
            # The expected value for the exponent of this prime is k s.t.
            # p^k || r, plus 1/(p-1) if p^k || m as well
            cof_r = gcd (r - 1, m)
            k = 0
            while cof_r % p == 0:
                cof_r //= p
                k += 1
            smoothness += k * math.log(float(p))

            cof_r = m
            i = 0
            while cof_r % p == 0:
                cof_r //= p
                i += 1

            if i == k:
                smoothness += (1./(p - 1.) * math.log(float(p)))

            while cof % p == 0:
                cof //= p
            logging.getLogger().debug(f"pm1prob_rm: p = {p}, k = {k}, i = {i}, new smoothness = {smoothness}")
        p += 1

    return prob(B1, B2, N, nr, S, smoothness)


# def memory_use(dF, sp_num, Ftreelvl, modulus):
#     mem = 9.0  # F:1, T:3*2, invF:1, G:1
#     mem += float(Ftreelvl)
#     mem *= float(dF)
#     mem += 2.0 * list_mul_mem(dF)  # Also in T
#     # estimated memory for list_mult_n /
#     # wrap-case in PrerevertDivision respectively
#     mem += (24.0 + 1.0) * float(sp_num if sp_num else dF)
#     mem *= float(mpz_size(modulus.orig_modulus)) * sizeof(mp_limb_t) + sizeof(mpz_t)
#
#     if sp_num:
#         mem += (
#             # peak malloc in ecm_ntt.c
#             4.0 * dF * sp_num * sizeof(sp_t)
#             # mpzspv_normalise
#             + MPZSPV_NORMALISE_STRIDE * (float(sp_num) * sizeof(sp_t) + 6.0 * sizeof(sp_t) + sizeof(float))
#             # sp_F, sp_invF
#             + (1.0 + 2.0) * dF * sp_num * sizeof(sp_t)
#         )
#
#     return mem
#
# # Helper functions
# def list_mul_mem(dF):
#     # Implement the list_mul_mem function
#     pass
#
# def mpz_size(modulus):
#     # Implement the mpz_size function
#     gmpy2.bit_length()
#
# def sizeof(obj):
#     # Implement the sizeof function
#     pass


# ret int
def bestD(B2min, B2, po2, use_ntt, maxmem, treefile, modulus, k=0):
    # the following list contains successive values of b with
    # increasing values of eulerphi(b). It was generated by the following
    # Maple program:
    # l := [[1,1]]:
    # for b from 12 by 6 do
    #    d:=numtheory[phi](b)/2;
    #    while d <= l[nops(l)][2] do l:=subsop(nops(l)=NULL, l) od;
    #    n := nops(l);
    #    if b>1.1*l[n][1] then l := [op(l), [b,d]]; lprint(l) fi;
    # od:
    N = 109
    l = [12, 18, 30, 42, 60, 90, 120, 150, 210, 240, 270, 330, 420, 510, 630, 840, 1050, 1260, 1470, 1680, 1890, 2310, 2730, 3150, 3570, 3990, 4620, 5460, 6090, 6930, 8190, 9240, 10920, 12180, 13860, 16170, 18480, 20790, 23100, 30030, 34650, 39270, 43890, 48510, 60060, 66990, 78540, 90090, 99330, 120120, 133980, 150150, 180180, 210210, 240240, 270270, 300300, 334950, 371280, 420420, 510510, 570570, 600600, 630630, 746130, 870870, 1021020, 1141140, 1291290, 1531530, 1711710, 1891890, 2081310, 2312310, 2552550, 2852850, 3183180, 3573570, 3993990, 4594590, 5105100, 5705700, 6322470, 7147140, 7987980, 8978970, 10210200, 11741730, 13123110, 14804790, 16546530, 19399380, 21411390, 23993970, 26816790, 29609580, 33093060, 36606570, 40330290, 44414370, 49639590, 54624570, 60090030, 67897830, 77597520, 87297210, 96996900, 107056950, 118107990]
    Npo2 = 23
    lpo2 = [12, 30, 60, 120, 240, 510, 1020, 2310,
            4620, 9240, 19110, 39270, 79170, 158340,
            324870, 690690, 1345890, 2852850, 5705700,
            11741730, 23130030, 48498450, 96996900]

    d1 = 0
    d2 = 0
    dF = 0

    if B2 < B2min:
        # No stage 2. Set relevant parameters to 0. Leave B2, B2min the same
        return dF, k, B2

    # /* Look for largest dF we can use while satisfying the maxmem parameter */
    maxN = Npo2 if po2 else N
    if maxmem != 0.:
        raise ValueError("maxmem not supported at this time")
    # if maxmem != 0.:
    #     i = 0
    #     while i < maxN; i++:
    #       lg_dF = 0
    #       sp_num = 0
    #
    #       d1 = lpo2[i] if po2 else l[i]
    #       phid = eulerphi (d1) // 2
    #       dF = 1 << math.ceil(math.log2(phid)) if po2 else phid
    #       lg_dF = math.ceil(math.log2(dF))
    #
    #       if use_ntt:
    #           sp_num = (2 * gmpy2.num_digits(modulus, 2) + lg_dF) / SP_NUMB_BITS + 4;
    #
    #       memory = memory_use (dF, sp_num, (treefile) ? 0 : lg_dF, modulus);
    #       outputf (OUTPUT_DEVVERBOSE,
    #                "Estimated mem for dF = %.0d, sp_num = %d: %.0f\n",
    #                dF, sp_num, memory);
    #       if (memory > maxmem)
    #         break;
    #     }
    #   maxN = i;
    # }

    i = 0
    while i < maxN:
        d1 = lpo2[i] if po2 else l[i]
        phid = eulerphi(d1) // 2
        dF = 1 << math.ceil(math.log2(phid)) if po2 else phid
        # /* Look for smallest prime < 25 that does not divide d1 */
        # /* The caller can force d2 = 1 by setting root_params->d2 != 0 */
        d2 = 1
        # not sure when this happens
        if True:  # root_params.d2 == 0:
            d2 = 3
            while d2 < 25:
                d2 += 2
                if d2 % 3 == 0:
                    continue
                if d1 % d2 > 0:
                    break

        if d2 >= 25 or d2 - 1 > dF:
            d2 = 1

        # #if 0
        #       /* The code to init roots of G can handle negative i0 now. */
        #       if (d2 > 1 && mpz_cmp_ui (B2min, (d1 - 1) * d2 - d1) <= 0)
        #         d2 = 1; /* Would make i0 < 0 */
        # #endif

        i0 = gmpy2.mpz(d1 - 1)  # mpz_set_ui (i0, d1 - 1);
        i0 = i0 * gmpy2.mpz(d2)  # mpz_mul_ui (i0, i0, d2);
        j = gmpy2.mpz(B2)  # mpz_set (j, B2);
        i1 = j + i0  # mpz_add (i1, j, i0); /* i1 = B2 + (d1 - 1) * d2 */
        j = gmpy2.mpz(B2min)  # mpz_set (j, B2min);
        i0 = j - i0  # mpz_sub (i0, j, i0); /* i0 = B2min - (d1 - 1) * d2 */
        i0 = gmpy2.c_div(i0, d1)  # mpz_cdiv_q_ui (i0, i0, d1); /* i0 = ceil ((B2min - (d1 - 1) * d2) / d1) */
        i1 = gmpy2.f_div(i1, d1)  # mpz_fdiv_q_ui (i1, i1, d1); /* i1 = floor ((B2 + (d1 - 1) * d2) / d1) */

        # How many roots of G will we need ?
        j = i1 - i0  # mpz_sub (j, i1, i0);
        j = j + 1  # mpz_add_ui (j, j, 1);

        # Integer multiples of d2 are skipped (if d2 > 1)
        if d2 > 1:
            t = gmpy2.f_div(i1, d2)  # mpz_fdiv_q_ui (t, i1, d2)
            j = j - t  # mpz_sub (j, j, t)
            t = gmpy2.f_div(i0, d2)  # mpz_fdiv_q_ui (t, i0, d2)
            j = j + t  # mpz_add (j, j, t); /* j -= floor (i1 / d2) - floor (i0 / d2) */

        # How many blocks will we need ? Divide lines by dF, rounding up
        j = gmpy2.c_div(j, dF)  # mpz_cdiv_q_ui (j, j, dF);

        if (k != ECM_DEFAULT_K and j <= k) or (k == ECM_DEFAULT_K and j <= (6 if po2 else 2)):
            break
        i += 1

    if i == maxN:
        if k != ECM_DEFAULT_K:
            # The user asked for a specific k and we couldn't satisfy the condition. Nothing we can do ...
            print("Error: too large step 2 bound, increase -k", file=sys.stderr)
            return ECM_ERROR
        elif gmpy2.num_digits(j, 2) > 32:
            # Can't fit the number of blocks in an unsigned long. Nothing we can do ...
            print("Error: stage 2 interval too large, cannot generate suitable parameters.\nTry a smaller B2 value.",
                  file=sys.stderr)
            return ECM_ERROR
        if maxN == 0:
            # We can't do a stage 2 at all with the memory the user allowed.
            # Nothing we can do ...
            print("Error: stage 2 not possible with memory allowed by -maxmem.", file= sys.stderr)
            return ECM_ERROR
        # else: We can fit the number of blocks into an unsigned int, so we use
        # it. This may be a very large value for huge B2-B2min, the user
        # is going to notice sooner or later

    # If the user specified a number of blocks, we'll use that no matter what.
    # Since j may be smaller than k, this may increase the B2 limit
    if k == ECM_DEFAULT_K:
        k = int(j)  # mpz_get_ui (j)

    # Now that we have the number of blocks, compute real i1. There will be
    # k * dF roots of G computed, starting at i0, skipping all that are not
    # coprime to d2. While d2 is prime, that means: are not multiples of d2.
    # Hence we want i1 so that
    #   i1 - floor(i1 / d2) - i0 + ceil((i0 / d2) == k * dF
    #   i1 - floor(i1 / d2) == k * dF + i0 - ceil((i0 / d2)

    j = gmpy2.mpz(k)  # mpz_set_ui (j, k);
    j = j * dF   # mpz_mul_ui (j, j, dF);
    if d2 == 1:
        i1 = i0 + j  # mpz_add (i1, i0, j);
        i1 = i1 - 1  # mpz_sub_ui (i1, i1, 1);
    else:
        j = j + i0  # mpz_add (j, j, i0);
        t = gmpy2.c_div(i0, d2)  # mpz_cdiv_q_ui (t, i0, d2);
        j = j - t  # mpz_sub (j, j, t); /* j = k * dF + i0 - ceil((i0 / d2) */
        j, t = gmpy2.f_divmod(j, d2 - 1)  # mpz_fdiv_qr_ui (j, t, j, d2 - 1);
        j = j * d2  # mpz_mul_ui (j, j, d2);
        i1 = j + t  # mpz_add (i1, j, t);

    # root_params->d1 = d1;
    # root_params->d2 = d2;
    # mpz_set (root_params->i0, i0);
    # *finaldF = dF;
    # *finalk = k;

    # We want B2' the largest integer that satisfies
    # i1 = floor ((B2' + (d1 - 1) * d2) / d1)
    #    = floor ((B2'-d2)/d1) + d2
    # i1 - d2 = floor ((B2'-d2)/d1)
    # (B2'-d2)/d1 < i1-d2+1
    # B2'-d2 < (i1-d2+1) * d1
    # B2' < (i1-d2+1) * d1 + d2
    # B2' = (i1-d2+1) * d1 + d2 - 1

    i1 = i1 - (d2 - 1)  # mpz_sub_ui (i1, i1, d2 - 1);
    B2 = i1 * d1  # mpz_mul_ui (B2, i1, d1);
    B2 = B2 + (d2 - 1)  # mpz_add_ui (B2, B2, d2 - 1);

    # need to return d1, d2, i0?
    return dF, k, B2


if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("Usage: rho <B1> <B2> <N> <nr> <S> [<r> <m>]\n\n")
        print(" Calculate the probability of ECM/P-1 finding a factor near N\n"
              " with B1/B2, evaluating nr random distinct points in stage 2,\n"
              " with a degree -S Dickson polynomial (if S < 0) or\n"
              " S'th power as the Brent-Suyama function\n")
        print(" <B1>        B1 limit.")
        print(" <B2>        B2 limit.")
        print(" <N>         N of similiar size, or number of bits in factor (if < 50).")
        print(" <nr>        Number of random points evaluated in stage 2.")
        print(" <S>         Degree of Brent-Suyama polynomial in stage 2.")
        print(" [<r> <m>]   Limit P-1 to primes p == r (mod m).")
        exit(1)
    if len(sys.argv) < 6:
        print("Need 5 or 7 arguments: B1 B2 N nr S [r m]")
        exit(1)

    # logging to stdout only, useful for turning off debug for api calls
    root = logging.getLogger()
    root.setLevel(logging.DEBUG)
    handler = logging.StreamHandler(sys.stdout)
    handler.setLevel(logging.DEBUG)
    formatter = logging.Formatter("%(message)s")
    handler.setFormatter(formatter)
    root.addHandler(handler)

    B1 = float(sys.argv[1])
    B2 = float(sys.argv[2])
    N = float(sys.argv[3])
    nr = float(sys.argv[4])
    S = int(sys.argv[5])
    r = 0
    m = 1
    if len(sys.argv) > 7:
        r = int(sys.argv[6])
        m = int(sys.argv[7])

    rhoinit(256, 10)
    if N < 50.:
        sum_ = _ecmprob(B1, B2, exp2(N), nr, S)
        sum_ += 4. * _ecmprob(B1, B2, 3./2. * exp2(N), nr, S)
        sum_ += _ecmprob(B1, B2, 2. * exp2(N), nr, S)
        sum_ *= 1./6.
        print("ECM: {sum_:.16f}")

        sum_ = pm1prob_rm(B1, B2, exp2(N), nr, S, r, m)
        sum_ += 4. * pm1prob_rm (B1, B2, 3./2. * exp2(N), nr, S, r, m)
        sum_ += pm1prob_rm (B1, B2, 2. * exp2(N), nr, S, r, m)
        sum_ *= 1./6.
        print(f"P-1: {sum_:.16f}")
    else:
        print(f"ECM: {_ecmprob(B1, B2, N, nr, S):.16f}", )
        print(f"P-1: {pm1prob_rm(B1, B2, N, nr, S, r, m):.16f}")
    rhoinit(0, 0)  # probably not necessary in python lol
    exit(0)



