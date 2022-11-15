import math
import mpmath

import matrix_exp
import null_space


def generalized_factorial_cumulants(wzMatrix, t, s0, mmax, dzdigts=5):
    """ Calculate the generalized factorial cumulants for a system at a given timepoint.

    :param wzMatrix: The generator of the system's dynamics
    :param t: Time point
    :param s0: Shift of zeros for generalized factorial cumulants
    :param mmax: Maximal cumulant order
    :param dzdigts: Grid precision exponent. Defaults to 5.
    :return:
    """
    assert mmax >= 1
    assert dzdigts >= 1
    mmax, dzdigts = int(mmax), int(dzdigts)

    # grid step size
    dz = 1 / mpmath.mpf(math.prod(dzdigts * [10]))

    # get stationary occupation
    p0 = null_space.null_space(wzMatrix(1))

    # get function values on mesh
    w_list = [wzMatrix((i / 2.0) * dz + s0) for i in range(-mmax, mmax + 1)]
    m_list = [mpmath.log(sum(mpmath.expm(t * w) * p0)) for w in w_list]

    # construct derivatives
    sol = [(-1) ** (m - 1) / (dz ** m) * sum(
        [(-1) ** i * mpmath.binomial(m, i) * m_list[mmax + m - 2 * i] for i in range(0, m + 1)]) for m in
           range(1, mmax + 1)]

    return sol


def generalized_factorial_cumulants2(wzMatrix, t, s0, mmax, dzdigts=5):
    """ Calculate the generalized factorial cumulants for a system at a given timepoint.
    Uses the unstable fast_matrix_exponential method so use with caution.

    :param wzMatrix: The generator of the system's dynamics
    :param t: Time point
    :param s0: Shift of zeros for generalized factorial cumulants
    :param mmax: Maximal cumulant order
    :param dzdigts: Grid precision exponent. Defaults to 5.
    :return:
    """
    assert mmax >= 1
    assert dzdigts >= 1
    mmax, dzdigts = int(mmax), int(dzdigts)

    # grid step size
    dz = 1 / mpmath.mpf(math.prod(dzdigts * [10]))

    # get stationary occupation
    p0 = null_space.null_space(wzMatrix(1))

    # get function values on mesh
    w_list = [wzMatrix((i / 2.0) * dz + s0) for i in range(-mmax, mmax + 1)]
    m_list = [mpmath.log(sum(matrix_exp.matrix_exp(t * w) * p0)) for w in w_list]

    # construct derivatives
    sol = [(-1) ** (m - 1) / (dz ** m) * sum(
        [(-1) ** i * mpmath.binomial(m, i) * m_list[mmax + m - 2 * i] for i in range(0, m + 1)]) for m in
           range(1, mmax + 1)]

    return sol
