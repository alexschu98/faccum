import mpmath


def matrix_exp(matrix):
    """
    Calculate the exponential of a given matrix. Not always stable, gets confused with complex values sometimes.
    TODO: make stable and error resilient

    :param matrix: The matrix whose exponential shall be calculated
    :return: The matrix exponential of the given matrix
    """
    eigvals, eigvecs = mpmath.eig(matrix)
    return eigvecs * mpmath.diag([mpmath.exp(e) for e in eigvals]) * mpmath.inverse(eigvecs)
