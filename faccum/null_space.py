import mpmath


def null_space(matrix, normalized=True):
    """
    Returns (normalized) null space of a give input matrix

    :param matrix: The matrix whose null space is to be determined
    :param normalized: Return the normalized null space. Default True:
    :return: The (normalized) null space of the given matrix
    """
    u, s, v = mpmath.svd(matrix)
    column = matrix.cols - 1
    solution = v.T.column(column)

    if normalized:
        return solution / sum(solution)

    return solution
