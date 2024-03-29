# -------------------------------------------------------------------
# Moment tensor
#
# refs
# Tape, W., & Tape, C. (2015). A uniform parametrization of moment tensors.
# Geophysical Journal International, 202(3), 2074–2081. https://doi.org/10.1093/gji/ggv262
#
# Author: Liang Ding
# Email: myliang.ding@mail.utoronto.ca
# -------------------------------------------------------------------


import numpy as np

def DCreate_Triple(colatitude, longitude):
    '''
    * Create the normalized triple according to the colat. and long. on the Lune.
    * Equ. (7) in Carl's paper (2015).
    :param colatitude: The colatitude in the Lune in radians (Beta)
    :param longitude:  The longitude in the Lune in radians (Gamma)
    :return: normalized eigenvalues in descent order (a triple on the LUNE) [e1>e2>e3]
    '''

    sin_beta = np.sin(colatitude)
    cos_beta = np.cos(colatitude)
    sin_gamma= np.sin(longitude)
    cos_gamma = np.cos(longitude)

    sqrt3 = 1.73205081
    sqrt2 = 1.41421356
    sqrt6 = 2.44948974

    eigenvalues = np.dot(sqrt6 * np.array([[sqrt3, -1.0, sqrt2],
                                           [0, 2.0, sqrt2],
                                           [-1.0 * sqrt3, -1, sqrt2]]),
                         np.array([sin_beta * cos_gamma, sin_beta * sin_gamma, cos_beta]).transpose())

    return eigenvalues / np.linalg.norm(eigenvalues)