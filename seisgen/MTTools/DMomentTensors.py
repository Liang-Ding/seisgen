# -------------------------------------------------------------------
# Full moment tensor
#
# refs
# Tape, W., & Tape, C. (2015). A uniform parametrization of moment tensors.
# Geophysical Journal International, 202(3), 2074–2081.
#
# Author: Liang Ding
# Email: myliang.ding@mail.utoronto.ca
# -------------------------------------------------------------------

from seisgen.MTTools.from_mtuq import to_mij, to_v_w
from seisgen.seismic.DSyn import RTP_to_DENZ
import numpy as np


def DMT_enz(strike, dip, rake, colatitude, longitude, b_ENZ=True):
    '''
    Calculate the FULL moment tensor for the given strike, dip, rake and the colatitude, longitude in Lune.

    :param strike:  Strike in degree, ranging [0, 180]
    :param dip:     Dip in degree, ranging [0, 90]
    :param rake:    Slip (rake) in degree, ranging [-90, 90]
    :param colatitude: The colatitude in degree in the Lune diagram, (Beta=no.deg2rad(colatitude), [0, PI])
    :param longitude:  The longitude in degree in the Lune diagram(Gamma=np.deg2rad(longitude), [-PI/6, PI/6])
    :return: The full moment tensor in ENZ or USE(RTP).
    '''

    rho = 1
    v, w = to_v_w(colatitude-90, longitude)
    mt = to_mij(rho, v, w, kappa=strike, sigma=rake, h=np.cos(dip))
    if b_ENZ:
        return RTP_to_DENZ(mt)
    else:
        return mt
