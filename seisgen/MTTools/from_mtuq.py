# Functions copied from MTUQ project (https://github.com/uafgeotools/mtuq)

import numpy as np

def to_mij(rho, v, w, kappa, sigma, h):
    """ Converts from lune parameters to moment tensor parameters
    (up-south-east convention)
    """
    kR3 = np.sqrt(3.)
    k2R6 = 2. * np.sqrt(6.)
    k2R3 = 2. * np.sqrt(3.)
    k4R6 = 4. * np.sqrt(6.)
    k8R6 = 8. * np.sqrt(6.)

    m0 = rho / np.sqrt(2.)

    delta, gamma = to_delta_gamma(v, w)
    beta = 90. - delta

    gamma = np.deg2rad(gamma)
    beta = np.deg2rad(90. - delta)
    kappa = np.deg2rad(kappa)
    sigma = np.deg2rad(sigma)
    theta = np.arccos(h)

    Cb = np.cos(beta)
    Cg = np.cos(gamma)
    Cs = np.cos(sigma)
    Ct = np.cos(theta)
    Ck = np.cos(kappa)
    C2k = np.cos(2.0 * kappa)
    C2s = np.cos(2.0 * sigma)
    C2t = np.cos(2.0 * theta)

    Sb = np.sin(beta)
    Sg = np.sin(gamma)
    Ss = np.sin(sigma)
    St = np.sin(theta)
    Sk = np.sin(kappa)
    S2k = np.sin(2.0 * kappa)
    S2s = np.sin(2.0 * sigma)
    S2t = np.sin(2.0 * theta)

    mt0 = m0 * (1. / 12.) * \
          (k4R6 * Cb + Sb * (kR3 * Sg * (-1. - 3. * C2t + 6. * C2s * St * St) + 12. * Cg * S2t * Ss))

    mt1 = m0 * (1. / 24.) * \
          (k8R6 * Cb + Sb * (-24. * Cg * (Cs * St * S2k + S2t * Sk * Sk * Ss) + kR3 * Sg * \
                             ((1. + 3. * C2k) * (
                                         1. - 3. * C2s) + 12. * C2t * Cs * Cs * Sk * Sk - 12. * Ct * S2k * S2s)))

    mt2 = m0 * (1. / 6.) * \
          (k2R6 * Cb + Sb * (kR3 * Ct * Ct * Ck * Ck * (1. + 3. * C2s) * Sg - k2R3 * Ck * Ck * Sg * St * St +
                             kR3 * (1. - 3. * C2s) * Sg * Sk * Sk + 6. * Cg * Cs * St * S2k +
                             3. * Ct * (-4. * Cg * Ck * Ck * St * Ss + kR3 * Sg * S2k * S2s)))

    mt3 = m0 * (-1. / 2.) * Sb * (k2R3 * Cs * Sg * St * (Ct * Cs * Sk - Ck * Ss) +
                                  2. * Cg * (Ct * Ck * Cs + C2t * Sk * Ss))

    mt4 = -m0 * (1. / 2.) * Sb * (Ck * (kR3 * Cs * Cs * Sg * S2t + 2. * Cg * C2t * Ss) +
                                  Sk * (-2. * Cg * Ct * Cs + kR3 * Sg * St * S2s))

    mt5 = -m0 * (1. / 8.) * Sb * (4. * Cg * (2. * C2k * Cs * St + S2t * S2k * Ss) +
                                  kR3 * Sg * ((1. - 2. * C2t * Cs * Cs - 3. * C2s) * S2k + 4. * Ct * C2k * S2s))

    if type(mt0) is np.ndarray:
        return np.column_stack([mt0, mt1, mt2, mt3, mt4, mt5])
    else:
        return np.array([mt0, mt1, mt2, mt3, mt4, mt5])


def to_delta_gamma(v, w):
    """ Converts from Tape2015 parameters to lune coordinates
    """
    return to_delta(w), to_gamma(v)


def to_gamma(v):
    """ Converts from Tape2015 parameter v to lune longitude
    """
    gamma = (1. / 3.) * np.arcsin(3. * v)
    return np.rad2deg(gamma)


def to_delta(w):
    """ Converts from Tape2015 parameter w to lune latitude
    """
    beta0 = np.linspace(0, np.pi, 100)
    u0 = 0.75 * beta0 - 0.5 * np.sin(2. * beta0) + 0.0625 * np.sin(4. * beta0)
    beta = np.interp(3. * np.pi / 8. - w, u0, beta0)
    delta = np.rad2deg(np.pi / 2. - beta)
    return delta


def to_v_w(delta, gamma):
    """ Converts from lune coordinates to Tape2015 parameters
    """
    return to_v(gamma), to_w(delta)


def to_v(gamma):
    """ Converts from lune longitude to Tape2015 parameter v
    """
    v = (1. / 3.) * np.sin(3. * np.deg2rad(gamma))
    return v


def to_w(delta):
    """ Converts from lune latitude to Tape2015 parameter w
    """
    beta = np.deg2rad(90. - delta)
    u = (0.75 * beta - 0.5 * np.sin(2. * beta) + 0.0625 * np.sin(4. * beta))
    w = 3. * np.pi / 8. - u
    return w


def to_M0(Mw):
    """ Converts from moment magnitude to scalar moment
    """
    return 10. ** (1.5 * float(Mw) + 9.1)


def to_rho(Mw):
    """ Converts from moment magnitude to Tape2012 magnitude parameter
    """
    return to_M0(Mw) * np.sqrt(2.)


def to_Mw(rho):
    """ Converts from Tape2012 magnitude parameter to moment magnitude
    """
    return ((np.log10(rho / np.sqrt(2)) - 9.1) / 1.5)

