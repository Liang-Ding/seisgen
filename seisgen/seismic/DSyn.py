# -------------------------------------------------------------------
# Tools to operate moment tensor and create synthetic waveforms
#
# Author: Liang Ding
# Email: myliang.ding@mail.utoronto.ca
# -------------------------------------------------------------------

import numpy as np
from obspy.core.stream import Stream
from obspy.core.trace import Trace

def RTP_to_DENZ(mt):
    ''' MT from RTP to ENZ   '''
    new_mt = np.zeros_like(mt)
    new_mt[0] = mt[2]
    new_mt[1] = mt[1]
    new_mt[2] = mt[0]
    new_mt[3] = -1.0 * mt[5]
    new_mt[4] = mt[4]
    new_mt[5] = -1.0 * mt[3]
    return new_mt


def DSyn(mt, sgt, element, b_GF=False):
    '''
    :param mt:      The moment tensor (MT) in ENZ
    :param SGT:     The strain Greens tensor.
                    The force order: N-E-Z
    :param element: The string of MT component.
    :param b_GF:    Calculate the Greens function.
    :return:
    '''

    new_mt = mt.copy()
    if not b_GF:
        new_mt[3:] = 2.0 * new_mt[3:]

    n_force = 3
    n_element = 6
    stream = Stream()
    channels = ['N', 'E', 'Z']
    for i in range(n_force):
        trace = Trace(np.dot(sgt[:, i, :], new_mt))
        trace.stats.channel = str(element)+str(channels[i])
        stream.append(trace)
    return stream



