# -------------------------------------------------------------------
# Acquire the waveform from the
# pre-stored 3D receiver-side strain Greens function database
#
# Author: Liang Ding
# Email: myliang.ding@mail.utoronto.ca
# -------------------------------------------------------------------


from seisgen.DSGTMgr import DSGTMgr
from obspy.core.util.attribdict import AttribDict
import numpy as np


def acquire_waveform():
    '''Acquire the waveform '''

    sgt_database_folder = "path-to-the-SGT-database/"
    model3D_folder="path-to-the-model-files(the *.bin files)/"
    point_info_file= "path-to-the-HDF5-file-storing-the-point-location-and-interpolation-parameters(xi, eta, gamma)"

    # origin example
    origin = AttribDict({'time': '2019-07-04T18:39:44.0000Z',
                     'latitude': 35.601333,
                     'longitude': -117.597,
                     'depth_in_m': 2810.0,
                     'id': 'EVT1'})

    # station example: CI.SLA
    station = AttribDict({
        'latitude': 35.8909,
        'longitude': -117.2833,
        'network': 'CI',
        'station': 'SLA',
        'location': '',
        'id': 'SLA'})

    sgt_mgr = DSGTMgr(sgt_database_folder, model3D_folder, point_info_file)

    # moment tensor example [Mrr, Mtt, Mpp, Mrt, Mrp, Mtp]
    mt_RTP = np.array([8.6E20, -1.3E22, 1.2E22, -2.1E21, -3.2E21, 6.8E21])

    print('Reading waveform ...\n')
    st = sgt_mgr.get_waveform(station, origin, mt_RTP)
    print(st)


if __name__ == '__main__':
    acquire_waveform()