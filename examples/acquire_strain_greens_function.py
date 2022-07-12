# -------------------------------------------------------------------
# Acquiring the strain Greens function
# from the 3D pre-stored receiver-side strain Greens function database
#
# Author: Liang Ding
# Email: myliang.ding@mail.utoronto.ca
# -------------------------------------------------------------------


from seisgen.DSGTMgr import DSGTMgr
from obspy.core.util.attribdict import AttribDict


def acquire_greens_function():
    '''Acquire the Greens function '''

    sgt_database_folder = "path-to-the-SGT-database/"
    model3D_folder="path-to-the-model-files(the *.bin files)/"
    point_info_file= "path-to-the-HDF5-file-storing-the-point-location-and-interpolation-parameters(xi, eta, gamma)"

    # origin example
    origin = AttribDict({'time': '2019-07-04T18:39:44.0000Z',
                     'latitude': 35.601333,
                     'longitude': -117.597,
                     'depth_in_m': 2810.0,
                     'id': 'evt1'})

    # station example: CI.SLA
    station = AttribDict({
        'latitude': 35.8909,
        'longitude': -117.2833,
        'network': 'CI',
        'station': 'SLA',
        'location': '',
        'id': 'SLA'})

    sgt_mgr = DSGTMgr(sgt_database_folder, model3D_folder, point_info_file)
    print('Reading Greens functions...\n')
    greens = sgt_mgr.get_greens_function(station, origin)
    print(greens)


if __name__ == '__main__':
    acquire_greens_function()


