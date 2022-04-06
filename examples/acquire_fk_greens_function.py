# -------------------------------------------------------------------
# Acquire the FK-type Green's function from the SGT database.
# Author: Liang Ding
# Email: myliang.ding@mail.utoronto.ca
# -------------------------------------------------------------------


from seisgen.DSGTMgr import DSGTMgr
from obspy.core.util.attribdict import AttribDict


def acquire_fk_greens_function():
    '''Acquire FK-type Green's function '''

    sgt_database_folder = "path-to-the-SGT-database/"
    model3D_folder="path-to-the-model-files(the *.bin files)/"
    point_info_file= "path-to-the-HDF5-file-storing-the-point-location-and-interpolation-parameters(xi, eta, gamma)"

    # Path to store greens function
    greens_path = "path-to-store-the-fk-type-greens-function/greens/"

    # Store fk-type greens function to pickle file.
    pkl_greens_file_path = "path-to-store-the-fk-type-greens-function/fk_greens.pkl"

    # origin example (example 1)
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
    print('Acquiring Greens functions...\n')
    greens = sgt_mgr.get_fk_greens_function(station, origin, b_save=True, greens_path=greens_path)
    print(greens)

    import pickle
    with open(pkl_greens_file_path, "wb") as f:
        pickle.dump(greens, f)


if __name__ == '__main__':
    acquire_fk_greens_function()


