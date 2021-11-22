# -------------------------------------------------------------------
# Tools to manage the pre-computed grids in the 3D background model.
# Author: Liang Ding
# Email: myliang.ding@mail.utoronto.ca
# -------------------------------------------------------------------

import h5py
import numpy as np

POINT_KEYS = ["latitude",
             "longitude",
             "z",
             "depth",
             "utm_x",
             "utm_y",
             "utm_z",
             "slice_index",
             "element_index",
             "xi",
             "eta",
             "gamma"]


class DPointCloud():
    ''' The pre-computed information (see POINT_KEYS) of user-selected points in the 3D background model. '''

    def __init__(self, file_path):
        '''
        :param file_path:  The HDF5 file storing the pre-computed information (see POINT_KEYS)
                           of user-selected points in the 3D background model.
        '''

        self.b_pointcloud_initial = False
        if file_path is None:
            return

        try:
            with h5py.File(file_path, 'r') as f:
                self.mesh_lat           = f[POINT_KEYS[0]][:]
                self.mesh_long          = f[POINT_KEYS[1]][:]
                self.mesh_z             = f[POINT_KEYS[2]][:]    # in meter
                self.mesh_depth         = f[POINT_KEYS[3]][:]    # in meter.
                self.mesh_utm_x         = f[POINT_KEYS[4]][:]
                self.mesh_utm_y         = f[POINT_KEYS[5]][:]
                self.mesh_utm_z         = f[POINT_KEYS[6]][:]
                self.mesh_slice_index   = f[POINT_KEYS[7]][:]
                self.mesh_element_index = f[POINT_KEYS[8]][:]
                self.mesh_xi            = f[POINT_KEYS[9]][:]
                self.mesh_eta           = f[POINT_KEYS[10]][:]
                self.mesh_gamma         = f[POINT_KEYS[11]][:]

            self.n_grid = len(self.mesh_lat)
            self.b_pointcloud_initial = True
        except:
            print("!!! Point cloud not found")
            raise Exception


    def _check(self):
        if self.b_pointcloud_initial is False:
            print("!!! Point cloud not initialized")
            raise Exception


    def find(self, x, y, z, n=1, mode='LATLONGZ', b_depth=True):
        '''
        Determine the closest N points to (x, y, z) in the point cloud and return the information.
         :param x:       Either the latitude or UTMX.
         :param y:       Either the longitude or UTMY.
         :param z:       The depth or elevation in meter. eg: depth: 2000, elevation: -1200
         :param n:       The number of point enquired.
         :return:        The information (see POINT_KEYS) of the closest N points.
        '''

        self._check()

        n = int(n)
        n_point = len(self.mesh_lat)
        distance_arr = np.zeros(n_point)

        if str(mode).upper() == str('LATLONGZ'):
            if b_depth:
                # use depth.
                for i in range(n_point):
                    _dist_H = 111.0 * 1000.0 * np.sqrt(np.square(self.mesh_lat[i] - x) + np.square(self.mesh_long[i]-y))
                    distance_arr[i] = np.sqrt(np.power(_dist_H, 2) + np.power(self.mesh_depth[i] - z, 2))
            else:
                # use elevation.
                for i in range(n_point):
                    _dist_H = 111.0 * 1000.0 * np.sqrt(np.square(self.mesh_lat[i] - x) + np.square(self.mesh_long[i]-y))
                    distance_arr[i] = np.sqrt(np.power(_dist_H, 2) + np.power(self.mesh_z[i] - z, 2))

        elif str(mode).upper() == str('UTM'):
            if b_depth:
                # use depth
                for i in range(n_point):
                    distance_arr[i] = np.sqrt(np.power(self.mesh_utm_x[i] - x, 2)
                                              + np.power(self.mesh_utm_y[i] - y, 2)
                                              + np.power(self.mesh_depth[i] - z, 2))
            else:
                # use elevation.
                for i in range(n_point):
                    distance_arr[i] = np.sqrt(np.power(self.mesh_utm_x[i] - x, 2)
                                              + np.power(self.mesh_utm_y[i] - y, 2)
                                              + np.power(self.mesh_z[i] - z, 2))
        else:
            print("!!! Undefined mode!")
            raise NotImplementedError

        if 1 == n:
            idx = np.argmin(distance_arr)
        else:
            idx = np.argpartition(distance_arr, n)[:n]

        return self.mesh_lat[idx], self.mesh_long[idx], self.mesh_z[idx], self.mesh_depth[idx],\
               self.mesh_utm_x[idx], self.mesh_utm_y[idx], self.mesh_utm_z[idx], \
               self.mesh_slice_index[idx], self.mesh_element_index[idx], \
               self.mesh_xi[idx], self.mesh_eta[idx], self.mesh_gamma[idx]

