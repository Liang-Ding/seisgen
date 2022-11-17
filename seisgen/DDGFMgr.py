# -------------------------------------------------------------------
# Green's Function (Displacement, DGF) database Manager
#
# Author: Liang Ding
# Email: myliang.ding@mail.utoronto.ca
# -------------------------------------------------------------------

import os.path
from seisgen.DPointMgr import DPointCloud
from seisgen.util_SPECFEM3D import get_proc_name
from seisgen.util_SPECFEM3D.ibool_reader import DEnquire_Element
from seisgen.util_SPECFEM3D.xyz_reader import DEnquire_XYZ_GLLs_Element
from seisgen.greens_function.dgf_reader import DEnquire_DGF, read_DGF_header
from seisgen.math.interp_tools import DCreate_anchors_xi_eta_gamma, DLagrange_interp_sgt, DLagrange_any3D

from obspy.core.util.attribdict import AttribDict
from obspy.core import Stream, Trace
from obspy.clients.iris import Client

import numpy as np
import time


class DDGFMgr(DPointCloud):
    '''Greens Function (Displacement, DGF) database Manager'''

    def __init__(self, dgf_database_folder, model3D_folder, point_cloud_file, DLite=False):
        '''
        :param dgf_database_folder:     The directory to the greens function database.
        :param model3D_folder:          The directory to the 3D background model.
        :param point_cloud_file:        The hdf5 file storing the information of user-selected points in the 3D model.
        '''

        self.dgf_database_folder = dgf_database_folder
        self.model3D_folder = model3D_folder
        self.idx_element = -1
        self.dt = 0
        self.NSPEC = 0
        self.nGLL_global = 0
        # initial parameters of the GF database.
        if DLite:
            pass
        else:
            self.__initial_paras()
        super().__init__(point_cloud_file)

    def __initial_paras(self):
        '''Initializing the parameters.'''
        hdf5_files = []
        for (dirpath, dirnames, filenames) in os.walk(self.dgf_database_folder):
            for file in filenames:
                if str(file).lower().endswith('.hdf5') or str(file).lower().endswith('.h5'):
                    hdf5_files.append(os.path.join(dirpath, file))

        for file in hdf5_files:
            try:
                self.dt, self.NSPEC, self.nGLL_global = read_DGF_header(file)
                break
                
            except:
                continue

        if 0 == self.dt or 0 == self.NSPEC or 0 == self.nGLL_global:
            raise Exception


    def _initial_element_frame(self):
        ''' return the gll information (index, location) at one selected element. '''
        ibool_file = os.path.join(str(self.model3D_folder), str(self.proc_name) + str("_ibool.bin"))
        self.idx_glls = DEnquire_Element(ibool_file, self.idx_element, self.NSPEC)
        x_glls, y_glls, z_glls = DEnquire_XYZ_GLLs_Element(self.model3D_folder,
                                                           self.idx_processor,
                                                           self.idx_element,
                                                           self.NSPEC)
        self.xyz_glls = np.transpose(np.vstack([x_glls, y_glls, z_glls]))


    def _initial_DGFs_N_station(self):
        '''Return the GFs between origin and station. '''
        dir_string = os.path.join(str(self.dgf_database_folder),
                                  str(self.station.network),
                                  str(self.station.station),
                                  str(self.proc_name))
        data_path = dir_string + str("_dgf_data.bin")
        hder_path = dir_string + str("_dgf_header.hdf5")
        self.dgfs = DEnquire_DGF(data_path, hder_path, self.idx_glls)


    def interp_Lagrange(self, xi, eta, gamma):
        '''Using lagrange method to interpolate the GF data.'''
        ngll_x = 3
        ngll_y = 3
        ngll_z = 3

        xi_gll, eta_gll, gamma_gll = DCreate_anchors_xi_eta_gamma(ngll_xyz=3)
        h_xi_arr, h_eta_arr, h_gamma_arr = DLagrange_any3D(xi, eta, gamma, xi_gll, eta_gll, gamma_gll)

        self.dgf_interp = DLagrange_interp_sgt(h_xi_arr, h_eta_arr, h_gamma_arr, self.dgfs,
                                               ngll_x=ngll_x, ngll_y=ngll_y, ngll_z=ngll_z)

    def set_dt(self, dt):
        '''
        Set time interval.
        Function observed for SeisClient.
        '''
        self.dt = dt
        return self


    def get_dgf(self, station, origin, b_new_origin=True, b_verbose=False):
        '''
        Get the interpolated GF between the station-origin pair.
        Unit: m/N.m

        :param station: An instance of the obspy AttribDict class. For example:
                        station = AttribDict({ 'latitude': 34.0210,
                                                'longitude': -118.287,
                                                'network': 'CI',
                                                'station': 'USC',
                                                'location': '',
                                                'id': 'USC'})

        :param origin:  An instance of the obspy AttribDict class. For example:
                        origin = Origin({'time': '2019-07-04T18:39:44.0000Z',
                                              'latitude': 35.601333,
                                              'longitude': -117.597,
                                              'depth_in_m': 2810.0,
                                              'id': 'evt11056825'})

        :param b_new_origin: True, If acquiring GFs at multiple stations for a same origin to save time.
        '''

        if b_verbose:
            t0 = time.time()

        self.station = station

        if not self.b_pointcloud_initial:
            raise Exception

        try:
            z = origin.depth_in_m
            b_depth = True
        except:
            z = origin.elevation_in_m
            b_depth = False

        if b_new_origin:
            _, _, _, _, \
            _, _, _, \
            self.idx_processor, self.element_index, \
            self.xi, self.eta, self.gamma = self.find(x=origin.latitude, y=origin.longitude, z=z, n=1, b_depth=b_depth)

            # MUST subtract 1, the element_index in the point_cloud_file starts from 1.
            self.idx_element = self.element_index - 1
            self.proc_name = get_proc_name(self.idx_processor)
            self._initial_element_frame()
        self._initial_DGFs_N_station()
        # The Lagrange interpolation to get interpolated GF.
        self.interp_Lagrange(self.xi, self.eta, self.gamma)
        if b_verbose:
            print("Station={}, time cost={} s.".format(station.id, time.time() - t0))

        return self.dgf_interp


    def get_greens_function(self, station, origin, b_new_origin=True):
        '''
        Get Greens Function between the station-origin pair.
        Unit: m/N.m

        :param station: An instance of the obspy AttribDict class. For example:
                        station = AttribDict({ 'latitude': 34.0210,
                                                'longitude': -118.287,
                                                'network': 'CI',
                                                'station': 'USC',
                                                'location': '',
                                                'id': 'USC'})

        :param origin:  An instance of the obspy AttribDict class. For example:
                        origin = Origin({'time': '2019-07-04T18:39:44.0000Z',
                                              'latitude': 35.601333,
                                              'longitude': -117.597,
                                              'depth_in_m': 2810.0,
                                              'id': 'evt11056825'})

        :param b_new_origin: True, If acquiring GFs at multiple stations for a same origin to save time.
        '''

        greens = self.get_dgf(station, origin, b_new_origin=b_new_origin)
        return self.get_greens_function_next(greens, station)


    def get_greens_function_next(self, greens, station):
        ''' The next step of get_greens_function()'''
        _, n_force, n_chans = np.shape(greens)
        element_order = ['E', 'N', 'Z']
        st = Stream()

        for i in range(n_force):
            for j in range(n_chans):
                tr = Trace(greens[:, i, j])
                tr.stats.channel = '%s%s' % (element_order[i], element_order[j])
                tr.stats.delta = self.dt
                tr.stats.sampling_rate = int(1.0 / self.dt)
                try:
                    tr.stats.network = station.network
                    tr.stats.station = station.station
                except:
                    pass

                st.append(tr)
        st.id = station.id
        return st


    def get_waveform(self, station, origin, force_enz, b_RTZ=False, b_new_origin=True):
        '''
        Return the synthetic waveform in ENZ (default) or RTZ.
        Unit: m

        :param station: An instance of the obspy AttribDict class. For example:
                        station = AttribDict({ 'latitude': 34.0210,
                                                'longitude': -118.287,
                                                'network': 'CI',
                                                'station': 'USC',
                                                'location': '',
                                                'id': 'USC'})

        :param origin:  An instance of the obspy AttribDict class. For example:
                        origin = Origin({'time': '2019-07-04T18:39:44.0000Z',
                                              'latitude': 35.601333,
                                              'longitude': -117.597,
                                              'depth_in_m': 2810.0,
                                              'id': 'evt11056825'})

        :param force: the force in ENZ [E, N, Z]
        :param b_new_origin: True, If acquiring GFs at multiple stations for a same origin to save time.
        '''

        greens = self.get_dgf(station, origin, b_new_origin=b_new_origin)
        return self.get_waveform_next(greens, station, origin, force_enz, b_RTZ)


    def get_waveform_next(self, greens, station, origin, force_enz, b_RTZ=False):
        ''' The next step of get_waveform()'''

        syn_enz = np.matmul(greens, force_enz).transpose()
        n_chans = 3
        channel_order = ['E', 'N', 'Z']
        channel_name = 'SY'

        # synthetic waveform
        st = Stream()
        for i in range(n_chans):
            tr = Trace(syn_enz[i])
            tr.stats.channel = '%s%s' % (channel_name, channel_order[i])
            tr.stats.delta = self.dt
            tr.stats.sampling_rate = int(1.0 / self.dt)
            st.append(tr)
            try:
                tr.stats.network = station.network
                tr.stats.station = station.station
            except:
                pass

        st.id = station.id

        # Rotation if needed.
        if b_RTZ:
            client = Client()
            res = client.distaz(stalat=station.latitude, stalon=station.longitude,
                            evtlat=origin.latitude, evtlon=origin.longitude)
            back_azimuth = res['backazimuth']
            distance_deg = res['distance']
            distance_m = res['distancemeters']
            st.rotate(method='NE->RT', back_azimuth=back_azimuth)
        return st
