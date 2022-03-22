# -------------------------------------------------------------------
# Strain Green's Tensor (SGT) database Manager
#
# Author: Liang Ding
# Email: myliang.ding@mail.utoronto.ca
# -------------------------------------------------------------------

import os.path
from seisgen.DPointMgr import DPointCloud
from seisgen.seismic.DSyn import DSyn, RTP_to_DENZ
from seisgen.util_SPECFEM3D import get_proc_name
from seisgen.util_SPECFEM3D.ibool_reader import DEnquire_Element
from seisgen.util_SPECFEM3D.xyz_reader import DEnquire_XYZ_GLLs_Element
from seisgen.sgt.sgt_reader import DEnquire_SGT, read_header_info
from seisgen.math.interp_tools import DCreate_anchors_xi_eta_gamma, DLagrange_interp_sgt, DLagrange_any3D
from obspy.core import Stream
from obspy.geodetics.base import gps2dist_azimuth
import numpy as np

import time

MT_ELEMENTS = [
    'Mrr',
    'Mtt',
    'Mpp',
    'Mrt',
    'Mrp',
    'Mtp'
]

class DSGTMgr(DPointCloud):
    '''Strain Green's Tensor (SGT) database Manager'''

    def __init__(self, sgt_database_folder, model3D_folder, point_cloud_file):
        '''
        :param sgt_database_folder:     The directory to the SGT database.
        :param model3D_folder:          The directory to the 3D background model.
        :param point_cloud_file:        The hdf5 file storing the information of user-selected points in the 3D model.
        '''

        self.sgt_database_folder = sgt_database_folder
        self.model3D_folder      = model3D_folder
        self.idx_element         = -1
        self.dt                  = 0
        self.NSPEC               = 0
        # initial parameters of the SGT database.
        self.__initial_paras()
        super().__init__(point_cloud_file)

    def __initial_paras(self):
        '''Function to initial parameters.'''
        hdf5_files = []
        for (dirpath, dirnames, filenames) in os.walk(self.sgt_database_folder):
            for file in filenames:
                if str(file).lower().endswith('.hdf5') or str(file).lower().endswith('.h5'):
                    hdf5_files.append(os.path.join(dirpath, file))

        # Read the Header file to get the information of SGT database.
        for file in hdf5_files:
            try:
                self.dt, self.NSPEC = read_header_info(file)
                break
            except:
                continue

        if 0 == self.dt or 0 == self.NSPEC:
            raise Exception

    def _initial_element_frame(self):
        ''' return the gll information (index, location) at one selected element. '''
        ibool_file = os.path.join(str(self.model3D_folder), str(self.proc_name)+str("_ibool.bin"))
        self.idx_glls = DEnquire_Element(ibool_file,  self.idx_element, self.NSPEC)

        x_glls, y_glls, z_glls = DEnquire_XYZ_GLLs_Element(self.model3D_folder,
                                                           self.idx_processor,
                                                           self.idx_element,
                                                           self.NSPEC)
        self.xyz_glls = np.transpose(np.vstack([x_glls, y_glls, z_glls]))

    def _initial_SGTs_N_station(self):
        '''Return the SGT between origin and station. '''
        dir_string = os.path.join(str(self.sgt_database_folder),
                                     # str(self.station.network),
                                     str(self.station.station),
                                     str(self.proc_name))
        sgt_data_path = dir_string + str("_sgt_data.bin")
        sgt_hder_path = dir_string + str("_header.hdf5")
        self.sgts = DEnquire_SGT(sgt_data_path, sgt_hder_path, self.idx_glls)

    def interp_sgt_Lagrange(self, xi, eta, gamma):
        '''Using lagrange method to interpolate the SGT data.'''
        ngll_x = 3
        ngll_y = 3
        ngll_z = 3

        xi_gll, eta_gll, gamma_gll = DCreate_anchors_xi_eta_gamma(ngll_xyz=3)
        h_xi_arr, h_eta_arr, h_gamma_arr = DLagrange_any3D(xi, eta, gamma, xi_gll, eta_gll, gamma_gll)

        self.sgt_interp = DLagrange_interp_sgt(h_xi_arr, h_eta_arr, h_gamma_arr, self.sgts,
                                          ngll_x=ngll_x, ngll_y=ngll_y, ngll_z=ngll_z)


    def get_sgt(self, station, origin, b_new_origin=True, b_verbose=False):
        '''
        Get the interpolated SGT between the station-origin pair.
        :param station: instance of the Station class in MTUQ.
        :param origin:  the dict of origin.
                        eg:  origin = Origin({'time': '2019-07-04T18:39:44.0000Z',
                                              'latitude': 35.601333,
                                              'longitude': -117.597,
                                              'depth_in_m': 2810.0,
                                              'id': 'evt11056825'})

        :param b_new_origin: Accelerating the extraction of SGT data for multiple stations with the same origin.
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

        self._initial_SGTs_N_station()
        # The Lagrange interpolation to get interpolated SGT.
        self.interp_sgt_Lagrange(self.xi, self.eta, self.gamma)

        if b_verbose:
            print("Station={}, time cost={} s.".format(station.id, time.time()-t0))

        return self.sgt_interp


    def get_greens_function(self, station, origin, b_new_origin=True):
        '''
        Get Green's Function between the station-origin pair.
        :param station: instance of the Station class in MTUQ.
        :param origin:  the dict of origin.
                        eg:  origin = Origin({'time': '2019-07-04T18:39:44.0000Z',
                                              'latitude': 35.601333,
                                              'longitude': -117.597,
                                              'depth_in_m': 2810.0,
                                              'id': 'evt11056825'})

        :param b_new_origin: Accelerating the extraction of SGT data for multiple stations with the same origin.
        '''

        sgt = self.get_sgt(station, origin, b_new_origin=b_new_origin)
        (distance, azimuth, back_azimuth) = gps2dist_azimuth(station.latitude, station.longitude,
                                                             origin.latitude, origin.longitude)
        stream = self.SGT2GF(sgt, back_azimuth)
        stream.id = station.id
        return stream


    def get_waveform(self, station, origin, mt_RTP, b_new_origin=True):
        '''
        Return the synthetic 3-C waveform in RTZ.
        :param station: instance of the Station class in MTUQ.
        :param origin:  the dict of origin.
                        eg:  origin = Origin({'time': '2019-07-04T18:39:44.0000Z',
                                              'latitude': 35.601333,
                                              'longitude': -117.597,
                                              'depth_in_m': 2810.0,
                                              'id': 'evt11056825'})

        :param mt: the moment tensor in RTP [Mrr, Mtt, Mpp, Mrt, Mrp, Mtp]
        :param b_new_origin: Accelerating the extraction of SGT data for multiple stations with the same origin.
        '''

        sgt = self.get_sgt(station, origin, b_new_origin=b_new_origin)
        (distance, azimuth, back_azimuth) = gps2dist_azimuth(station.latitude, station.longitude,
                                                             origin.latitude, origin.longitude)
        # synthetic
        element='SY'
        mt_enz = RTP_to_DENZ(mt_RTP)
        _st = DSyn(mt_enz, sgt, element)
        _st.rotate(method='NE->RT', back_azimuth=back_azimuth)
        return _st


    def SGT2GF(self, sgt, ba):
        '''Generate Green's Function (GF) from Strain Green's Tensor (SGT) database.'''
        MTs_rtp = np.identity(6)
        stream = Stream()

        for i, mt_rtp in enumerate(MTs_rtp):
            mt_enz = RTP_to_DENZ(mt_rtp)
            # Synthetic waveform in ENZ
            _st = DSyn(mt_enz, sgt, MT_ELEMENTS[i], b_GF=True)
            # Rotation (ENZ => RTZ)
            _st.rotate(method='NE->RT', back_azimuth=ba)
            for _tr in _st:
                ch = _tr.stats.channel
                _tr.stats.channel = '%s.%s' % (ch[-1], ch[:3])
                _tr.stats._component = ch[-1]
                _tr.stats.delta = self.dt
                _tr.stats.sampling_rate = int(1.0/self.dt)

            stream += _st
        return stream

