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
from seisgen.MTTools.DMomentTensors import DMT_enz
from seisgen.util_SPECFEM3D.ibool_reader import DEnquire_Element
from seisgen.util_SPECFEM3D.xyz_reader import DEnquire_XYZ_GLLs_Element
from seisgen.greens_function.sgt_reader import DEnquire_SGT, read_header_info
from seisgen.math.interp_tools import DCreate_anchors_xi_eta_gamma, DLagrange_interp_sgt, DLagrange_any3D

from obspy.core.util.attribdict import AttribDict
from obspy.core import Stream, Trace
from obspy.clients.iris import Client
from obspy.taup import TauPyModel

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

# fundamental faults
FF_ELEMENT = [
            "MEP",
            "MDD",
            "MDS",
            "MSS"
        ]


class DSGTMgr(DPointCloud):
    '''Strain Green's Tensor (SGT) database Manager'''

    def __init__(self, sgt_database_folder, model3D_folder, point_cloud_file, DLite=False):
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
        if DLite:
            pass
        else:
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
                                     str(self.station.network),
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

    def set_dt(self, dt):
        '''
        Set time interval.
        Function observed for SeisClient.
        '''
        self.dt = dt
        return self


    def get_sgt(self, station, origin, b_new_origin=True, b_verbose=False):
        '''
        Get the interpolated SGT between the station-origin pair.
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
        :param b_new_origin: True, If acquiring SGF at multiple stations for a same origin to save time.
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

        :param b_new_origin: Accelerating the extraction of SGT data for multiple stations with the same origin.
        '''

        sgt = self.get_sgt(station, origin, b_new_origin=b_new_origin)
        return self.get_greens_function_next(sgt, station, origin, b_new_origin)


    def get_greens_function_next(self, sgt, station, origin):
        '''The next step of get_greens_function()'''
        client = Client()
        res = client.distaz(stalat=station.latitude, stalon=station.longitude,
                            evtlat=origin.latitude, evtlon=origin.longitude)
        back_azimuth = res['backazimuth']
        azimuth = res['azimuth']
        distance_deg = res['distance']
        distance_m = res['distancemeters']
        stream = self._SGT2GF(sgt, azimuth, back_azimuth)
        stream.id = station.id
        return stream


    def get_fk_greens_function(self, station, origin, b_new_origin=True, b_save=False, greens_path=None):
        '''
        Get FK-type Greens Function between the station-origin pair.
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

        :param b_new_origin: Accelerating the extraction of SGT data for multiple stations with the same origin.
        :param b_save:  Whether to save the greens function to file?
        :para greens_path: the path to store the Greens function
        '''
        sgt = self.get_sgt(station, origin, b_new_origin=b_new_origin)
        return self.get_fk_greens_function_next(sgt, station, origin, b_save, greens_path)


    def get_fk_greens_function_next(self, sgt, station, origin, b_save=False, greens_path=None):
        '''The next step of get_fk_greens_function()'''

        client = Client()
        res = client.distaz(stalat=station.latitude, stalon=station.longitude,
                            evtlat=origin.latitude, evtlon=origin.longitude)

        azimuth = res['azimuth']
        back_azimuth = res['backazimuth']
        distance_deg = res['distance']
        distance_m = res['distancemeters']
        stream = self._SGT2FKGF(sgt, azimuth, back_azimuth)
        stream.id = station.id

        # Add SAC header
        b_sacHdeader = True
        try:
            depth_src_km = origin.depth_in_m / 1000.0
        except:
            b_sacHdeader = False

        if b_sacHdeader:
            model = TauPyModel(model="ak135")
            try:
                arrivals = model.get_travel_times(source_depth_in_km=depth_src_km,
                                                  distance_in_degree=distance_deg,
                                                  phase_list=["p", "s"])
                sac = AttribDict()
                sac.o = 0.
                sac.dist, sac.az, sac.baz = distance_m / 1000, azimuth, back_azimuth
                sac.stla = station.latitude
                sac.stlo = station.longitude
                sac.evla = origin.latitude
                sac.evlo = origin.longitude
                sac.evdp = depth_src_km  # in km
                sac.b = 0.0
                try:
                    sac.t1 = arrivals[0].time
                except:
                    pass
                try:
                    sac.t2 = arrivals[1].time
                except:
                    pass

                for tr in stream.traces:
                    tr.stats.sac = sac
            except:
                pass

        if b_save:
            if not os.path.exists(greens_path):
                os.makedirs(greens_path)

        # save FK-type Greens function to files.
        if b_save and greens_path is not None:
            chs = ['ZDD', 'RDD', 'TDD',
                   'ZDS', 'RDS', 'TDS',
                   'ZSS', 'RSS', 'TSS',
                   'ZEP', 'REP', 'TEP']

            fk_chs = ['0', '1', '2',
                      '3', '4', '5',
                      '6', '7', '8',
                      'a', 'b', 'c']

            for i, ch in enumerate(chs):
                file_path = os.path.join(greens_path, "%d.grn.%s" % (np.ceil(distance_m/1000.0), fk_chs[i]))
                stream.select(channel="%s" % ch).write(file_path, format='SAC')

        return stream


    def get_waveform(self, station, origin, mt_RTP, b_RTZ=False, b_new_origin=True):
        '''
        Return the synthetic 3-C waveform in RTZ.
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

        :param mt: the moment tensor in RTP [Mrr, Mtt, Mpp, Mrt, Mrp, Mtp]
        :param b_new_origin: Accelerating the extraction of SGT data for multiple stations with the same origin.
        '''

        sgt = self.get_sgt(station, origin, b_new_origin=b_new_origin)
        return self.get_waveform_next(sgt, station, origin, mt_RTP, b_RTZ)


    def get_waveform_next(self, sgt, station, origin, mt_RTP, b_RTZ=False):
        '''The next step of get_waveform()'''

        client = Client()
        res = client.distaz(stalat=station.latitude, stalon=station.longitude,
                            evtlat=origin.latitude, evtlon=origin.longitude)
        back_azimuth = res['backazimuth']
        distance_deg = res['distance']
        distance_m = res['distancemeters']

        # synthetic
        element='SY'
        mt_enz = RTP_to_DENZ(mt_RTP)
        _st = DSyn(mt_enz, sgt, element)
        for _tr in _st:
            _tr.stats.delta = self.dt
            _tr.stats.sampling_rate = int(1.0 / self.dt)

        # waveform in ENZ or RTZ convention
        if b_RTZ:
            _st.rotate(method='NE->RT', back_azimuth=back_azimuth)
        return _st


    def _SGT2GF(self, sgt, azi, ba, b_USE=True):
        '''
        Get 3D MT Greens functions
        ( Up-South-East convention (USE) by default, which is compatible with MTUQ,
        otherwise North-East-Down convention)
        '''

        fk_grn_st = self._SGT2FKGF(sgt, azi, ba)
        stream = Stream()

        az = np.deg2rad(azi)
        sa = np.sin(az)
        ca = np.cos(az)
        sa2 = np.sin(2.0*az)
        ca2 = np.cos(2.0*az)
        fk0 = fk_grn_st.select(channel='ZDD').traces[0].data
        fk1 = fk_grn_st.select(channel='RDD').traces[0].data
        fk2 = fk_grn_st.select(channel='TDD').traces[0].data
        fk3 = fk_grn_st.select(channel='ZDS').traces[0].data
        fk4 = fk_grn_st.select(channel='RDS').traces[0].data
        fk5 = fk_grn_st.select(channel='TDS').traces[0].data

        fk6 = fk_grn_st.select(channel='ZSS').traces[0].data
        fk7 = fk_grn_st.select(channel='RSS').traces[0].data
        fk8 = fk_grn_st.select(channel='TSS').traces[0].data
        fk9 = fk_grn_st.select(channel='ZEP').traces[0].data
        fk10 = fk_grn_st.select(channel='REP').traces[0].data
        fk11 = fk_grn_st.select(channel='TEP').traces[0].data

        # [a-r] = [0-17]
        n_component = 18
        npts = len(fk0)
        mt = np.zeros([n_component, npts])
        mt[0] = fk9 / 3. - fk0 / 6. - fk6 * ca2 / 2.
        mt[1] = fk10 / 3. - fk1 / 6. - fk7 * ca2 / 2.
        mt[2] = -1. * fk8 * sa2 / 2.
        mt[3] = fk6 * sa2
        mt[4] = fk7 * sa2
        mt[5] = -1. * fk8 * ca2
        mt[6] = -1. * fk3 * ca
        mt[7] = -1. * fk4 * ca
        mt[8] = -1. * fk5 * sa
        mt[9] = fk9 / 3. - fk0 / 6. + fk6 * ca2 / 2.
        mt[10] = fk10 / 3. - fk1 / 6. + fk7 * ca2 / 2.
        mt[11] = -1. * mt[2]
        mt[12] = fk3 * sa
        mt[13] = fk4 * sa
        mt[14] = -1. * fk5 * ca
        mt[15] = (fk9 + fk0) / 3.0
        mt[16] = (fk10 + fk1) / 3.0
        mt[17] = np.zeros(npts)
        grn_chs = ['Z.Mtt', 'R.Mtt', 'T.Mtt',
                   'Z.Mtp', 'R.Mtp', 'T.Mtp',
                   'Z.Mrt', 'R.Mrt', 'T.Mrt',
                   'Z.Mpp', 'R.Mpp', 'T.Mpp',
                   'Z.Mrp', 'R.Mrp', 'T.Mrp',
                   'Z.Mrr', 'R.Mrr', 'T.Mrr']

        scales = np.array([1, 1, 1,
                           -1, -1, -1,
                           1, 1, 1,
                           1, 1, 1,
                           -1, -1, -1,
                           1, 1, 0])

        if b_USE:
            # Up - South - East convention
            for i in range(n_component):
                _tr = Trace(mt[i])
                _tr.stats.channel = grn_chs[i]
                _tr.stats.delta = self.dt
                _tr.stats.sampling_rate = int(1.0/self.dt)
                stream.append(_tr)
        else:
            # North - East - Down convention
            for i in range(n_component):
                _tr = Trace(mt[i] * scales[i])
                _tr.stats.channel = grn_chs[i]
                _tr.stats.delta = self.dt
                _tr.stats.sampling_rate = int(1.0/self.dt)
                stream.append(_tr)
        return stream


    def _SGT2FKGF(self, sgt, azi, ba):
        '''Generate fk-type Greens functions. [*.0-8] from DC, [*.a-c] from EP.'''
        # Fundamental faults:
        # EP: Explosion source
        # DD: 45-degree-dip slip
        # DS: Vertical dip slip
        # SS: vertical strike slip

        # calculate the moment tensor (ENZ) of fundamental faults (EXP, DD, DS, SS)
        mt_enz_ff = []
        strike = [np.mod(azi + 270, 360)-45, np.mod(azi + 270, 360)-45, np.mod(azi + 270, 360)-67.5]
        dip_arr = np.array([45, 90, 90])
        rake_arr = np.array([90, 90, 0])
        colatitude, lune_longitude = 90, 0
        _mt_EXP = np.array([1, 1, 1, 0, 0, 0])  # EP
        mt_enz_ff.append(_mt_EXP)

        # DD, DS, SS
        for i in range(3):
            _mt_enz = DMT_enz(np.deg2rad(strike[i]), np.deg2rad(dip_arr[i]), np.deg2rad(rake_arr[i]),
                              np.deg2rad(colatitude), np.deg2rad(lune_longitude))
            _mt_enz[3:] *= 2
            mt_enz_ff.append(_mt_enz)
        mt_enz_ff = np.asarray(mt_enz_ff)
        sqrt2 = np.sqrt(2)
        scaling = np.array([
            [1, 1, 0],
            [2, 2, 0],
            [sqrt2, sqrt2, sqrt2],
            [sqrt2, sqrt2, sqrt2],
        ])

        stream = Stream()
        for i, mt_enz in enumerate(mt_enz_ff):
            _st = DSyn(mt_enz, sgt, FF_ELEMENT[i])
            _st.rotate(method='NE->RT', back_azimuth=ba)
            for _tr in _st:
                ch = _tr.stats.channel
                if ch[-1] == 'Z':
                    _tr.data *= scaling[i][0]
                elif ch[-1] == 'R':
                    _tr.data *= scaling[i][1]
                elif ch[-1] == 'T':
                    _tr.data *= scaling[i][2]

                _tr.stats.channel = '%s%s' % (ch[-1], ch[1:3])
                _tr.stats._component = ch[-1]
                _tr.stats.delta = self.dt
                _tr.stats.sampling_rate = int(1.0 / self.dt)

            stream += _st
        return stream
