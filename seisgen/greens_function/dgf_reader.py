# -------------------------------------------------------------------
# Tools to acquire Greens Function (the displacement, DGF)
# from the 3D pre-stored receiver-side Greens function database
#
# Author: Liang Ding
# Email: myliang.ding@mail.utoronto.ca
# -------------------------------------------------------------------


import numpy as np
import h5py
import zlib
from seisgen.util_SPECFEM3D import SGT_ENCODING_LEVEL


DGF_KEYS = [
    'index',
    'start',
    'length',
    'offset',
    'scale',
]

DGF_ATTRS = [
    'ngll',
    'nstep',
    'nforce',
    'nparas',
    'dt',
    'nspec',
    'nGLL_global',
    'type',
]


def read_DGF_header(header_path):
    '''Acquire parameters in the header file. '''
    # load the information of about the GF database.
    with h5py.File(header_path, 'r') as f:
        dt      = f.attrs[DGF_ATTRS[4]]
        NSPEC   = f.attrs[DGF_ATTRS[5]]
        nGLL_global = f.attrs[DGF_ATTRS[6]]
    return dt, NSPEC, nGLL_global


def DEnquire_DGF(data_path, header_path, g_indx_GLL_points):
    '''
    * Enquire the GF from the database (*.bin files).

    :param data_path:           The path to the data file (.bin)
    :param header_path:         The path to the header file (.hdf5)
    :param g_indx_GLL_points:   The global index of the GLL points acquired in the slice.
    :return:                    The list of GF array acquired.
    '''

    with h5py.File(header_path, 'r') as f:
        names_GLL_arr       = f[DGF_KEYS[0]][:]
        data_start_array    = f[DGF_KEYS[1]][:]
        data_length_array   = f[DGF_KEYS[2]][:]
        data_offset_array   = f[DGF_KEYS[3]][:]
        data_scale_array    = f[DGF_KEYS[4]][:]

        n_gll   = f.attrs[DGF_ATTRS[0]]
        n_step  = f.attrs[DGF_ATTRS[1]]
        n_dim   = f.attrs[DGF_ATTRS[2]]
        n_paras = f.attrs[DGF_ATTRS[3]]
        dt      = f.attrs[DGF_ATTRS[4]]
        NSPEC   = f.attrs[DGF_ATTRS[5]]
        nGLL_global = f.attrs[DGF_ATTRS[6]]

    dgf_arr_list = []
    with open(data_path, "rb") as fr:
        for gll in g_indx_GLL_points:
            idx_gll = (np.where(names_GLL_arr == gll))[0][0]

            offset_min = data_offset_array[idx_gll]
            normal_factor = data_scale_array[idx_gll]
            sgt_begin_bytes = data_start_array[idx_gll]
            sgt_length_bytes = data_length_array[idx_gll]

            # extract the compressed data.
            fr.seek(sgt_begin_bytes)
            data = fr.read(sgt_length_bytes)

            # uncompress the data to bytes.
            data = zlib.decompress(data)

            # recover bytes to uint8
            data = np.frombuffer(data, dtype=np.uint8)

            # recover waveform
            data = data / (2 ** SGT_ENCODING_LEVEL - 1) * normal_factor + offset_min

            # data container
            unzip_sgt = np.zeros([n_step, n_dim, n_paras]).astype(np.float32)
            count = 0
            for j in range(n_dim):
                for k in range(n_paras):
                    unzip_sgt[:, k, j] = data[count * n_step:(count + 1) * n_step]
                    count += 1

            dgf_arr_list.append(unzip_sgt)

    return dgf_arr_list