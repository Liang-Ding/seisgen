# -------------------------------------------------------------------
# Widgets making life easier
#
# Author: Liang Ding
# Email: myliang.ding@mail.utoronto.ca
# -------------------------------------------------------------------

# DO NOT CHANGE
# if you utilize the pyCAPLunar package to create the SGT database.
SGT_ENCODING_LEVEL = 8
NGLL_SPEC_COMPRESSION = 27

NGLLX = 5
NGLLY = NGLLX
NGLLZ = NGLLX
CONSTANT_INDEX_27_GLL = [0, 2, 4, 10, 12, 14, 20, 22, 24,
                         50, 52, 54, 60, 62, 64, 70, 72, 74,
                         100, 102, 104, 110, 112, 114, 120, 122, 124]


def get_proc_name(idx_processor):
    '''Return the processor name.'''
    return str('proc')+str(idx_processor).rjust(6, '0')
