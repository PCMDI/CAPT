#! /usr/bin/env python

#==============================================================================
#   filter_field.py
"""
Python module.

This module is to applied numerical filtering/smoothing on the reanalysis fields.

There are two options:

    1. bleck:

       Bleck, R., 1974: Short-range prediction in isentropic coordinates with
           filterred and unfilterred numerical models. Mon. Wea. Rev., 102, 
           813-829.

    2. gerrity:

       Gerrity, J. P. and R. D. McPherson, 1970: Noise analysis of a limited-area
           fine-mesh prediction model. ESSA Technical Memoranda, WBTM NMC 46, 
           PB-191-188. 81pp. 

Adapted for numpy/ma/cdms2 by convertcdms.py.
"""
#==============================================================================

import numpy
n = numpy

import filterG
import filter25


def filter_field(aa, filter_type, filter_index, xwrap=True):
#-----------------------------------------------------------
    """
    See file header.
    """

    ret_val = ''

    if xwrap:
        wrap_col0 = aa[:,-5:]
        wrap_colx = aa[:,0:5]

        aa_wrap = n.concatenate((wrap_col0, aa, wrap_colx), 1)

    else:
        aa_wrap = aa

    work = aa_wrap * 1.0

    if filter_type == 'bleck':
        aa_wrap = \
            n.transpose(filter25.filter25(n.transpose(aa_wrap), \
                                          n.transpose(work), filter_index))

        ret_val = aa_wrap[:,5:-5]

    elif filter_type == 'gerrity':
        aa_wrap = \
            n.transpose(filterG.filterg(n.transpose(aa_wrap), \
                                        n.transpose(work), filter_index))

        ret_val = aa_wrap[:,5:-5]

    else:
        print 'FILTER_TYPE UNKNOWN!: ', filter_type

        ret_val = None

    return ret_val
