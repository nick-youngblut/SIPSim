#!/usr/bin/env python
# import
from __future__ import print_function
## batteries
import os
import sys
import pytest
## 3rd party
import numpy as np
import pandas as pd
## package
from SIPSim import Utils

# data dir
test_dir = os.path.join(os.path.dirname(__file__))
data_dir = os.path.join(test_dir, 'data')


# tests
def test_file_status():
    # real file
    f = os.path.join(data_dir, 'ampFrag_skewN90-25-n5-nS_kde.pkl')
    Utils.is_file(f)
    Utils.checkExists(f)
    Utils.checkEmpty(f)
    # fake file
    f = os.path.join(data_dir, 'DOES_NOT_EXIST')
    with pytest.raises(IOError):
        Utils.is_file(f)
    with pytest.raises(IOError):
        Utils.checkExists(f)
    
def test_kde():
    f = os.path.join(data_dir, 'ampFrag_skewN90-25-n5-nS_kde.pkl')
    # load kde
    kde = Utils.load_kde(f)
    assert isinstance(kde, dict)
    # check type
    kde_type = Utils.KDE_type(kde)
    assert kde_type == 2
    
# def test_exp_design():
#     f = os.path.join(data_dir, 'comm-n2-unif.txt')
#     x = Utils.load_exp_design(f)
#     assert isinstance(x, pd.dataframe)
