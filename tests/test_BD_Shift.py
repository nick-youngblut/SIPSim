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
from SIPSim.Commands import BD_Shift as BD_Shift_CMD
from SIPSim import Utils

# data dir
test_dir = os.path.join(os.path.dirname(__file__))
data_dir = os.path.join(test_dir, 'data')

# load kde
@pytest.fixture(scope="module")
def KDE1():
    kde = Utils.load_kde('ampFrag_skewN90-25-n5-nS.pkl')
    return kde

def KDE2():
    kde = Utils.load_kde('ampFrag_skewN90-25-n5-nS_kde.pkl')
    return kde

def KDE3():
    kde = Utils.load_kde('ampFrag_skewN90-25-n5-nS_dif_kde_DBL_incorp.pkl')
    return kde

def test_main():
    f1 = os.path.join(data_dir, 'ampFrag_skewN90-25-n5-nS_dif_kde.pkl')
    f2 = os.path.join(data_dir, 'ampFrag_skewN90-25-n5-nS_dif_kde_DBL.pkl')
    args = [f1, f2, '--debug']
    BD_Shift_CMD.opt_parse(args)
