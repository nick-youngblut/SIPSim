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
from SIPSim.Commands import Fragment_KDE_cat as Fragment_KDE_cat_CMD
from SIPSim import Utils

# data dir
test_dir = os.path.join(os.path.dirname(__file__))
data_dir = os.path.join(test_dir, 'data')

# load kde
def test_cmd():
    f1 = os.path.join(data_dir, 'ampFrag_skewN90-25-n5-nS_dif_kde.pkl')
    f2 = os.path.join(data_dir, 'ampFrag_skewN90-25-n5-nS_dif_kde_DBL.pkl')
    args = [f1, f2]
    Fragment_KDE_cat_CMD.opt_parse(args)
