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
from SIPSim.Commands import Isotope_incorp as Isotope_incorp_CMD

# data dir
test_dir = os.path.join(os.path.dirname(__file__))
data_dir = os.path.join(test_dir, 'data')


# tests
def test_cmd():
    kde_file = os.path.join(data_dir, 'ampFrag_skewN90-25-n5-nS_dif_kde.pkl')
    config_file = os.path.join(data_dir, 'incorp50.config')
    shift_file = os.path.join(data_dir, 'BD_shift_stats.txt')
    args = ['--shift', shift_file, kde_file, config_file]
    Isotope_incorp_CMD.opt_parse(args)
