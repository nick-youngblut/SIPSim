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
from SIPSim.Commands import OTU_table as OTU_table_CMD

# data dir
test_dir = os.path.join(os.path.dirname(__file__))
data_dir = os.path.join(test_dir, 'data')


# tests
def test_cmd():
    kde_file = os.path.join(data_dir, 'ampFrag_skewN90-25-n5-nS_dif_kde.pkl')
    comm_file = os.path.join(data_dir, 'comm-n2-unif.txt')
    frac_file = os.path.join(data_dir, 'fracs-n2-unif.txt')
    args = [kde_file, comm_file, frac_file]
    OTU_table_CMD.opt_parse(args)
