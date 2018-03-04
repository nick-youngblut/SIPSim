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
from SIPSim.Commands import Diffusion as Diffusion_CMD

# data dir
test_dir = os.path.join(os.path.dirname(__file__))
data_dir = os.path.join(test_dir, 'data')


# tests
def test_main():
    kde_file = os.path.join(data_dir, 'ampFrag_skewN90-25-n5-nS_kde.pkl')
    diff_idx_file = os.path.join(data_dir, 'diffusion_index.txt')
    args = ['--index_out', diff_idx_file, kde_file]
    Diffusion_CMD.opt_parse(args)
