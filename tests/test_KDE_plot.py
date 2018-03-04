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
from SIPSim.Commands import KDE_plot as KDE_plot_CMD

# data dir
test_dir = os.path.join(os.path.dirname(__file__))
data_dir = os.path.join(test_dir, 'data')


# tests
def test_cmd():
    kde_file = os.path.join(data_dir, 'ampFrag_skewN90-25-n5-nS_dif_kde.pkl')
    plot_file = os.path.join(data_dir, 'kde_plot.png')
    args = ['-o', plot_file, kde_file]
    KDE_plot_CMD.opt_parse(args)
