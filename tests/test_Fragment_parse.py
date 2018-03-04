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
from SIPSim.Commands import Fragment_parse as Fragment_parse_CMD
from SIPSim import Utils

# data dir
test_dir = os.path.join(os.path.dirname(__file__))
data_dir = os.path.join(test_dir, 'data')

# load kde
def test_cmd():
    frag_file = os.path.join(data_dir, 'ampFrag_skewN90-25-n5-nS.pkl')
    genome_file = os.path.join(data_dir, 'genome_parse.txt')
    log_file = os.path.join(data_dir, 'genome_parse.log')
    args = ['--log', log_file, frag_file, genome_file]
    Fragment_parse_CMD.opt_parse(args)
