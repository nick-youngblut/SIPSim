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
from SIPSim.Commands import OTU_sampleData as OTU_sampleData_CMD

# data dir
test_dir = os.path.join(os.path.dirname(__file__))
data_dir = os.path.join(test_dir, 'data')


# tests
def test_cmd():
    otu_file = os.path.join(data_dir, 'ampFrag_OTU_n2_abs1e9.txt')
    args = [otu_file]
    OTU_sampleData_CMD.opt_parse(args)
