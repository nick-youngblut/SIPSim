#!/usr/bin/env python
# import
from __future__ import print_function
## batteries
import os
import sys
import configobj
## 3rd party
import pytest
import numpy as np
import pandas as pd
## package
from SIPSim import Config

# data dir
test_dir = os.path.join(os.path.dirname(__file__))
data_dir = os.path.join(test_dir, 'data')


# tests
def test_example_config():
    basicConfig = Config.get_basicConfig()
    cfg = Config.ExampleConfig(basicConfig)
    assert isinstance(cfg, Config.ExampleConfig)

def test_load_config():
    f = os.path.join(data_dir, 'incorp.config')
    config_spec = Config.get_configspec()
    cfg = Config.Config(configspec=config_spec)
    config = cfg.load_config(f)
    assert isinstance(config, configobj.ConfigObj)


    
