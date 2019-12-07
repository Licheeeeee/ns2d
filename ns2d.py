""" Main function of the ns2d solver """

import numpy as np
from map import *


#
#   User settings
#
setting = {
            'dim'       :   [4,5],
            'delta'     :   [1.0,1.0],
            'T'         :   [0.5,100],
            'pBC'       :   [1100.0,1000.0]
}

map = Map(setting)
