""" Main function of the ns2d solver """

import numpy as np
from map import *
from variable import *


#
#   User settings
#
setting = {
            'dim'       :   [4,5],
            'delta'     :   [1.0,1.0],
            'T'         :   [0.5,100],
            'pBC'       :   [1100.0,1000.0],
            'pIC'       :   1000.0,
            'rho'       :   1000.0,
            'nu'        :   0.0001
}

#
#   Main function
#

def main(setting):
    map = Map(setting)
    data = Variable(setting)
    data.enforceBoundary(map, setting)
    for tstep in range(setting['T'][1]):
        #   Execute one time step
        data.explicitSource(map, setting)
        data.matrixRHS(map, setting)
        data.buildMatrix(map, setting)
        data.solve()
        data.updateVelocity(map, setting)
        data.enforceBoundary(map, setting)
        print('Time step ',str(tstep),' executed!')
    return data
    
main(setting)
        
        
    
