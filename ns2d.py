""" Main function of the ns2d solver """

import numpy as np
import matplotlib.pyplot as plt
from map import *
from variable import *
from scipy.sparse import lil_matrix


#
#   User settings
#
plot = True
setting = {
            'dim'       :   [100,150],
            'delta'     :   [1.0,1.0],
            'T'         :   [1.0,30,20],
            'pBC'       :   10000.0,
            'vBC'       :   1e-1,
            'sBC'       :   100.0,
            'sIC'       :   0.0,
            'rho'       :   1000.0,
            'nu'        :   0.0001,
            'kappa'     :   0.01
}

#
#   Main function
#

def main(setting):
    map = Map(setting)
    data = Variable(map, setting)
    data.enforcePressureBC(map, setting)
    data.enforceVelocityBC(map, setting)
    for tstep in range(setting['T'][1]):
        #   Execute one time step
        data.explicitSource(map, setting)
        data.matrixRHS(map, setting)
        data.buildMatrix(map, setting)
        data.solve()
        data.enforcePressureBC(map, setting)
        data.updateVelocity(map, setting)
        data.enforceVelocityBC(map, setting)
        data.transport(map, setting)
        data.updateVariables(setting)
        print('Time step ',str(tstep),' executed!')
        # if tstep % setting['T'][2] == 0:
        #     makePlot(data, setting)
    return data

#
#   Plot function
#
def makePlot(data, setting):
    p = np.transpose(np.reshape(data.pp[:data.Ni], (setting['dim'][1],setting['dim'][0])))
    u = np.transpose(np.reshape(data.up[:data.Ni], (setting['dim'][1],setting['dim'][0])))
    v = np.transpose(np.reshape(data.vp[:data.Ni], (setting['dim'][1],setting['dim'][0])))
    s = np.transpose(np.reshape(data.sp[:data.Ni], (setting['dim'][1],setting['dim'][0])))
    u_max = np.amax(abs(u))
    v_max = np.amax(abs(v))
    fig, ax = plt.subplots()
    ax1 = plt.subplot(2,2,1)
    plt.imshow(p)
    cbar = plt.colorbar()
    ax1 = plt.subplot(2,2,2)
    plt.imshow(s, vmin=setting['sIC'], vmax=setting['sBC'])
    cbar = plt.colorbar()
    ax2 = plt.subplot(2,2,3)
    plt.imshow(u, vmin=-u_max, vmax=u_max)
    cbar = plt.colorbar()
    ax3 = plt.subplot(2,2,4)
    plt.imshow(v, vmin=-v_max, vmax=v_max)
    cbar = plt.colorbar()
    plt.show()

#
# Execution
#
if __name__=="__main__":
    data = main(setting)
    if plot == True:
        makePlot(data, setting)
