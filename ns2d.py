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
            'dim'       :   [20,30],
            'delta'     :   [1.0,1.0],
            'T'         :   [1.0,50,20],
            'pBC'       :   10000.0,
            'vBC'       :   1e-3,
            'sBC'       :   100.0,
            'sIC'       :   0.0,
            'rho'       :   1000.0,
            'nu'        :   0.00001,
            'kappa'     :   0.00001
}
block = np.zeros((setting['dim']), dtype=int)
inner = block[:,int(setting['dim'][0]/5):int(4*setting['dim'][0]/5)]
indim = inner.shape
for jj in range(indim[1]):
    inner[:min(int(indim[0]/3),int(jj*indim[0]/indim[1])),jj] = 1
    inner[max(int(2*indim[0]/3),int(jj*indim[0]/indim[1])):,jj] = 1
for jj in range(int(indim[1]/3)):
    inner[int(2*indim[0]/3):int(indim[0] - jj*indim[0]/indim[1]),jj] = 0
for jj in range(int(2*indim[1]/3),indim[1]):
    inner[int(indim[0] - jj*indim[0]/indim[1]):int(indim[0]/3),jj] = 0
block[:,int(setting['dim'][0]/5):int(4*setting['dim'][0]/5)] = inner

#
#   Main function
#

def main(block, setting):
    map = Map(block, setting)
    data = Variable(map, setting)
    data.enforcePressureBC(map, setting)
    data.enforceVelocityBC(map, setting)
    for tstep in range(setting['T'][1]):
        #   Execute one time step
        data.explicitSource(map, setting)
        data.matrixRHS(map, setting)
        data.buildMatrix(map, setting)
        data.solve()
        # Ad = lil_matrix.todense(data.A)
        # for ii in range(data.Ni):
        #     print(Ad[ii,:])
        # for ii in range(data.Ni):
        #     print([ii, data.Ex[ii], data.Ex[map.iMjc[ii]], data.Ey[ii], data.Ey[map.icjM[ii]], data.B[ii])
        data.enforcePressureBC(map, setting)
        data.updateVelocity(map, setting)
        data.enforceVelocityBC(map, setting)
        data.transport(map, setting)
        data.updateVariables(setting)
        print('Time step ',str(tstep),' executed!')
        # if tstep % setting['T'][2] == 0:
        #     makePlot(data, setting)
    return map, data

#
#   Plot function
#
def makePlot(data, map, setting):
    p = np.transpose(np.reshape(data.pp[:data.Ni], (setting['dim'][1],setting['dim'][0])))
    u = np.transpose(np.reshape(data.up[:data.Ni], (setting['dim'][1],setting['dim'][0])))
    v = np.transpose(np.reshape(data.vp[:data.Ni], (setting['dim'][1],setting['dim'][0])))
    s = np.transpose(np.reshape(data.sp[:data.Ni], (setting['dim'][1],setting['dim'][0])))
    block = np.transpose(np.reshape(map.block[:data.Ni], (setting['dim'][1],setting['dim'][0])))
    u_max = np.amax(abs(u))
    v_max = np.amax(abs(v))
    p_max = np.amax(abs(p))
    fig, ax = plt.subplots()
    ax1 = plt.subplot(2,2,1)
    plt.imshow(p, vmin=setting['pBC'], vmax=p_max)
    ax1.set_title('Pressure')
    cbar = plt.colorbar()
    ax2 = plt.subplot(2,2,2)
    plt.imshow(s, vmin=setting['sIC'], vmax=setting['sBC'])
    ax2.set_title('Tracer concentration')
    cbar = plt.colorbar()
    ax3 = plt.subplot(2,2,3)
    plt.imshow(u, vmin=-u_max, vmax=u_max)
    ax3.set_title('y-Velocity (positive downward)')
    cbar = plt.colorbar()
    ax4 = plt.subplot(2,2,4)
    plt.imshow(v, vmin=-v_max, vmax=v_max)
    ax4.set_title('x-Velocity (positive rightward)')
    cbar = plt.colorbar()
    plt.show()

#
# Execution
#
if __name__=="__main__":
    map, data = main(block, setting)
    if plot == True:
        makePlot(data, map, setting)
