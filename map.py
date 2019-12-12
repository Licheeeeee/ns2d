""" A class that generates connection maps """
import numpy as np

class Map():
    """ Define cell connections """
    def __init__(self, block, setting):
        self.Nx = setting['dim'][0]
        self.Ny = setting['dim'][1]
        self.dx = setting['dim'][0]
        self.dy = setting['dim'][1]
        self.Ni = self.Nx * self.Ny
        self.Nt = (self.Nx+2) * (self.Ny+2)
        self.block = np.ravel(block,order='F')
        #   Build 2D (row,col) index of 1D cell
        self.row = []
        self.col = []
        for ii in range(self.Ny):
            for jj in range(self.Nx):
                self.row.append(jj)
                self.col.append(ii)
        #   Build center map
        self.cntr = []
        for ii in range(self.Ni):
            self.cntr.append(ii)
        #   Build iP map
        self.iPjc = []
        for ii in range(self.Ni):
            if (ii+1) % self.Nx == 0:
                self.iPjc.append(self.Ni + 2*self.Nx + self.col[ii])
            else:
                self.iPjc.append(ii+1)
        #    Build iM map
        self.iMjc = []
        for ii in range(self.Ni):
            if ii % self.Nx == 0:
                self.iMjc.append(self.Ni + 2*self.Nx + self.Ny + self.col[ii])
            else:
                self.iMjc.append(ii-1)
        #    Build jP map
        self.icjP = []
        for ii in range(self.Ni):
            if ii >= self.Nx*(self.Ny-1):
                self.icjP.append(self.Ni + self.row[ii])
            else:
                self.icjP.append(ii+self.Nx)
        #    Build jM map
        self.icjM = []
        for ii in range(self.Ni):
            if ii < self.Nx:
                self.icjM.append(self.Ni + self.Nx + self.row[ii])
            else:
                self.icjM.append(ii-self.Nx)
