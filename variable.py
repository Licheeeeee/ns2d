""" Create the data structure that stores the variables """
import numpy as np
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import spsolve

class Variable():
    def __init__(self, setting):
        self.Nx = setting['dim'][0]
        self.Ny = setting['dim'][1]
        self.dx = setting['dim'][0]
        self.dy = setting['dim'][1]
        self.Ni = self.Nx * self.Ny
        self.Nt = (self.Nx+2) * (self.Ny+2)
        self.dt = setting['T'][0]
        self.rho = setting['rho']
        self.nu = setting['nu']
        self.pBC = setting['pBC']
        self.pIC = setting['pIC']
        #   Variables that remain constants for now
        self.V = self.dx * self.dy
        self.Ax = self.dy
        self.Ay = self.dx
        #   Variables on cell centers
        self.pn = self.pIC * np.ones(self.Nt)
        self.pp = self.pIC * np.ones(self.Nt)
        self.B = np.zeros((self.Ni,1))
        #   Variables on cell faces
        self.un = np.zeros(self.Nt)
        self.up = np.zeros(self.Nt)
        self.vn = np.zeros(self.Nt)
        self.vp = np.zeros(self.Nt)
        self.Gxp = (self.dt * self.Ax / (self.rho * self.V)) * np.ones(self.Ni)
        self.Gxm = (self.dt * self.Ax / (self.rho * self.V)) * np.ones(self.Ni)
        self.Gyp = (self.dt * self.Ax / (self.rho * self.V)) * np.ones(self.Ni)
        self.Gym = (self.dt * self.Ax / (self.rho * self.V)) * np.ones(self.Ni)
        self.Gct = self.Gxp + self.Gxm + self.Gyp + self.Gym
        #   Matrix coefficients
        self.Gx = self.dt * self.Ax / (self.rho * self.V)
        self.Gy = self.dt * self.Ay / (self.rho * self.V)
        self.Gc = 2.0 * self.Gx + 2.0 * self.Gy
        
    def explicitSource(self, map, setting):
        """ Calculate the explicit terms """
        self.Ex = np.zeros(self.Ni)
        self.Ey = np.zeros(self.Ni)
        for ii in range(self.Ni):
            #   Explicit velocity
            self.Ex[ii] = self.un[ii]
            self.Ey[ii] = self.vn[ii]
            #   Advection terms
            self.Ex[ii] -= (0.5 * self.dt / self.dx) * ((self.un[ii] + abs(self.un[ii])) * (self.un[ii] - self.un[map.iMjc[ii]]) + (self.un[ii] - abs(self.un[ii])) * (self.un[ii] - self.un[map.iPjc[ii]]))
            self.Ex[ii] -= (0.5 * self.dt / self.dy) * ((self.vn[ii] + abs(self.vn[ii])) * (self.un[ii] - self.un[map.icjM[ii]]) + (self.vn[ii] - abs(self.vn[ii])) * (self.un[ii] - self.un[map.icjP[ii]]))
            self.Ey[ii] -= (0.5 * self.dt / self.dx) * ((self.un[ii] + abs(self.un[ii])) * (self.vn[ii] - self.vn[map.iMjc[ii]]) + (self.un[ii] - abs(self.un[ii])) * (self.vn[ii] - self.vn[map.iPjc[ii]]))
            self.Ey[ii] -= (0.5 * self.dt / self.dy) * ((self.vn[ii] + abs(self.vn[ii])) * (self.vn[ii] - self.vn[map.icjM[ii]]) + (self.vn[ii] - abs(self.vn[ii])) * (self.vn[ii] - self.vn[map.icjP[ii]]))
            #   Diffusion terms
            self.Ex[ii] += (self.dt * self.nu * self.Ax / (self.dx * self.V)) * (self.un[map.iPjc[ii]] - 2.0*self.un[ii] + self.un[map.iMjc[ii]])
            self.Ex[ii] += (self.dt * self.nu * self.Ay / (self.dy * self.V)) * (self.un[map.icjP[ii]] - 2.0*self.un[ii] + self.un[map.icjM[ii]])
            self.Ey[ii] += (self.dt * self.nu * self.Ax / (self.dx * self.V)) * (self.vn[map.iPjc[ii]] - 2.0*self.vn[ii] + self.vn[map.iMjc[ii]])
            self.Ey[ii] += (self.dt * self.nu * self.Ay / (self.dy * self.V)) * (self.vn[map.icjP[ii]] - 2.0*self.vn[ii] + self.vn[map.icjM[ii]])
    
    def matrixRHS(self, map, setting):
        """ Calculate right-hand-side of the linear system """
        self.B = np.zeros(self.Ni)
        for ii in range(self.Ni):
            self.B[ii] = -(self.Ex[ii] + self.Ex[map.iMjc[ii]] + self.Ey[ii] + self.Ey[map.icjM[ii]])
            
    def updateVelocity(self, map, setting):
        """ Update velocity based on pressure """
        for ii in range(self.Ni):
            self.up[ii] = -self.Gx * (self.pp[map.iPjc[ii]] - self.pp[ii]) + self.Ex[ii]
            self.vp[ii] = -self.Gy * (self.pp[map.icjP[ii]] - self.pp[ii]) + self.Ey[ii]
            
    def enforceBoundary(self, map, setting):
        """ Enforce BC on pressure and velocity """
        for ii in range(self.Nx):
            jj = self.Nx * (self.Ny-1) + ii
            self.up[map.icjM[ii]] = self.up[ii]
            self.up[map.icjP[jj]] = self.up[jj]
            self.vp[map.icjM[ii]] = self.vp[ii]
            self.vp[map.icjP[jj]] = self.vp[jj]
            self.pp[map.icjM[ii]] = self.pp[ii]
            self.pp[map.icjP[jj]] = self.pp[jj]
        for ii in range(self.Ny):
            jj = ii * self.Nx
            kk = (ii+1) * self.Nx - 1
            self.up[map.iMjc[jj]] = self.up[jj]
            self.up[map.iPjc[kk]] = self.up[kk]
            self.vp[map.iMjc[jj]] = self.vp[jj]
            self.vp[map.iPjc[kk]] = self.vp[kk]
            self.pp[map.iMjc[jj]] = self.pp[jj]
            self.pp[map.iPjc[kk]] = self.pp[kk]
            
    def buildMatrix(self, map, setting):
        """ Build sparse matrix A """
        self.A = csc_matrix((self.Ni, self.Ni), dtype=float)
        for ii in range(self.Ni):
            self.A[ii,ii] = self.Gc
            #   For Neumann boundary
            if map.iMjc[ii] >= self.Ni:
                self.A[ii,ii] -= self.Gx
            else:
                self.A[ii,ii-1] = -self.Gx
            if map.iPjc[ii] >= self.Ni:
                self.A[ii,ii] -= self.Gx
            else:
                self.A[ii,ii+1] = -self.Gx
            #   For Direchlet boundary
            if map.icjM[ii] >= self.Ni:
                self.A[ii,ii] = 1.0
                self.B[ii] = self.pBC[0]
            else:
                self.A[ii,ii-self.Nx] = -self.Gy
            if map.icjP[ii] >= self.Ni:
                self.A[ii,ii] = 1.0
                self.B[ii] = self.pBC[1]
            else:
                self.A[ii,ii+self.Nx] = -self.Gy
            
    def solve(self):
        """ Solve the linear system Ax = B """
        self.x = spsolve(self.A, self.B)
        for ii in range(self.Ni):
            self.un[ii] = self.up[ii]
            self.vn[ii] = self.vp[ii]
            self.pn[ii] = self.pp[ii]
            self.pp[ii] = self.x[ii]
        
