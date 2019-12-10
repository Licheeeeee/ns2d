""" Create the data structure that stores the variables """
import numpy as np
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve

class Variable():
    def __init__(self, map, setting):
        self.Nx = setting['dim'][0]
        self.Ny = setting['dim'][1]
        self.dx = setting['delta'][0]
        self.dy = setting['delta'][1]
        self.Ni = self.Nx * self.Ny
        self.Nt = (self.Nx+2) * (self.Ny+2)
        self.dt = setting['T'][0]
        self.rho = setting['rho']
        self.nu = setting['nu']
        self.pBC = setting['pBC']
        self.pIC = setting['pIC']
        self.sBC = setting['sBC']
        self.sIC = setting['sIC']
        #   Variables that remain constants for now
        self.V = self.dx * self.dy
        self.Ax = self.dy
        self.Ay = self.dx
        #   Variables on cell centers
        self.pn = self.pIC * np.ones(self.Nt)
        self.pp = self.pIC * np.ones(self.Nt)
        for ii in range(self.Nx):
            jj = self.Nx * (self.Ny-1) + ii
            self.pn[ii] = self.pBC[0]
            self.pp[ii] = self.pBC[0]
            self.pn[jj] = self.pBC[1]
            self.pp[jj] = self.pBC[1]
        # use linear initial Pressure
        for ii in range(self.Ni):
            col = map.col[ii]
            self.pn[ii] = self.pBC[0] - col / (self.pBC[0]-self.pBC[1])
        self.sn = self.sIC * np.ones(self.Nt)
        self.sp = self.sIC * np.ones(self.Nt)
        self.B = np.zeros((self.Ni,1))
        #   Variables on cell faces
        self.un = np.zeros(self.Nt)
        self.up = np.zeros(self.Nt)
        self.vn = np.zeros(self.Nt)
        self.vp = np.zeros(self.Nt)
        #   Matrix coefficients
        self.Gx = self.dt * self.Ax / (self.rho * self.V)
        self.Gy = self.dt * self.Ay / (self.rho * self.V)
        self.Gc = 2.0 * self.Gx + 2.0 * self.Gy

    def explicitSource(self, map, setting):
        """ Calculate the explicit terms """
        self.Ex = np.zeros(self.Nt)
        self.Ey = np.zeros(self.Nt)
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
        for ii in range(self.Nx):
            self.Ey[map.icjM[ii]] = self.Ey[ii]
        for ii in range(self.Ny):
            jj = ii * self.Nx
            self.Ex[map.iMjc[jj]] = self.Ex[jj]

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

    def enforcePressureBC(self, map, setting):
        """ Enforce BC on pressure """
        for ii in range(self.Nx):
            jj = self.Nx * (self.Ny-1) + ii
            self.pp[map.icjM[ii]] = self.pp[ii]
            # self.pp[map.icjP[jj]] = 2.0 * self.pp[jj] - self.pp[map.icjM[jj]]
            self.pp[map.icjP[jj]] = self.pp[jj]
        for ii in range(self.Ny):
            jj = ii * self.Nx
            kk = (ii+1) * self.Nx - 1
            self.pp[map.iMjc[jj]] = self.pp[jj]
            self.pp[map.iPjc[kk]] = self.pp[kk]

    def enforceVelocityBC(self, map, setting):
        """ Enforce BC on velocity """
        #   BC on ghost cells
        for ii in range(self.Nx):
            jj = self.Nx * (self.Ny-1) + ii
            self.up[map.icjM[ii]] = self.up[ii]
            self.up[map.icjP[jj]] = self.up[jj]
            self.vp[jj] = self.vp[map.icjM[jj]] + (self.up[map.iMjc[jj]]-self.up[jj]) * self.Ax / self.Ay
            self.vp[map.icjM[ii]] = self.vp[ii] + (self.up[ii]-self.up[map.iMjc[ii]]) * self.Ax / self.Ay
            self.vp[map.icjP[jj]] = self.vp[jj]
            # self.vp[map.icjP[jj]] = self.vp[jj] + (self.up[map.iMjc[jj]]-self.up[jj]) * self.Ax / self.Ay
        for ii in range(self.Ny):
            jj = ii * self.Nx
            kk = (ii+1) * self.Nx - 1
            self.up[kk] = 0.0
            self.up[map.iMjc[jj]] = 0.0
            self.up[map.iPjc[kk]] = self.up[kk]
            self.vp[map.iMjc[jj]] = self.vp[jj]
            self.vp[map.iPjc[kk]] = self.vp[kk]

    def buildMatrix(self, map, setting):
        """ Build sparse matrix A """
        self.A = lil_matrix((self.Ni, self.Ni), dtype=float)
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
                self.A[ii,ii] -= self.Gy
            else:
                self.A[ii,ii-self.Nx] = -self.Gy
            if map.icjP[ii] >= self.Ni:
                self.A[ii,ii] -= self.Gy
                # self.B[ii] += self.Gy * self.pBC[1]
            else:
                self.A[ii,ii+self.Nx] = -self.Gy
        #   Enforce pressure BC
        for ii in range(self.Ni):
            if map.icjM[ii] >= self.Ni:
                self.A[ii,ii] = 1.0
                self.A[ii,ii+1] = 0.0
                self.A[ii,ii+self.Nx] = 0.0
                self.B[ii] = self.pBC[0]
                if ii > 0:
                    self.A[ii,ii-1] = 0.0
            if map.icjP[ii] >= self.Ni:
                self.A[ii,ii] = 1.0
                self.A[ii,ii-1] = 0.0
                self.A[ii,ii-self.Nx] = 0.0
                self.B[ii] = self.pBC[1]
                if ii < self.Ni-1:
                    self.A[ii,ii+1] = 0.0

    def solve(self):
        """ Solve the linear system Ax = B """
        self.x = spsolve(self.A, self.B)
        for ii in range(self.Ni):
            self.un[ii] = self.up[ii]
            self.vn[ii] = self.vp[ii]
            self.pn[ii] = self.pp[ii]
            self.pp[ii] = self.x[ii]
