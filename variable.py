""" Create the data structure that stores the variables """
import numpy as np
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve
from decimal import *

class Variable():
    def __init__(self, map, setting):
        self.Nx = setting['dim'][0]
        self.Ny = setting['dim'][1]
        self.dx = setting['delta'][0]
        self.dy = setting['delta'][1]
        self.Ni = self.Nx * self.Ny
        self.Nt = (self.Nx+2) * (self.Ny+2)
        self.dt = setting['T'][0]
        self.nu = setting['nu']
        self.ka = setting['kappa']
        self.pBC = setting['pBC']
        self.vBC = setting['vBC']
        self.sBC = setting['sBC']
        self.sIC = setting['sIC']
        #   Variables that remain constants for now
        self.V = self.dx * self.dy
        self.Ax = self.dy
        self.Ay = self.dx
        #   Variables on cell centers
        self.pn = self.pBC[1] * np.ones(self.Nt)
        self.pp = self.pBC[1] * np.ones(self.Nt)
        self.sn = self.sIC * np.ones(self.Nt)
        self.sp = self.sIC * np.ones(self.Nt)
        self.rho = setting['rho'] * np.ones(self.Ni)
        for ii in range(self.Ni):
            if map.col[ii] <= int(self.Ny/5) and map.block[ii] == 0:
                self.sn[ii] = self.sBC
                self.sp[ii] = self.sBC
        self.B = np.zeros((self.Ni,1))
        #   Variables on cell faces
        self.un = np.zeros(self.Nt)
        self.up = np.zeros(self.Nt)
        self.vn = np.zeros(self.Nt)
        self.vp = np.zeros(self.Nt)
        for ii in range(self.Nx):
            # self.vn[map.icjM[ii]] = self.vBC
            # self.vp[map.icjM[ii]] = self.vBC
            self.sn[map.icjM[ii]] = self.sBC
            self.sp[map.icjM[ii]] = self.sBC
        #   Matrix coefficients
        self.Gx = (self.dt * self.Ax)**2.0 / (self.rho * self.V) * np.ones((self.Ni))
        self.Gy = (self.dt * self.Ay)**2.0 / (self.rho * self.V) * np.ones((self.Ni))
        self.Gc = (3e7 * self.V / self.rho) * np.ones((self.Ni))
        # self.Gc = 2.0 * self.Gx + 2.0 * self.Gy
        # self.Gc = np.zeros((self.Ni))
        for ii in range(self.Ni):
            if map.block[ii] == 1:
                self.Gx[ii] = 0.0
                self.Gy[ii] = 0.0
                if map.iMjc[ii] < self.Ni:
                    self.Gx[map.iMjc[ii]] = 0.0
                if map.icjM[ii] < self.Ni:
                    self.Gy[map.icjM[ii]] = 0.0
            if map.iPjc[ii] >= self.Ni:
                self.Gx[ii] = 0.0
        for ii in range(self.Ni):
            self.Gc[ii] += self.Gx[ii] + self.Gy[ii]
            # print('ii : ',ii,' Gc,Gx,Gy=',self.Gc[ii],self.Gx[ii],self.Gy[ii])
            if map.iMjc[ii] < self.Ni:
                self.Gc[ii] += self.Gx[map.iMjc[ii]]
            if map.icjM[ii] < self.Ni:
                self.Gc[ii] += self.Gy[map.icjM[ii]]
            # if map.icjP[ii] >= self.Ni:
            #     self.Gc[ii] -= self.Gy[ii]
            # print(self.Gc[ii])

    def explicitSource(self, map, setting):
        """ Calculate the explicit terms """
        self.Ex = np.zeros(self.Nt)
        self.Ey = np.zeros(self.Nt)
        for ii in range(self.Ni):
            #   Transverse velocity
            if map.icjP[ii] < self.Ni:
                indU = map.iMjc[map.icjP[ii]]
            elif map.iMjc[ii] < self.Ni:
                indU = map.icjP[map.iMjc[ii]]
            else:
                indU = map.icjP[ii]
            if map.iPjc[ii] < self.Ni:
                indV = map.icjM[map.iPjc[ii]]
            elif map.icjM[ii] < self.Ni:
                indV = map.iPjc[map.icjM[ii]]
            else:
                indV = map.icjM[ii]

            u4 = np.array((self.un[ii],self.un[map.iMjc[ii]],self.un[map.icjP[ii]],self.un[indU]))
            is0 = u4 != 0
            if np.nansum(is0) == 0:
                ut = 0.0
            else:
                ut = np.nansum(u4) / np.nansum(is0)
            v4 = np.array((self.vn[ii],self.vn[map.iPjc[ii]],self.vn[map.icjM[ii]],self.vn[indV]))
            is0 = v4 != 0
            if np.nansum(is0) == 0:
                vt = 0.0
            else:
                vt = np.nansum(v4) / np.nansum(is0)
            #   Explicit velocity
            self.Ex[ii] = self.un[ii]
            self.Ey[ii] = self.vn[ii]
            #   Advection terms
            self.Ex[ii] -= (0.5 * self.dt / self.dx) * ((self.un[ii] + abs(self.un[ii])) * (self.un[ii] - self.un[map.iMjc[ii]]) + (self.un[ii] - abs(self.un[ii])) * (self.un[ii] - self.un[map.iPjc[ii]]))
            self.Ex[ii] -= (0.5 * self.dt / self.dy) * ((vt + abs(vt)) * (self.un[ii] - self.un[map.icjM[ii]]) + (vt - abs(vt)) * (self.un[ii] - self.un[map.icjP[ii]]))
            self.Ey[ii] -= (0.5 * self.dt / self.dx) * ((ut + abs(ut)) * (self.vn[ii] - self.vn[map.iMjc[ii]]) + (ut - abs(ut)) * (self.vn[ii] - self.vn[map.iPjc[ii]]))
            self.Ey[ii] -= (0.5 * self.dt / self.dy) * ((self.vn[ii] + abs(self.vn[ii])) * (self.vn[ii] - self.vn[map.icjM[ii]]) + (self.vn[ii] - abs(self.vn[ii])) * (self.vn[ii] - self.vn[map.icjP[ii]]))
            #   Diffusion terms
            self.Ex[ii] += (self.dt * self.nu * self.Ax / (self.dx * self.V)) * (self.un[map.iPjc[ii]] - 2.0*self.un[ii] + self.un[map.iMjc[ii]])
            self.Ex[ii] += (self.dt * self.nu * self.Ay / (self.dy * self.V)) * (self.un[map.icjP[ii]] - 2.0*self.un[ii] + self.un[map.icjM[ii]])
            self.Ey[ii] += (self.dt * self.nu * self.Ax / (self.dx * self.V)) * (self.vn[map.iPjc[ii]] - 2.0*self.vn[ii] + self.vn[map.iMjc[ii]])
            self.Ey[ii] += (self.dt * self.nu * self.Ay / (self.dy * self.V)) * (self.vn[map.icjP[ii]] - 2.0*self.vn[ii] + self.vn[map.icjM[ii]])
            #   Zero on blocked faces
            if self.Gx[ii] == 0.0:
                self.Ex[ii] = 0.0
            if self.Gy[ii] == 0.0:
                self.Ey[ii] = 0.0
        for ii in range(self.Nx):
            jj = map.icjM[ii]
            self.Ey[jj] = self.vn[jj]
            self.Ey[jj] += (self.dt * self.nu * self.Ay / (self.dy * self.V)) * (self.vn[ii] - self.vn[jj])
        for ii in range(self.Ny):
            jj = ii * self.Nx
            self.Ex[map.iMjc[jj]] = 0.0

    def matrixRHS(self, map, setting):
        """ Calculate right-hand-side of the linear system """
        self.B = np.zeros(self.Ni)
        for ii in range(self.Ni):
            self.B[ii] += (3e7 * self.V / self.rho[ii]) * self.pn[ii]
            self.B[ii] += self.dt * (self.Ax * (-self.Ex[ii] + self.Ex[map.iMjc[ii]]) \
                - self.Ay * (self.Ey[ii] + self.Ey[map.icjM[ii]]))

    def updateVelocity(self, map, setting):
        """ Update velocity based on pressure """
        for ii in range(self.Ni):
            # print('ii, P, Pip, pjp = ',ii, self.pp[ii], self.pp[map.iPjc[ii]], self.pp[map.icjP[ii]])
            self.up[ii] = -(self.dt*self.Ax / (self.rho[ii]*self.V)) * (self.pp[map.iPjc[ii]] - self.pp[ii]) + self.Ex[ii]
            self.vp[ii] = -(self.dt*self.Ay / (self.rho[ii]*self.V)) * (self.pp[map.icjP[ii]] - self.pp[ii]) + self.Ey[ii]

    def enforcePressureBC(self, map, setting):
        """ Enforce BC on pressure """
        for ii in range(self.Nx):
            jj = self.Nx * (self.Ny-1) + ii
            self.pp[map.icjM[ii]] = self.pp[ii]
            # self.pp[jj] = self.pBC
            self.pp[map.icjP[jj]] = self.pp[jj]
        for ii in range(self.Ny):
            jj = ii * self.Nx
            kk = (ii+1) * self.Nx - 1
            self.pp[map.iMjc[jj]] = self.pp[jj]
            self.pp[map.iPjc[kk]] = self.pp[kk]

    def enforceVelocityBC(self, map, setting):
        """ Enforce BC on velocity """
        #   BC on blocked cells
        for ii in range(self.Ni):
            if map.block[ii] == 1:
                self.up[ii] = 0.0
                self.up[map.iMjc[ii]] = 0.0
                self.vp[ii] = 0.0
                self.vp[map.icjM[ii]] = 0.0
        #   BC on ghost cells
        for ii in range(self.Nx):
            jj = self.Nx * (self.Ny-1) + ii
            self.up[map.icjM[ii]] = self.up[ii]
            self.up[map.icjP[jj]] = self.up[jj]
            # self.vp[map.icjM[ii]] = self.vBC
            self.vp[map.icjM[ii]] = self.vp[ii]
            self.vp[jj] = self.vp[map.icjM[jj]] + (self.up[map.iMjc[jj]]-self.up[jj]) * self.Ax / self.Ay
        for ii in range(self.Ny):
            jj = ii * self.Nx
            kk = (ii+1) * self.Nx - 1
            self.up[kk] = 0.0
            self.up[map.iMjc[jj]] = 0.0
            self.up[map.iPjc[kk]] = self.up[kk]
            self.vp[map.iMjc[jj]] = self.vp[jj]
            self.vp[map.iPjc[kk]] = self.vp[kk]

    def updateVariables(self, setting):
        for ii in range(self.Nt):
            self.un[ii] = self.up[ii]
            self.vn[ii] = self.vp[ii]
            self.pn[ii] = self.pp[ii]
            self.sn[ii] = self.sp[ii]

    def buildMatrix(self, map, setting):
        """ Build sparse matrix A """
        self.A = lil_matrix((self.Ni, self.Ni), dtype=float)
        for ii in range(self.Ni):
            self.A[ii,ii] = self.Gc[ii]
            #   For Neumann boundary
            if map.iMjc[ii] >= self.Ni:
                # self.A[ii,ii] -= self.Gx[ii]
                self.A[ii,ii] = self.roundTo(self.A[ii,ii])
            else:
                self.A[ii,ii-1] = -self.Gx[map.iMjc[ii]]
            if map.iPjc[ii] >= self.Ni:
                # self.A[ii,ii] -= self.Gx[ii]
                self.A[ii,ii] = self.roundTo(self.A[ii,ii])
            else:
                self.A[ii,ii+1] = -self.Gx[ii]
            #   For Direchlet boundary
            if map.icjM[ii] >= self.Ni:
                self.A[ii,ii] = self.roundTo(self.A[ii,ii])
            elif map.icjM[ii] < self.Nx:
                self.B[ii] += self.Gy[map.icjM[ii]] * self.pBC[0]
            else:
                self.A[ii,ii-self.Nx] = -self.Gy[map.icjM[ii]]
            if map.icjP[ii] >= self.Ni:
                self.A[ii,ii] = self.roundTo(self.A[ii,ii])
            elif map.icjP[ii] >= self.Ni-self.Nx:
                self.B[ii] += self.Gy[ii] * self.pBC[1]
            else:
                self.A[ii,ii+self.Nx] = -self.Gy[ii]
        #   Enforce pressure BC
        for ii in range(self.Ni):
            if map.icjP[ii] >= self.Ni:
                self.A[ii,ii] = 1.0
                self.A[ii,ii-self.Nx] = 0.0
                self.B[ii] = self.pBC[1]
                if ii > 0:
                    self.A[ii,ii-1] = 0.0
                if ii < self.Ni-1:
                    self.A[ii,ii+1] = 0.0
            if map.icjM[ii] >= self.Ni:
                self.A[ii,ii] = 1.0
                self.A[ii,ii+self.Nx] = 0.0
                self.B[ii] = self.pBC[0]
                if ii > 0:
                    self.A[ii,ii-1] = 0.0
                if ii < self.Ni-1:
                    self.A[ii,ii+1] = 0.0
        #    Enforce blocked BC
        for ii in range(self.Ni):
            if map.block[ii] == 1:
                self.A[ii,ii] = 1.0
                self.B[ii] = 1.0
                if ii-1 >= 0:
                    self.A[ii,ii-1] = 0.0
                if ii+1 < self.Ni:
                    self.A[ii,ii+1] = 0.0
                if ii-self.Nx >= 0:
                    self.A[ii,ii-self.Nx] = 0.0
                if ii+self.Nx < self.Ni:
                    self.A[ii,ii+self.Nx] = 0.0

    def solve(self):
        """ Solve the linear system Ax = B """
        for ii in range(self.Ni):
            self.B[ii] = self.roundTo(self.B[ii])
        self.x = spsolve(self.A, self.B)
        # print('Ex,Ey,B,p,px,py = ',self.Ex[36],self.Ey[27],self.B[37],self.x[37],self.x[36],self.x[27])
        for ii in range(self.Ni):
            self.pp[ii] = self.x[ii]
            self.pp[ii] = self.roundTo(self.x[ii])
            # if ii % self.Nx == 5:
            #     print(self.x[ii])
            # print(self.pp[ii])

    def roundTo(self, n, m=8):
        return round(n, m)

    def transport(self, map, setting):
        """ Calculate the explicit terms """
        for ii in range(self.Ni):
            #   Transverse velocity
            if map.icjP[ii] < self.Ni:
                indU = map.iMjc[map.icjP[ii]]
            elif map.iMjc[ii] < self.Ni:
                indU = map.icjP[map.iMjc[ii]]
            else:
                indU = map.icjP[ii]
            if map.iPjc[ii] < self.Ni:
                indV = map.icjM[map.iPjc[ii]]
            elif map.icjM[ii] < self.Ni:
                indV = map.iPjc[map.icjM[ii]]
            else:
                indV = map.icjM[ii]
            ut = 0.25 * (self.up[ii] + self.up[map.iMjc[ii]] + self.up[map.icjP[ii]] + self.up[indU])
            vt = 0.25 * (self.vp[ii] + self.vp[map.iPjc[ii]] + self.vp[map.icjM[ii]] + self.vp[indV])
            #   Explicit velocity
            self.sp[ii] = self.sn[ii]
            self.sp[ii] = self.sn[ii]
            #   Advection terms
            self.sp[ii] -= (0.5 * self.dt / self.dx) * ((self.up[ii] + abs(self.up[ii])) * (self.sn[ii] - self.sn[map.iMjc[ii]]) + (self.up[ii] - abs(self.up[ii])) * (self.sn[ii] - self.sn[map.iPjc[ii]]))
            self.sp[ii] -= (0.5 * self.dt / self.dy) * ((vt + abs(vt)) * (self.sn[ii] - self.sn[map.icjM[ii]]) + (vt - abs(vt)) * (self.sn[ii] - self.sn[map.icjP[ii]]))
            self.sp[ii] -= (0.5 * self.dt / self.dx) * ((ut + abs(ut)) * (self.sn[ii] - self.sn[map.iMjc[ii]]) + (ut - abs(ut)) * (self.sn[ii] - self.sn[map.iPjc[ii]]))
            self.sp[ii] -= (0.5 * self.dt / self.dy) * ((self.vn[ii] + abs(self.vn[ii])) * (self.sn[ii] - self.sn[map.icjM[ii]]) + (self.vn[ii] - abs(self.vn[ii])) * (self.sn[ii] - self.sn[map.icjP[ii]]))
            #   Diffusion terms
            self.sp[ii] += (self.dt * self.ka * self.Ax / (self.dx * self.V)) * (self.sn[map.iPjc[ii]] - 2.0*self.sn[ii] + self.sn[map.iMjc[ii]])
            self.sp[ii] += (self.dt * self.ka * self.Ay / (self.dy * self.V)) * (self.sn[map.icjP[ii]] - 2.0*self.sn[ii] + self.sn[map.icjM[ii]])
            self.sp[ii] += (self.dt * self.ka * self.Ax / (self.dx * self.V)) * (self.sn[map.iPjc[ii]] - 2.0*self.sn[ii] + self.sn[map.iMjc[ii]])
            self.sp[ii] += (self.dt * self.ka * self.Ay / (self.dy * self.V)) * (self.sn[map.icjP[ii]] - 2.0*self.sn[ii] + self.sn[map.icjM[ii]])
        #   Boundary conditions
        for ii in range(self.Nx):
            jj = self.Nx * (self.Ny-1) + ii
            self.sp[map.icjM[ii]] = self.sBC
            self.sp[map.icjP[jj]] = self.sp[jj]
        for ii in range(self.Ny):
            jj = ii * self.Nx
            kk = (ii+1) * self.Nx - 1
            self.sp[map.iMjc[jj]] = self.sp[jj]
            self.sp[map.iPjc[kk]] = self.sp[kk]
