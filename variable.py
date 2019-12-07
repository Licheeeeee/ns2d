""" Create the data structure that stores the variables """
import numpy as np

class Variable():
    def __init__(self, ic, setting):
        self.Nx = setting['dim'][0]
        self.Ny = setting['dim'][1]
        self.dx = setting['dim'][0]
        self.dy = setting['dim'][1]
        self.Ni = self.Nx * self.Ny
        self.Nt = (self.Nx+2) * (self.Ny+2)
        self.dt = setting['T'][0]
        self.rho = setting['rho']
        #   Variables that remain constants for now
        self.V = self.dx * self.dy
        self.Ax = self.dy
        self.Ay = self.dx
        #   Variables on cell centers
        self.pn = ic['p'] * np.ones((self.Nt,1))
        self.pp = ic['p'] * np.ones((self.Nt,1))
        self.B = np.zeros((self.Ni,1))
        #   Variables on cell faces
        self.un = np.zeros((self.Nt,1))
        self.up = np.zeros((self.Nt,1))
        self.vn = np.zeros((self.Nt,1))
        self.vp = np.zeros((self.Nt,1))
        self.Gxp = (self.dt * self.Ax / (self.rho * self.V)) * np.ones((self.Ni,1))
        self.Gxm = (self.dt * self.Ax / (self.rho * self.V)) * np.ones((self.Ni,1))
        self.Gyp = (self.dt * self.Ax / (self.rho * self.V)) * np.ones((self.Ni,1))
        self.Gym = (self.dt * self.Ax / (self.rho * self.V)) * np.ones((self.Ni,1))
        self.Gct = self.Gxp + self.Gxm + self.Gyp + self.Gym
