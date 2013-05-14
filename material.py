#!/usr/bin/python
import numpy as np
class material(object):
#   ------- vacuum configuration

    def __init__(self,parameters,t=0,**kargs):
        self.vacuum_config = parameters.vacuum_config    
        cls_dict = self.__class__.__dict__
        for key in kargs:
            assert key in cls_dict
        self.__dict__.update(kargs)
        self.bkg_er = parameters.bkg_er
        self.bkg_mr = parameters.bkg_mr
        self.n = np.sqrt(self.bkg_er*self.bkg_mr)
        self.mat_shape = parameters.mat_shape
        self.vacuum_config = parameters.vacuum_config
        self.t_aux = t

    def vacuum(self):
        if self.vacuum_config=='real':
            self.eo = 8.854187817e-12            # vacuum permittivity   - [F/m]
            self.mo = 4e-7*np.pi                 # vacuum peremeability  - [V.s/A.m]
            self.co = 1/np.sqrt(self.eo*self.mo)           # vacuum speed of light - [m/s]
            self.zo = np.sqrt(self.eo/self.mo)
        else:
            self.eo = 1            # vacuum permittivity   - [F/m]
            self.mo = 1                 # vacuum peremeability  - [V.s/A.m]
            self.co = 1/np.sqrt(self.eo*self.mo)           # vacuum speed of light - [m/s]
            self.zo = np.sqrt(self.eo/self.mo)
        
        return self
    
    # def etar(da,ddx,ddy,t=0):
    #     """
    #     eta = etar(num_aux,xi,xf,yi,yf,ddx,ddy)

    #     Sets the auxiliary arrays for permittivity and permeability.

    #     Implemented mappings

    #     ..gaussian1dx:  stationary and moving gaussian shape for eps and mu
    #     ..homogeneous:  homogeneous refractive index in eps and mu
    #     ..interface:    simple interface (jump) acroos the 2d domain
    #     ..interfacex:   simple interface (jump) 1D in x-direction
    #     ..interfacey:   ibid in y-direction
    #     ..multilayer:   2D multilayers in x or y direction.


    #     y,x are the point coordinates of the grid.
    #     t is the time coordinate

    #     on output aux holds:

    #                                EM equivalent

    #          idim   curvilinear  |   TE      TM
    #          0:     eta_1        |   mu1     eps1
    #          1:     eta_2        |   mu2     eps2
    #          2:     eta_3        |   eps3    mu3

    #     """
    #     nx, ny = da.getSizes()
    #     (xi, xf), (yi, yf) = da.getRanges()
    #     X = np.linspace(xi,xf,xf-xi)*ddx
    #     Y = np.linspace(yi,yf,yf-yi)*ddy
    #     y,x = np.meshgrid(Y,X)
    #     eta = np.empty( [3,len(X),len(Y)], order='F')

    #     if mat_shape=='gaussian1dx':
    #         u_x_e = x - rip_vx_e*t - rip_xoff_e
    #         u_x_m = x - rip_vx_m*t - rip_xoff_m
    #         u_y_e = y - rip_vy_e*t - rip_yoff_e
    #         u_y_m = y - rip_vy_m*t - rip_yoff_m

    #         u_e = (u_x_e/rip_xsig_e)**2 + (u_y_e/rip_ysig_e)**2
    #         u_m = (u_x_m/rip_xsig_m)**2 + (u_y_m/rip_ysig_m)**2

    #         eta[0,:,:] = d_e*np.exp(-u_e) + bkg_er
    #         eta[1,:,:] = d_e*np.exp(-u_e) + bkg_er
    #         eta[2,:,:] = d_m*np.exp(-u_m) + bkg_mr
    #     elif mat_shape=='homogeneous':
    #         eta[0,:,:] = bkg_er
    #         eta[1,:,:] = bkg_er
    #         eta[2,:,:] = bkg_mr
    #     elif mat_shape=='vacuum':
    #         eta[0,:,:] = 1.0
    #         eta[1,:,:] = 1.0
    #         eta[2,:,:] = 1.0
    #     elif mat_shape=='interfacex':
    #         eta[0,:,:] = 1*(x<x_change) + 4*(x>=x_change)
    #         eta[1,:,:] = 1*(x<x_change) + 4*(x>=x_change)
    #         eta[2,:,:] = 1*(x<x_change) + 4*(x>=x_change)
    #     elif mat_shape=='interfacey':
    #         yy = y_upper-y_lower
    #         eta[0,:,:] = 1*(y<yy/2) + 4*(x>=yy/2)
    #         eta[1,:,:] = 1*(y<yy/2) + 4*(x>=yy/2)
    #         eta[2,:,:] = 1*(y<yy/2) + 4*(x>=yy/2)
    #     elif mat_shape=='multilayer':
    #         for n in range(0,N_layers):
    #             yi = n*tlp
    #             for m in range(0,n_layers):
    #                 if m==0:
    #                     eta[0,:,:] = layers[m,0]*(yi<y)*(y<=yi+layers[m,3])
    #                     eta[1,:,:] = layers[m,1]*(yi<y)*(y<=yi+layers[m,3])
    #                     eta[2,:,:] = layers[m,1]*(yi<y)*(y<=yi+layers[m,3])
    #                 else:
    #                     eta[0,:,:] = layers[m,0]*(yi+layers[m-1,3]<y)*(y<=yi+layers[m,3])
    #                     eta[1,:,:] = layers[m,1]*(yi+layers[m-1,3]<y)*(y<=yi+layers[m,3])
    #                     eta[2,:,:] = layers[m,1]*(yi+layers[m-1,3]<y)*(y<=yi+layers[m,3])


    #         eta[0,:,:] = layers[0,0]*(N_layers*tlp<y)*(y<=N_layers*tlp+layers[0,3])
    #         eta[1,:,:] = layers[0,1]*(N_layers*tlp<y)*(y<=N_layers*tlp+layers[0,3])
    #         eta[2,:,:] = layers[0,1]*(N_layers*tlp<y)*(y<=N_layers*tlp+layers[0,3])

    #     return eta
