#!/usr/bin/env python
# encoding: utf-8
# -------- GLOBAL SCALAR DEFINITIONS -----------------------------
# ======== all definitions are in m,s,g unit system.
import numpy as np 
class Params(object):
    """
    simulation parameters for fdtd simulation
    """

    n_frames = 30
    max_steps = 250000
    # ....... dimensions .............................................
    num_dim = 2
    dim_x_lower = 0.0e-6
    dim_x_upper = 100e-6                    # lenght [m]
    dim_y_lower = 0.0e-6
    dim_y_upper = 1.0e-6                    # notice that for multilayer this is value will be over-written
    dim_x_resolution = 60
    dim_y_resolution = 20
    # ....... aux config .............................................
    num_aux = 3
    vacuum_config = 'real'

    material_shape = 'homogeneous'
    # background configuration
    material_bkg_er = 1.5
    material_bkg_mr = 1.5

    # ........ Boundary settings .................

    bc_x_lower = 'scattering'
    bc_x_upper = 'none'
    bc_y_lower = 'scattering'
    bc_y_upper = 'none'

    aux_bc_x_lower = 'scattering'
    aux_bc_x_upper = 'pml'
    aux_bc_y_upper = 'metallic'
    aux_bc_y_upper = 'metallic'

    # parameters needed for pml calculation
    num_pml = 8
    pml_norder = 3
    pml_Ro = 1.0e-6

    # ........ excitation - initial conditoons .......................

    source_type  = 'plane'
    source_lambda  = 1e-6             # wavelength
    source_t_sig = 1.0*source_lambda          # width in time (pulse only)
    source_x_sig = 1.0*source_lambda          # width in the x-direction (pulse)
    source_y_sig = dim_y_upper-dim_y_lower
    source_toff  = 0.0                  # offset in time
    source_xoff  = 0.0                  # offset in the x-direction
    source_yoff  = dim_y_upper/2            # offset in the y-direction# frequency
    source_k        = 2.0*np.pi/source_lambda
    source_amp_Ex   = 0.
    source_amp_Ey   = 1.
    source_amp_Hz   = 1.


    def __init__(self, **kargs):
        cls_dict = self.__class__.__dict__
        for key in kargs:
            assert key in cls_dict
        self.__dict__.update(kargs)
        self.set_vacuum(vc_config=self.vacuum_config)
        self.set_material(shape=self.material_shape)
        self.set_source()
        self.set_dims(shape=self.material_shape)

    def set_source(self):
        self.source_omega = 2.0*np.pi*self.co/self.source_lambda 
        self.source_vr = 1./self.material_bkg_n
        self.source_v = self.co*self.source_vr
        source_vx = self.source_v
        source_vy = 0.0
        source_kx = self.source_k
        source_ky = 0.0        

    def set_dims(self,shape=material_shape):
        """
        calculates y_upper and my based on the material shape
        """
        self.dim_mx = np.floor(40*(self.dim_x_upper - self.dim_x_lower)/self.source_lambda)
        self.dim_ddx = (self.dim_x_upper-self.dim_x_lower)/self.dim_mx
        self.t_final = (self.dim_x_upper - self.dim_x_lower)/self.source_v
        self.dxdt = 1
        self.dydt = 1
        if shape=='multilayer':
            self.dim_y_upper = self.material_N_layers*np.sum(self.material_layers[:,3]) + self.material_layers[0,3]
            self.dim_tlp = np.sum(self.material_layers[:,3])
            self.dim_mlp = np.floor(self.tlp/1e-9)
            self.dim_my = np.floor(self.dim_y_resolution*(self.dim_y_upper - self.dim_y_lower)/1e-9)
        else:
            self.dim_my = np.floor(self.dim_y_resolution*(self.dim_y_upper - self.dim_y_lower)/self.source_lambda)

        self.dim_ddy = (self.dim_y_upper-self.dim_y_lower)/self.dim_my
        self.ddt = 0.90/(self.co*np.sqrt(1.0/(self.dim_ddx**2) + 1.0/(self.dim_ddy**2)))
        self.dt = self.ddt
        return self

    def set_material(self,shape=material_shape,vx=0,vy=0,x_offset=10e-6,y_offset=0,sigma=10e-6,n_layers=2,N_layers=5):
        """
        :Implemented mappings
        
        ..gaussian1dx:  stationary and moving gaussian shape for eps and mu
        ..homogeneous:  homogeneous refractive index in eps and mu
        ..interface:    simple interface (jump) acroos the 2d domain
        ..interfacex:   simple interface (jump) 1D in x-direction
        ..interfacey:   ibid in y-direction
        ..multilayer:   2D multilayers in x or y direction.
        
        """
        self.material_bkg_n = np.sqrt(self.material_bkg_er*self.material_bkg_mr)
        # if interface declare position
        if shape=='homogeneous':
            pass
        elif shape=='interfacex':
            self.material_x_change = self.dim_x_upper/2
        elif shape=='gaussian1dx' or shape=='gaussian':
            # set moving refractive index parameters
            self.material_rip_vx_e = vx*self.co   # replace here the value of x
            self.material_rip_vx_m    = self.material_rip_vx_e
            self.material_rip_vy_e    = vy*self.co
            self.material_rip_vy_m    = self.material_rip_vy_e
            # offset values
            self.material_rip_xoff_e  = x_offset
            self.material_rip_xoff_m  = self.material_rip_xoff_e
            self.material_rip_yoff_e  = 0.0
            self.material_rip_yoff_m  = self.material_rip_yoff_e
            # dispersion
            self.material_rip_xsig_e  = sigma
            self.material_rip_xsig_m  = self.material_rip_xsig_e
            self.material_rip_ysig_e  = .9*self.dim_y_upper/2
            self.material_rip_ysig_m  = self.material_rip_ysig_e
            self.material_s_x_e       = self.material_rip_xsig_e**2
            self.material_s_x_m       = self.material_rip_xsig_m**2
            self.material_s_y_e       = self.material_rip_ysig_e**2
            self.material_s_y_m       = self.material_rip_ysig_m**2
            # perturbation length
            self.material_prip        = 0.1
            self.material_deltan      = self.material_prip*(self.material_bkg_n) # assumes epsilon = mu
            self.material_d_e         = self.material_deltan #*(2.0*1.5+deltan)
            self.material_d_m         = self.material_deltan #*(2.0*1.5+deltan)
        elif shape=='multilayer':
            # multilayered definition
            self.material_n_layers = n_layers
            self.material_layers = np.zeros([n_layers,7]) # _layer:  eps mu N t chi2e chi2m chi3e chi3m
            self.material_layers[0,0] = 1.5
            self.material_layers[0,1] = 1.5
            self.material_layers[0,2] = 10
            self.material_layers[0,3] = 15e-9
            self.material_layers[1,0] = 2.5
            self.material_layers[1,1] = 2.5
            self.material_layers[1,2] = self.material_layers[0,2] - 1
            self.material_layers[1,3] = 50e-9
            self.material_N_layers = N_layers
        
        return self

    def set_vacuum(self,vc_config=vacuum_config):
        if vc_config=='real':
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
    