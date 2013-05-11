# excitation calculations
# ........ excitation - initial conditoons .......................

ex_type  = 'plane'
alambda  = 1e-6             # wavelength
ex_t_sig = 1.0*alambda          # width in time (pulse only)
ex_x_sig = 1.0*alambda          # width in the x-direction (pulse)
ex_y_sig = y_upper-y_lower
ex_toff  = 0.0                  # offset in time
ex_xoff  = 0.0                  # offset in the x-direction
ex_yoff  = y_upper/2            # offset in the y-direction
omega    = 2.0*np.pi*co/alambda # frequency
k        = 2.0*np.pi/alambda
amp_Ex   = 0.
amp_Ey   = 1.
amp_Hz   = 1.
