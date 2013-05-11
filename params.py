
class Params(object):

    n_frames = 30
    x_lower = 0.0e-6
    x_upper = 100e-6
    y_lower = 0.0e-6
    y_upper = 1.0e-6
    num_aux = 3
    vacuum_config = 'real'
    bkg_er = 1.5
    bkg_mr = 1.5

    def __init__(self, **kargs):
        cls_dict = self.__class__.__dict__
        for key in kargs:
            assert key in cls_dict
        self.__dict__.update(kargs)
