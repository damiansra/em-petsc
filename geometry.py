import numpy as np
# Dimension object based on PyClaw dimension object
# Default mapc2p function

## try:
##    from petsc4py.PETSc import DA
## except ImportError:
##    from fakeda import DA


class Dimension(object):
    r"""
    basic class to represent a dimension
    
    x = geometry.dimension(name,x_lower,x_upper,n)

    where:

     - *name*   - (string) string Name of dimension
     - *lower*  - (float) Lower extent of dimension
     - *upper*  - (float) Upper extent of dimension
     - *n*      - (int) Number of points

    Example:

    >>> from geometry import dimension
    >>> x = dimension('x',0.,1.,100)
    >>> print x
    Dimension x:  (num_cells,delta,[lower,upper]) = (100,0.01,[0.0,1.0])
    >>> x.name
    'x'
    >>> x.num_points
    100
    >>> x.delta
    0.01
    >>> x.grid_points[0]
    0.0
    >>> x.grid_points[1]
    0.01
    >>> x.grid_points[-1]
    1.0
    >>> x.centers[-1]
    0.995
    >>> len(x.centers)
    100
    >>> len(x.edges)
    101
    """
    
    # ========== Property Definitions ========================================
    @property
    def delta(self):
        r"""(float) - distance between grid points"""
        return (self.upper-self.lower) / float(self.num_points)
    @property
    def grid_points(self):
        r"""(ndarrary(:)) - Location of grid points"""
        if self._grid_points is None:
            self._grid_points = np.linspace(self.lower,self.upper,self.num_points+1)
        return self._grid_points
    _grid_points = None
    @property
    def centers(self):
        r"""(ndarrary(:)) - Location of all cell center coordinates
        for this dimension"""
        if self._centers is None:
            self._centers = np.linspace(self.lower+self.delta/2,self.upper-self.delta/2,self.num_points)
        return self._centers
    _centers = None

    def __init__(self, *args, **kargs):
        r"""
        Creates a Dimension object
        
        See :class:`Dimension` for full documentation
        """
        
        # ========== Class Data Attributes ===================================
        self.name = 'x'
        r"""(string) Name of this coordinate dimension (e.g. 'x')"""
        self.num_points = None
        r"""(int) - Number of cells in this dimension :attr:`units`"""
        self.lower = 0.0
        r"""(float) - Lower computational dimension extent"""
        self.upper = 1.0
        r"""(float) - Upper computational dimension extent"""
        
        # Parse args
        if isinstance(args[0],float):
            self.lower = float(args[0])
            self.upper = float(args[1])
            self.num_points = int(args[2]+1)
            self.num_cells = int(args[2])
        elif isinstance(args[0],basestring):
            self.name = args[0]
            self.lower = float(args[1])
            self.upper = float(args[2])
            self.num_points = int(args[3]+1)
            self.num_cells = int(args[3])
        else:
            raise Exception("Invalid initializer for dimension.")
        
        for (k,v) in kargs.iteritems():
            setattr(self,k,v)

    def __str__(self):
        output = "Dimension %s" % self.name
        if self.units:
            output += " (%s)" % self.units
        output += ":  (num_cells,delta,[lower,upper]) = (%s,%s,[%s,%s])" \
            % (self.num_points,self.delta,self.lower,self.upper)
        return output

class Grid(object):
    r"""
    Basic representation of a single grid in
    
    :Dimension information:
    
        Each dimension has an associated name with it that can be accessed via
        that name such as ``grid.x.num_cells`` which would access the x dimension's
        number of cells.
    
    :Properties:

        If the requested property has multiple values, a list will be returned
        with the corresponding property belonging to the dimensions in order.
         
    :Initialization:
    
        Input:
         - *dimensions* - (list of :class:`Dimension`) Dimensions that are to 
           be associated with this grid
            
        Output:
         - (:class:`grid`) Initialized grid object

    \grid is usually constructed from a tuple of PyClaw Dimension objects:

    >>> from clawpack.pyclaw.geometry import Dimension, Grid      
    >>> x = Dimension('x',0.,1.,10)
        >>> y = Dimension('y',-1.,1.,25)
        >>> grid = Grid((x,y))
        >>> print grid
        Dimension x:  (num_cells,delta,[lower,upper]) = (10,0.1,[0.0,1.0])
        Dimension y:  (num_cells,delta,[lower,upper]) = (25,0.08,[-1.0,1.0])
        >>> grid.num_dim
        2
        >>> grid.num_cells
        [10, 25]
        >>> grid.lower
        [0.0, -1.0]
        >>> grid.delta
        [0.1, 0.08]

    A grid can be extended to higher dimensions using the add_dimension() method:

        >>> z=Dimension('z',-2.0,2.0,21)
        >>> grid.add_dimension(z)
        >>> grid.num_dim
        3
        >>> grid.num_cells
        [10, 25, 21]
        >>> grid.c_edges[0][0,0,0]
        0.0
        >>> grid.c_edges[1][0,0,0]
        -1.0
        >>> grid.c_edges[2][0,0,0]
        -2.0
    """

    # ========== Property Definitions ========================================
    @property
    def num_dim(self):
        r"""(int) - Number of dimensions"""
        return len(self._dimensions)
    @property
    def dimensions(self):
        r"""(list) - List of :class:`Dimension` objects defining the 
                grid's extent and resolution"""
        return [getattr(self,name) for name in self._dimensions]
    @property
    def num_cells(self): 
        r"""(list) - List of the number of cells in each dimension"""
        return self.get_dim_attribute('num_cells')
    @property
    def lower(self):
        r"""(list) - Lower coordinate extents of each dimension"""
        return self.get_dim_attribute('lower')
    @property
    def upper(self):
        r"""(list) - Upper coordinate extends of each dimension"""
        return self.get_dim_attribute('upper')
    @property
    def delta(self):
        r"""(list) - List of computational cell widths"""
        return self.get_dim_attribute('delta')
    @property
    def grid_points(self):
        r"""(list) - List of grid_points"""
        return self.get_dim_attribute('grid_points')
    @property
    def centers(self):
        r"""(list) - List of center coordinate arrays"""
        return self.get_dim_attribute('centers')
    @property
    def num_points(self):
        "number of points in each dimension"
        return self.get_dim_attribute('num_points')
    @property
    def curvilinear(self):
        r"""(list of ndarray(...)) - List containing the arrays locating
                  the physical locations of cell centers, see 
                  :meth:`compute_p_centers` for more info."""
        self.compute_curvilinear(self)
        return self._curvilinear
    _curvilinear = None

       
    
    # ========== Class Methods ===============================================
    def __init__(self,dimensions):
        r"""
        Instantiate a Grid object
        
        See :class:`Grid` for more info.
        """
        
        # ========== Attribute Definitions ===================================
        # Dimension parsing
        if isinstance(dimensions,Dimension):
            dimensions = [dimensions]
        self._dimensions = []
        for dim in dimensions:
            self.add_dimension(dim)

        #super(Grid,self).__init__()
    
    # ========== Dimension Manipulation ======================================
    def add_dimension(self,dimension):
        r"""
        Add the specified dimension to this patch
        
        :Input:
         - *dimension* - (:class:`Dimension`) Dimension to be added
        """

        # Add dimension to name list and as an attribute
        if dimension.name in self._dimensions:
            raise Exception('Unable to add dimension. A dimension'\
             +' of the same name: {name}, already exists.'\
             .format(name=dimension.name))

        self._dimensions.append(dimension.name)
        setattr(self,dimension.name,dimension)
        
        
    def get_dim_attribute(self,attr):
        r"""
        Returns a tuple of all dimensions' attribute attr
        """
        return [getattr(getattr(self,name),attr) for name in self._dimensions]
    
    
    # ========== Copy functionality ==========================================
    def __copy__(self):
        return self.__class__(self)
        
    # ========== Grid Operations =============================================
    def curvilinear(self):
        r"""
        grid transformation in curvilinear coordinates
        """
        pass
    def cart_to_curv(self):
        r"""
        Apply transformation from cartesian to curvilinar coordinates as specificied in function curvilinear
        """
        pass
    def curv_to_cart(self):
        r"""
        back transformation from curvilinear to cartesian
        """

class PETSc_grid(Grid,object):
    r"""
    creates the petsc_grid object
    """

    @property
    def staggered_edge(self,nvec):
        r"""
        returns the value of the staggered coordinates in the local grid for the given nvec.
        """
        pass
    @property
    def staggered_centers(self):
        r"""
        values of the grid at staggered coordinate values
        """
    @property
    def grid_points(self):
        r"""
        returns the grid points 
        """
        pass
    @property
    def grid_cell_centers(self):
        r"""
        return the centers of the inner cells of the grid
        """
        pass
    @property
    def grid_edge(self):
        pass
    @property
    def grid_with_ghosts(self):
        pass
