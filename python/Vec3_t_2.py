#from https://gist.github.com/mcleonard/5351452
import math

class Vec3_t(object):
    def __init__(self, *args):
        """ Create a Vec3_t, example: v = Vec3_t(1,2) """
        if len(args)==0: self.values = (0,0)
        else: self.values = args
        
    def norm(self):
        """ Returns the norm (length, magnitude) of the Vec3_t """
        return math.sqrt(sum( x*x for x in self ))
        
    def argument(self, radians=False):
        """ Returns the argument of the Vec3_t, the angle clockwise from +y. In degress by default, 
            set radians=True to get the result in radians. This only works for 2D Vec3_ts. """
        arg_in_rad = math.acos(Vec3_t(0, 1)*self/self.norm())
        if radians:
            return arg_in_rad
        arg_in_deg = math.degrees(arg_in_rad)
        if self.values[0] < 0: 
            return 360 - arg_in_deg
        else: 
            return arg_in_deg

    def normalize(self):
        """ Returns a normalized unit Vec3_t """
        norm = self.norm()
        normed = tuple( x / norm for x in self )
        return self.__class__(*normed)
    
    def rotate(self, theta):
        """ Rotate this Vec3_t. If passed a number, assumes this is a 
            2D Vec3_t and rotates by the passed value in degrees.  Otherwise,
            assumes the passed value is a list acting as a matrix which rotates the Vec3_t.
        """
        if isinstance(theta, (int, float)):
            # So, if rotate is passed an int or a float...
            if len(self) != 2:
                raise ValueError("Rotation axis not defined for greater than 2D Vec3_t")
            return self._rotate2D(theta)
        
        matrix = theta
        if not all(len(row) == len(self) for row in matrix) or not len(matrix)==len(self):
            raise ValueError("Rotation matrix must be square and same dimensions as Vec3_t")
        return self.matrix_mult(matrix)
        
    def _rotate2D(self, theta):
        """ Rotate this Vec3_t by theta in degrees.
            
            Returns a new Vec3_t.
        """
        theta = math.radians(theta)
        # Just applying the 2D rotation matrix
        dc, ds = math.cos(theta), math.sin(theta)
        x, y = self.values
        x, y = dc*x - ds*y, ds*x + dc*y
        return self.__class__(x, y)
        
    def matrix_mult(self, matrix):
        """ Multiply this Vec3_t by a matrix.  Assuming matrix is a list of lists.
        
            Example:
            mat = [[1,2,3],[-1,0,1],[3,4,5]]
            Vec3_t(1,2,3).matrix_mult(mat) ->  (14, 2, 26)
         
        """
        if not all(len(row) == len(self) for row in matrix):
            raise ValueError('Matrix must match Vec3_t dimensions') 
        
        # Grab a row from the matrix, make it a Vec3_t, take the dot product, 
        # and store it as the first component
        product = tuple(Vec3_t(*row)*self for row in matrix)
        
        return self.__class__(*product)
    
    def inner(self, Vec3_t):
        """ Returns the dot product (inner product) of self and another Vec3_t
        """
        if not isinstance(Vec3_t, Vec3_t):
            raise ValueError('The dot product requires another Vec3_t')
        return sum(a * b for a, b in zip(self, Vec3_t))
    
    def __mul__(self, other):
        """ Returns the dot product of self and other if multiplied
            by another Vec3_t.  If multiplied by an int or float,
            multiplies each component by other.
        """
        if isinstance(other, Vec3_t):
            return self.inner(other)
        elif isinstance(other, (int, float)):
            product = tuple( a * other for a in self )
            return self.__class__(*product)
        else:
            raise ValueError("Multiplication with type {} not supported".format(type(other)))
    
    def __rmul__(self, other):
        """ Called if 4 * self for instance """
        return self.__mul__(other)
            
    def __truediv__(self, other):
        if isinstance(other, Vec3_t):
            divided = tuple(self[i] / other[i] for i in range(len(self)))
        elif isinstance(other, (int, float)):
            divided = tuple( a / other for a in self )
        else:
            raise ValueError("Division with type {} not supported".format(type(other)))
        
        return self.__class__(*divided)
    
    def __add__(self, other):
        """ Returns the Vec3_t addition of self and other """
        if isinstance(other, Vec3_t):
            added = tuple( a + b for a, b in zip(self, other) )
        elif isinstance(other, (int, float)):
            added = tuple( a + other for a in self )
        else:
            raise ValueError("Addition with type {} not supported".format(type(other)))
        
        return self.__class__(*added)

    def __radd__(self, other):
        """ Called if 4 + self for instance """
        return self.__add__(other)
    
    def __sub__(self, other):
        """ Returns the Vec3_t difference of self and other """
        if isinstance(other, Vec3_t):
            subbed = tuple( a - b for a, b in zip(self, other) )
        elif isinstance(other, (int, float)):
            subbed = tuple( a - other for a in self )
        else:
            raise ValueError("Subtraction with type {} not supported".format(type(other)))
        
        return self.__class__(*subbed)
    
    def __rsub__(self, other):
        """ Called if 4 - self for instance """
        return self.__sub__(other)
    
    def __iter__(self):
        return self.values.__iter__()
    
    def __len__(self):
        return len(self.values)
    
    def __getitem__(self, key):
        return self.values[key]
        
    def __repr__(self):
        return str(self.values)
        