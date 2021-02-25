#From operator overloading
#https://gist.github.com/gmalmquist/6b0ee1c38d7263f6158f

class Vec3_t(object):

  @classmethod
  def Lerp(cls, a, b, s):
    return (1.0-s)*a + s*b

  @classmethod
  def Cross(cls, a, b):
    return a^b

  def __init__(self, *components):
    if len(components) == 1:
      self.components = list(components[0])
    else:
      self.components = list(map(float, components))
    if len(self.components) < 3:
      self.components.extend((0,)*3)
      
  @property
  def x(self):
    return self.components[0]

  @x.setter
  def x(self, value):
    self.components[0] = value

  @property
  def y(self):
    return self.components[1]

  @y.setter
  def y(self, value):
    self.components[1] = value

  @property
  def z(self):
    return self.components[2]

  @z.setter
  def z(self, value):
    self.components[2] = value

  def __getitem__(self, index):
    return self.components[index]
  
  def __call__(self, index):
    return self.components[index]
    
  def __setitem__(self, index, value):
    self.components[index] = value

  def __len__(self):
    return len(self.components)

  def __iter__(self):
    return iter(self.components)

  def __add__(self, other):
    return Vector(*(a+b for a,b in zip(self, other)))

  def __mul__(self, other):
    try:
      return sum(a*b for a,b in zip(self, other))
    except:
      return Vector(*(a*other for a in self))

  def __rmul__(self, other):
    return self.__mul__(other)

  def __radd__(self, other):
    return self.__add__(other)

  def __sub__(self, other):
    return Vector(*(a-b for a,b in zip(self, other)))

  def __rsub__(self, other):
    return other + (-self)

  def __neg__(self, other):
    return Vector(*(-a for a in self))

  def __str__(self):
    return '<{}>'.format(', '.join(map(str, self)))

  def __eq__(self, other):
    return tuple(self) == tuple(other)

  def __ne__(self, other):
    return not self.__eq__(other)

  def __hash__(self):
    return hash(tuple(self))

  def __repr__(self):
    return str(self)

  def __xor__(a, b):
    return Vector(a.y*b.z - a.z*b.y,
                  a.z*b.x - a.x*b.z,
                  a.x*b.y - a.y*b.x)

  #@property
  def mag2(self):
    return float(self*self)

  #@property
  def norm(self):
    return float(self*self)**0.5

  @property
  def normalized(self):
    mag2 = self*self
    if mag2 == 0:
      return Vec3_t()
    return Vec3_t(*self) / mag2**0.5
    
#def norm(Vec3_t):