from abc import ABC, abstractmethod
from .QuadratureRule import QuadratureRule
from .Triangle import Triangle
import numpy as np


class OneForm(ABC):
  '''
  Abstract interface for one forms
  '''

  def __init__(self, testID : int = 0):
    '''
    Constructor for abstract one-form. The only argument is the ID for
    the test functions appearing in this one-form. 
    '''
    if testID < 0:
      raise ValueError('Test function index should be non-negative')
    
    self._testID = testID

  def testID(self):
    '''
    Return the ID of the test function appearing in this
    one-form.
    '''
    return self._testID

  @abstractmethod
  def localVec(self, tri : Triangle, nodes : tuple):
    '''
    Interface for computing a local vector on a triangle. 
    '''
    pass


class QuadratureOneForm(OneForm):
  '''
  One form to be computed by quadrature 
  '''

  def __init__(self, quad : QuadratureRule, testID:int=0):
    super().__init__(testID = testID)
    self._quad = quad
    self._x = quad.X()[:,0].reshape((quad.n(),1))
    self._y = quad.X()[:,1].reshape((quad.n(),1))

    self._phiAtQuadPts = np.hstack((1 - self._x - self._y,
                                   self._x, self._y))
  
  def x(self):
    return self._x
  
  def y(self):
    return self._y
  
  def w(self):
    return self._quad.W()
  
  def xy(self):
    return self._quad.X()
  
  def phiAtQuadPts(self):
    return self._phiAtQuadPts
  

class ConstCoeffOneForm(OneForm):

  def __init__(self, coeff : float, testID:int=0):
    super().__init__(testID = testID)
    self._coeff = coeff
    self._vec = np.array([1.0,1.0,1.0])/6.0

  
  def localVec(self, tri : Triangle, nodes : tuple):
    return self._coeff * tri.detJ * self._vec
  

class VarCoeffOneForm(QuadratureOneForm):

  def __init__(self, quad, coeffFunc, dfs = None, testID:int=0):
    super().__init__(quad, testID = testID)
    self._coeffFunc = coeffFunc
    self._dfs = dfs
    self._phiWork = np.zeros_like(self.phiAtQuadPts())
    self._valsWork = np.zeros_like(self.w())

  
  def localVec(self, tri : Triangle, nodes : tuple):

    if self._dfs == None:
      dfVals = None
    else:
      dfVals = self._dfs.interpolate(nodes, self._phiAtQuadPts)

    XY = tri.ref_to_phys(self.xy())

    vals = self._coeffFunc(XY[:,0], XY[:,1], dfVals)

    # Compute w * vals elementwise, writing result into workspace
    np.multiply(vals, self.w(), self._valsWork)

    # Sum over quad points
    # TODO: figure out how to do this with numpy multiply
    for i,(v,f) in enumerate(zip(self._valsWork, self._phiAtQuadPts)):
      self._phiWork[i,:] = v*f

    sum = np.sum(self._phiWork, 0)

    return tri.area * sum
  


  


  


