from abc import ABC, abstractmethod
from .QuadratureRule import QuadratureRule
from .Triangle import Triangle
import numpy as np


class TwoForm(ABC):
  '''
  Abstract interface for two forms
  '''

  def __init__(self, testID:int=0, unkID:int=0):
    '''
    Constructor for abstract two-form. The only arguments are the ID for
    the test and unknown functions 
    '''
    if testID < 0:
      raise ValueError('Test function index should be non-negative')
    if unkID < 0:
      raise ValueError('Unknown function index should be non-negative')
    
    self._testID = testID
    self._unkID = unkID
    
  @abstractmethod
  def localMat(self, tri : Triangle, nodes : tuple):
    '''
    Interface for computing a local matrix on a triangle. 
    '''
    pass

  def testID(self):
    '''
    Return the ID of the test function appearing in this two-form.
    '''
    return self._testID

  def unkID(self):
    '''
    Return the ID of the unknown function appearing in this two-form.
    '''
    return self._unkID


class LaplacianTwoForm(TwoForm):
  '''
  Local Laplacian matrix builder
  '''

  def __init__(self, coeff=1.0, testID:int=0, unkID:int=0):
    super().__init__(testID=testID, unkID=unkID)
    if not isinstance(coeff, (int, float, np.double)):
      raise TypeError('Coefficient argument to LaplacianTwoForm '
                      'should be a constant; argument was {}'\
                      .format(coeff))
    
    self._coeff = coeff
  
  @staticmethod
  def gradPhiRef():
    '''
    The gradients of the three Lagrange basis functions on the reference
    triangle are constant vectors and can be computed once and for all.
    These are packed columnwise into the 2 by 3 matrix:
    [ -1 1 0 ]
    [ -1 0 1 ]
    '''
    return np.array([[-1, 1, 0], [-1, 0, 1]])

  
  def localMat(self, tri:Triangle, nodes:tuple):
    '''
    Form the local matrix. Since the gradients are constant, no quadrature
    is needed, but we do have to transform the gradients from reference to
    physical coordinates. That's done as follows:
    grad_phys(phi) = J^{-T} * grad_ref(phi).
    This can be done on all three basis functions at once since we've packed
    their gradients into columns of the matrix LocalLaplacian.gradPhiRef.
    '''
    # Transform the gradients
    gradPhi = np.matmul(tri.JtInv, LaplacianTwoForm.gradPhiRef())
    # Form the local matrix
    rtn = self._coeff*tri.area*np.matmul(np.transpose(gradPhi), gradPhi)
    return rtn
  

# --------------------------------------------------------------------------
#
# The LocalMass class provides a localMat() function for computing
# the local mass matrix (aka Gram matrix)
# \int_T phi_i*phi_j
# on a triangle with degree one Lagrange basis functions. In transforming
# from reference to physical triangle, the integral is simply scaled by the
# Jacobian determinant. We can compute all the integrals once and for all
# offline, and re-use the matrix for all triangles.
#
class MassTwoForm(TwoForm):

  def __init__(self, coeff=1.0, testID:int=0, unkID:int=0):
    super().__init__(testID=testID, unkID=unkID)
    if not isinstance(coeff, (int, float, np.double)):
      raise TypeError('Coefficient argument to LaplacianTwoForm '
                      'should be a constant; argument was {}'\
                      .format(coeff))
    
    self._coeff = coeff
  
  # Here's the local mass matrix, up to scaling by the Jacobian, computed
  # offline using Mathematica. Since it's the same for all triangles,
  # make it a class variable rather than an instance variable.
  M = np.array([[2, 1, 1], [1, 2, 1], [1, 1, 2]]) / 12.0

  # Compute the local mass matrix on a triangle.
  def localMat(self, tri : Triangle, nodes : tuple):
    return tri.area * MassTwoForm.M

