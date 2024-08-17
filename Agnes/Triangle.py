# --------------------------------------------------------------------------
# 
#
#
#
# Katharine Long, Sep 2020
# For Math 5344
# --------------------------------------------------------------------------

import numpy as np
from .QuadratureRule import GaussRule


# --------------------------------------------------------------------------
#
# Class representation of a triangle, which computes the Jacobian and related
# quantities needed in element integrations.
#
class Triangle:
  # Construct with the coordinate pairs of the three vertices
  def __init__(self, A, B, C):
    # Copy points into numpy vectors so we can add/subtract them
    self.A = np.array(A)
    self.B = np.array(B)
    self.C = np.array(C)

    # Now we're going to compute |det(J)| and J^{-T}. We'll need these
    # for element integrations.

    # Form the transpose of the Jacobian of the transformation from the
    # reference right triangle to this triangle
    self.Jt = np.array([self.B - self.A, self.C - self.A])
    # Compute det(J)
    self.detJ = self.Jt[0, 0] * self.Jt[1, 1] - self.Jt[0, 1] * self.Jt[1, 0]
    # The area of the physical triangle is |detJ| times the area of
    # the reference triangle
    self.area = 0.5 * np.abs(self.detJ)
    # Go ahead and form the inverse of J^T. With a 2x2 matrix this is as
    # fast and accurate as an LU factorization followed by triangular solves
    self.JtInv = np.array([[self.Jt[1, 1], -self.Jt[0, 1]],
                           [-self.Jt[1, 0], self.Jt[0, 0]]]) / self.detJ

  # Map quadrature points to physical coordinates
  def ref_to_phys(self, quadX):

    rtn = np.zeros_like(quadX)
    Jt = self.Jt

    for i, x_ref in enumerate(quadX):
      rtn[i, :] = self.A + self.Jt.transpose() @ x_ref

    return rtn
  
  def integrate(self, func, quad=GaussRule(2)):

    if isinstance(func, (float, int, np.double)):
      return func * self.area
    
    elif callable(func):

      xy = self.ref_to_phys(quad.X())
      x = xy[:,0]
      y = xy[:,1]

      fVals = func(x,y)
      return self.area * np.dot(fVals, quad.W())
    
    else:
      raise TypeError('argument {} is neither numeric nor callable'\
                      .format(func)
                      )
    

if __name__ == '__main__':

  from Quadrature import GaussRule

  T = Triangle((1, 1), (4, 2), (0, 4))
  print('Jt=\n', T.Jt)
  print('det(J)=', T.detJ)

  A = np.matmul(T.Jt, T.JtInv)
  I = np.array([[1,0],[0,1]])
  err = np.linalg.norm(A-I)
  print('Jacobian inversion error: {:12.5g}'.format(err))





  def testFunc(n, x, y):
    return (x+y)**n

  def integrationResult(n):
    '''
    Exact results computed with Mathematica,
    Integrate[(x + y)^n, {x, y} \\[Element]Triangle[{{1, 1}, {4, 2}, {0, 4}}]]
    '''
    return (5*2**n*(1 - 2**(3 + n) + 3**(2 + n))) / ((1 + n)*(2 + n))
  
  for n in range(1, 6):
    quad = GaussRule(n)
    x_phys = T.ref_to_phys(quad.X())

    def f(x,y):
      return testFunc(n,x,y)
    I = T.integrate(f, quad=quad)
    ans = integrationResult(n)
    err_a = np.abs(ans-I)
    err_r = err_a/np.abs(ans)

    print('Exact={:20.15g}, quad={:20.15g}, '
          'err_a={:10.5g}, err_r={:10.5g}'.format(ans, I, err_a, err_r))
