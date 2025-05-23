import numpy as np
import pprint

class QuadratureRule:
  def __init__(self, wQuad, xQuad, order=4, dim=2, name='Quad'):
    self._order = order
    self._xQuad = xQuad
    self._wQuad = wQuad
    self._dim = dim
    self._name = name

    assert(len(xQuad) == len(wQuad))

  def X(self):
    return self._xQuad
  
  def W(self):
    return self._wQuad
  
  def nPts(self):
    return len(self._wQuad)
  
  def order(self):
    return self._order
  
  def dim(self):
    return self._dim
  
  
  def show(self):
    pprint.pp('Quadrature rule: {}'.format(self._name))
    for x,w in zip(self.X(), self.W()):
      print('x={} w={}'.format(x,w)) 

  def test(self, tol=1.0e-14):

    dim = self.dim()    
    if dim==1:
      err = self._test1D()
    elif dim==2:
      err = self._test2D()
    
    passed = err <= tol

    if not passed:
      print('Test FAILED! Showing quadrature rule:')
      self.show()
    
    return passed, err

  def _test1D(self):
    
    import numpy.polynomial.polynomial as poly	
	
    a = np.ones(self.order()+1)
    P = poly.Polynomial(a)

    answer = P.integ()(1) - P.integ()(-1)

    sum = 0
    for w, x in zip(self.W(), self.X()):
      sum += w * P(x)

    err = np.abs(answer - sum)
    return err
              
  def _test2D(self):
    p = self.order()
    answer = _test_poly2D_exact_integral(p)

    sum = 0
    for w, xy in zip(self.W(), self.X()):
      sum += w * _test_poly2D(p,xy)
    sum = 0.5*sum

    err = np.abs(answer - sum)
    return err

    


  
  
def _test_poly2D_exact_integral(p : int):
  
  sum = 0

  for i in range(0, p+1):
    for j in range(0, p + 1 - i):
      sum += _fact(i)*_fact(j)/_fact(i+j+2)
  
  return sum

def _fact(n : int):
  
  assert(n>=0)

  if n==0:
    return 1
  return n*_fact(n-1)
  
    

    
def _test_poly2D(p : int, xy):
  sum = 0
  x = xy[0]
  y = xy[1]
  
  powx = 1
  for i in range(0, p+1):
    powy = 1
    for j in range(0, p + 1 - i):
      sum += powx * powy
      powy *= y
    powx *= x

  return sum




def main():
  maxP = 9

  for p in range(1, maxP+1):
    xy=(1/3,1/3)
    pVal = _test_poly2D(p, xy)
    Q = _test_poly2D_exact_integral(p)
    print('p={} P(1/3,1/3)={}, Q={}'.format(p,pVal,Q))

    



    
if __name__=='__main__':

  main()